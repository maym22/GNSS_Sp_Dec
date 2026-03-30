#include "antispoofing/AntiSpoofing.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_cdf.h>
#include <unordered_set>

// ============================================================================
// M7：载噪比（CN0）相关性欺骗检测（统计显著性检验版）
// 逻辑：当前有效卫星两两配对 → 共用有效时刻计算相关系数
//      → 基于Pfa和实际样本长度n计算显著性阈值 → 超过则判欺骗
// ============================================================================
void AntiSpoofing::detectMethod07(std::vector<SatSpoofResult>& sat_results)
{
    if (!method_enabled_[7]) 
        return;

    // 从配置读取参数（现在只需要 Pfa，不再需要固定阈值）
    int windowSize     = cfg_.anti_spoofing_.mchnCn0Corr.WindowSize;
    int minWindowSize  = cfg_.anti_spoofing_.mchnCn0Corr.MinWindowSize;
    double alpha       = cfg_.anti_spoofing_.mchnCn0Corr.Pfa; // 虚警率

    // ======================
    // 1. 获取【当前有效卫星】
    // ======================
    std::unordered_map<int, const SatHistory*> curr_valid_sat;

    for (const auto& sr : sat_results) {
        int prn = sr.prn;
        auto it = sat_history_.find(prn);
        if (it == sat_history_.end()) continue;

        const SatHistory& hist = it->second;
        if (hist.size() == 0) continue;

        // 最新时刻有效（effective ==1）才保留
        int lastIdx = hist.size() - 1;
        if (hist.effective[lastIdx] == 1) {
            curr_valid_sat[prn] = &hist;
        }
    }

    // 卫星列表
    std::vector<int> prnList;
    for (auto& p : curr_valid_sat) {
        prnList.push_back(p.first);
    }

    // 至少需要2颗卫星才能计算相关性
    if (prnList.size() < 2) {
        return;
    }

    // ======================
    // 2. 两两卫星计算 CN0 相关系数 + 统计显著性检验
    // ======================
    std::unordered_set<int> spoofedPrns;

    for (int i = 0; i < prnList.size(); ++i) {
        int prn_i = prnList[i];
        const SatHistory* h_i = curr_valid_sat[prn_i];

        for (int j = i + 1; j < prnList.size(); ++j) {
            int prn_j = prnList[j];
            const SatHistory* h_j = curr_valid_sat[prn_j];

            // ==============================
            // 【最终正确版】从最新时刻对齐取共同有效CN0
            // ==============================
            std::vector<double> xi, xj;

            int idx_i = h_i->size() - 1;  // 卫星i最后一位
            int idx_j = h_j->size() - 1;  // 卫星j最后一位

            // 从最新时刻向前一一匹配，直到一方到头 或 取够windowSize
            while (idx_i >= 0 && idx_j >= 0 && xi.size() < windowSize)
            {
                // 两者同时有效 → 加入序列
                if (h_i->effective[idx_i] == 1 && h_j->effective[idx_j] == 1)
                {
                    xi.push_back(h_i->cn0[idx_i]);
                    xj.push_back(h_j->cn0[idx_j]);
                }

                // 同时往前挪一位（时间对齐）
                idx_i--;
                idx_j--;
            }

            int n = xi.size();
            if (n < minWindowSize) {
                continue;
            }

            // ==============================
            // 计算皮尔逊相关系数 PCC
            // ==============================
            Eigen::Map<Eigen::VectorXd> ei(xi.data(), n);
            Eigen::Map<Eigen::VectorXd> ej(xj.data(), n);

            double mean_i = ei.mean();
            double mean_j = ej.mean();

            Eigen::VectorXd di = ei.array() - mean_i;
            Eigen::VectorXd dj = ej.array() - mean_j;

            double numerator   = di.dot(dj);
            double var_i       = di.squaredNorm();
            double var_j       = dj.squaredNorm();
            double pcc         = numerator / (std::sqrt(var_i * var_j) + 1e-9);

            // ==============================
            // 【关键修改】根据实际 n 和 Pfa 计算显著性阈值
            // ==============================
            double nu        = (double)(n - 2); // 阈值随着n的降低和Pfa的减小而升高
            double t_alpha   = gsl_cdf_tdist_Pinv(1.0 - alpha, nu);
            double gamma_thr = t_alpha / std::sqrt(nu + t_alpha * t_alpha);

            // ==============================
            // 超过统计阈值 → 判欺骗
            // ==============================
            if (pcc > gamma_thr) {
                spoofedPrns.insert(prn_i);
                spoofedPrns.insert(prn_j);
            }
        }
    }

    // ======================
    // 3. 输出结果到 sat_results
    // ======================
    for (auto& sr : sat_results) {
        bool isSpoof = spoofedPrns.count(sr.prn) > 0;
        sr.is_spoof[7] = isSpoof;
        sr.score[7]    = isSpoof ? 1.0 : 0.0;
    }
}


// =============================================================================
// M8: 多星多普勒单差检测（最终工程版）
// 满足所有要求：
// 1. 仅处理【当前历元最新时刻 effective=1】的卫星
// 2. 共同有效点 < 5 → 不跳过，RMSE 设为极大值（不参与线性判决）
// 3. 按窗长取最新点
// 4. 严格按 effective 取公共时刻
// =============================================================================
void AntiSpoofing::detectMethod08(std::vector<SatSpoofResult>& sat_results)
{

    double thr_rmse    = cfg_.anti_spoofing_.mchnDoppSD.Thr;
    int window_size    = cfg_.anti_spoofing_.mchnDoppSD.WindowSize;
    int MIN_PTS        = cfg_.anti_spoofing_.mchnDoppSD.MinWindowSize; // 最小有效点数，不足则认为非线性

    // ==========================
    // 【要求1 严格实现】
    // 只保留：当前历元【最新时刻 effective == 1】的卫星
    // ==========================
    std::unordered_map<int, const SatHistory*> curr_valid_sat;
    // 考虑一下这段话是否必要
    for (const auto& sr : sat_results) {
        int prn = sr.prn;
        auto it = sat_history_.find(prn);
        if (it == sat_history_.end()) continue;

        const SatHistory& hist = it->second;
        if (hist.size() == 0) continue;

        // 最新时刻 == 最后一位
        int latest_idx = hist.size() - 1;
        if (hist.effective[latest_idx] == 1) {
            curr_valid_sat[prn] = &hist;
        }
    }

    std::vector<int> prn_list;
    for (auto& p : curr_valid_sat) {
        prn_list.push_back(p.first);
    }

    if (prn_list.size() < 2) {
        return;
    }

    // ==========================
    // 两两卫星计算 DFD
    // ==========================
    std::map<std::pair<int, int>, double> rmse_map;

    for (uint i = 0; i < prn_list.size(); ++i) {
        int p1 = prn_list[i];
        const SatHistory& h1 = *curr_valid_sat[p1];

        for (uint j = i + 1; j < prn_list.size(); ++j) {
            int p2 = prn_list[j];
            const SatHistory& h2 = *curr_valid_sat[p2];

            // 取最新 window_size 长度区间
            int start_t = std::max(0, (int)h1.size() - window_size);
            std::vector<double> dfd_seq;

            // ==========================
            // 公共有效时刻（必须同时有效）
            // ==========================
            for (int t = start_t; t < (int)h1.size(); ++t) {
                bool eff1 = (t < (int)h1.effective.size()) && (h1.effective[t] == 1);
                bool eff2 = (t < (int)h2.effective.size()) && (h2.effective[t] == 1);
                if (eff1 && eff2) {
                    dfd_seq.push_back(h1.dopp[t] - h2.dopp[t]);
                }
            }

            double rmse;

            // ==========================
            // 【要求2 严格实现】
            // 点数不足 → 不跳过！RMSE = 极大值，表示不满足线性
            // ==========================
            if ((int)dfd_seq.size() < MIN_PTS) {
                rmse = 1e9;
            }
            else {
                // 最小二乘拟合 dfd = a * k + b
                int n = (int)dfd_seq.size();
                double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;

                for (int k = 0; k < n; ++k) {
                    double x = k;
                    double y = dfd_seq[k];
                    sum_x += x;
                    sum_y += y;
                    sum_xy += x * y;
                    sum_x2 += x * x;
                }

                double det = n * sum_x2 - sum_x * sum_x;
                if (fabs(det) < 1e-6) {
                    rmse = 1e9;
                } else {
                    double a = (n * sum_xy - sum_x * sum_y) / det;
                    double b = (sum_y * sum_x2 - sum_x * sum_xy) / det;

                    double sse = 0;
                    for (int k = 0; k < n; ++k) {
                        double res = dfd_seq[k] - (a * k + b);
                        sse += res * res;
                    }
                    rmse = sqrt(sse / n);
                }
            }

            rmse_map[{p1, p2}] = rmse;
            rmse_map[{p2, p1}] = rmse;
        }
    }

    // ==========================
    // 投票：有多少颗卫星满足 RMSE < 阈值
    // ==========================
    std::unordered_map<int, int> linear_count;
    for (int prn : prn_list) {
        linear_count[prn] = 0;
    }

    for (auto& entry : rmse_map) {
        int a = entry.first.first;
        int b = entry.first.second;
        double r = entry.second;

        if (r < thr_rmse) {
            linear_count[a]++;
            linear_count[b]++;
        }
    }

    // ==========================
    // 输出结果：≥2 个线性相关卫星 → 判骗
    // ==========================
    for (auto& sr : sat_results) {
        int prn = sr.prn;
        if (linear_count.count(prn)) {
            sr.is_spoof[7] = (linear_count[prn] >= 5);
            sr.score[7] = linear_count[prn];
        } else {
            sr.is_spoof[7] = false;
            sr.score[7] = 0;
        }
    }
}
