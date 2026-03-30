#include "antispoofing/AntiSpoofing.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_cdf.h>

// =============================================================================
// M1: 星历规范性检测 —— 对应 MATLAB chnEphemNormativity + checkGPSEph + checkBDSEph
// =============================================================================
void AntiSpoofing::detectMethod01(const EphemerisData& eph, SatSpoofResult& res)
{
    bool invalid = false;
    const double PI = M_PI;

    // ============================== GPS 星历检查 ==============================
    if (eph.sys == 'G')
    {
        bool Toc_TOW_result    = (eph.Toc_TOW < 0 || eph.Toc_TOW > 604784);
        bool Toc_WN_result     = (eph.Toc_WN < 0 || eph.Toc_WN > 8191);
        bool PRN_result        = (eph.PRN < 1 || eph.PRN > 32);
        bool WN_result         = (eph.WN < 0 || eph.WN > 8191);
        bool Health_result     = (eph.Health != 0);
        bool IODE_result       = (eph.IODE < 0 || eph.IODE > 255);

        bool af0_result        = (fabs(eph.af0) > 9.7656e-04);
        bool af1_result        = (fabs(eph.af1) > 3.7253e-09);
        bool af2_result        = (fabs(eph.af2) > 3.5527e-15);
        bool TGD_result        = (fabs(eph.TGD) > 5.9605e-08);

        bool sqrt_a_result     = (eph.sqrt_a < 2530 || eph.sqrt_a > 8192);
        bool e_result          = (eph.e < 0 || eph.e > 0.03);
        bool Toe_result        = (eph.Toe < 0 || eph.Toe > 604784);
        bool Toes_result       = (eph.Toes < 0 || eph.Toes > 604784);
        bool M0_result         = (fabs(eph.M_0) > 2 * PI);
        bool Omega0_result     = (fabs(eph.Omega_0) > 2 * PI);
        bool omega_result      = (fabs(eph.omega) > 2 * PI);
        bool i0_result         = (fabs(eph.i_0) > 2 * PI);

        bool Crs_result        = (fabs(eph.Crs) > 1024);
        bool Crc_result        = (fabs(eph.Crc) > 1024);
        bool Cuc_result        = (fabs(eph.Cuc) > 6.2e-5);
        bool Cus_result        = (fabs(eph.Cus) > 6.2e-5);
        bool Cic_result        = (fabs(eph.Cic) > 6.2e-5);
        bool Cis_result        = (fabs(eph.Cis) > 6.2e-5);

        bool Delta_n_result    = (fabs(eph.Delta_n) > 1.4901e-08);   // 2^-28 * PI
        bool Omega_dot_result  = (eph.Omega_dot < -2.0e-6 || eph.Omega_dot > 1e-9);
        bool I_dot_result      = (fabs(eph.I_dot) > 3e-9);

        // 合并所有异常（任意一项异常 → 星历非法）
        invalid = Toc_TOW_result || Toc_WN_result || PRN_result || WN_result || Health_result || IODE_result
               || af0_result || af1_result || af2_result || TGD_result
               || sqrt_a_result || e_result || Toe_result || Toes_result
               || M0_result || Omega0_result || omega_result || i0_result
               || Crs_result || Crc_result || Cuc_result || Cus_result || Cic_result || Cis_result
               || Delta_n_result || Omega_dot_result || I_dot_result;
    }

    // ============================== BDS 星历检查 ==============================
    else if (eph.sys == 'C')
    {
        int useGPS = cfg_.anti_spoofing_.chnEphemNormativity.useGPS;

        bool Toc_TOW_result    = (eph.Toc_TOW < 0 || eph.Toc_TOW > 604500);
        bool Toc_WN_result     = (eph.Toc_WN < 0 || eph.Toc_WN > 8191);
        bool PRN_result        = (eph.PRN < (1 + 32*useGPS) || eph.PRN > (30 + 32*useGPS));
        bool WN_result         = (eph.WN < 0 || eph.WN > 8191);
        bool Health_result     = (eph.Health != 0);
        bool IODE_result       = (eph.IODE < 0 || eph.IODE > 255);

        bool af0_result        = (fabs(eph.af0) > 9.7656e-04);
        bool af1_result        = (fabs(eph.af1) > 1.8626e-09);
        bool af2_result        = (fabs(eph.af2) > 1.3878e-17);
        bool TGD_result        = (fabs(eph.TGD) > 1.2e-7);

        bool sqrt_a_result     = (eph.sqrt_a < 5000 || eph.sqrt_a > 7000);
        bool e_result          = (eph.e < 0 || eph.e > 0.05);
        bool Toe_result        = (eph.Toe < 0 || eph.Toe > 604500);
        bool Toes_result       = (eph.Toes < 0 || eph.Toes > 604500);
        bool M0_result         = (fabs(eph.M_0) > 2 * PI);
        bool Omega0_result     = (fabs(eph.Omega_0) > 2 * PI);
        bool omega_result      = (fabs(eph.omega) > 2 * PI);
        bool i0_result         = (fabs(eph.i_0) > 2 * PI);

        bool Crs_result        = (fabs(eph.Crs) > 32768);
        bool Crc_result        = (fabs(eph.Crc) > 32768);
        bool Cuc_result        = (fabs(eph.Cuc) > 1e-3);
        bool Cus_result        = (fabs(eph.Cus) > 1e-3);
        bool Cic_result        = (fabs(eph.Cic) > 1e-3);
        bool Cis_result        = (fabs(eph.Cis) > 1e-3);

        bool Delta_n_result    = (fabs(eph.Delta_n) > 1.2e-8);
        bool Omega_dot_result  = (eph.Omega_dot < -3.0e-6 || eph.Omega_dot > 5e-8);
        bool I_dot_result      = (fabs(eph.I_dot) > 3e-9);

        invalid = Toc_TOW_result || Toc_WN_result || PRN_result || WN_result || Health_result || IODE_result
               || af0_result || af1_result || af2_result || TGD_result
               || sqrt_a_result || e_result || Toe_result || Toes_result
               || M0_result || Omega0_result || omega_result || i0_result
               || Crs_result || Crc_result || Cuc_result || Cus_result || Cic_result || Cis_result
               || Delta_n_result || Omega_dot_result || I_dot_result;
    }

    // ===================== 输出结果 =====================
    res.is_spoof[0] = invalid;  // 异常 = 欺骗
    res.score[0]    = invalid ? 1.0 : 0.0;
}

// ========================== M2 卫星可见性 ==========================
void AntiSpoofing::detectMethod02(const Eigen::Vector3d& ps, SatSpoofResult& res) {
    auto& pos = cfg_.anti_spoofing_.chnSatVisibility.rec_pos_lla;
    Eigen::Vector3d rec_lla(pos[0], pos[1], pos[2]);

    if (std::isnan(rec_lla(0))) {
        res.is_spoof[1] = false; 
        res.score[1] = 0; 
        return;
    }
    auto [az, el] = cal_azel(rec_lla, ps);
    double thr = cfg_.anti_spoofing_.chnSatVisibility.Thr;
    res.is_spoof[1] = (el < thr);
    res.score[1] = res.is_spoof[1] ? 1 : 0;
}

// =============================================================================
// M3: 载噪比(CN0)异常监测 — 稳健版
// 逻辑：取最新的 最多N个 有效数据 → 平均 → 与阈值比较
// 不足N个 → 有多少用多少，立刻判决
// =============================================================================
void AntiSpoofing::detectMethod03(const SatHistory& hist, SatSpoofResult& res)
{
    // 1. 从配置读取参数
    int window_size   = cfg_.anti_spoofing_.chnCn0.WindowSize;    // 滑窗长度 N
    double cn0_thresh = cfg_.anti_spoofing_.chnCn0.Thr;           // 载噪比阈值

    // 2. 倒序遍历：取【最新】的有效数据（最多N个）
    double sum = 0.0;
    int cnt = 0;

    // 从最新时刻往旧时刻遍历
    for (int i = hist.size() - 1; i >= 0 && cnt < window_size; --i) {
        if (hist.effective[i] == 1) {
            sum += hist.cn0[i];
            cnt++;
        }
    }

    // 3. 无任何有效数据 → 不判决
    if (cnt == 0) {
        res.is_spoof[2] = false;
        res.score[2] = 0.0;
        return;
    }

    // 4. 计算平均（有多少算多少）
    double avg = sum / cnt;

    // 5. 判决：超过阈值=欺骗
    res.is_spoof[2] = (avg > cn0_thresh);
    res.score[2] = avg;
}

// =============================================================================
// M4: 多普勒一致性检测 (Doppler Consistency)
// 规则：
// 1. 有效点数 < 2 → 不检测
// 2. 不足窗长 → 有多少用多少
// 3. 只使用【最新】的有效数据
// 4. 差分多普勒 = 载波多普勒 - 1540 * 码多普勒
// 5. 窗内平均后绝对值 > 阈值 → 判定欺骗
// =============================================================================
void AntiSpoofing::detectMethod04(const SatHistory& hist, SatSpoofResult& res)
{
    // 1. 从配置读取参数
    int window_size = cfg_.anti_spoofing_.chnDoppCon.WindowSize;    // 窗长
    double thr_dopp = cfg_.anti_spoofing_.chnDoppCon.Thr;    // 阈值
    double Fif      = cfg_.anti_spoofing_.chnDoppCon.Fif;    // 中频

    // 2. 倒序提取【最新N个有效数据】(时间、伪距、载波相位)
    std::vector<double> t_vec, rho_vec, acph_vec;
    for (int i = hist.size() - 1; i >= 0; i--)
    {
        if (hist.effective[i] == 1)
        {
            t_vec.push_back(hist.time[i]);
            rho_vec.push_back(hist.rho[i]);
            acph_vec.push_back(hist.acph[i]);

            // 够窗长就停止
            if ((int)t_vec.size() >= window_size)
                break;
        }
    }

    // 3. 有效数据 < 2 → 无法差分，直接返回
    if (t_vec.size() < 2)
    {
        res.is_spoof[3] = false;
        res.score[3] = 0.0;
        return;
    }

    // 4. 反转 → 变成【旧→新】时间顺序（方便差分）
    std::reverse(t_vec.begin(), t_vec.end());
    std::reverse(rho_vec.begin(), rho_vec.end());
    std::reverse(acph_vec.begin(), acph_vec.end());

    int n = (int)t_vec.size();
    std::vector<double> diff_dopp;

    // 5. 逐历元计算差分多普勒（和MATLAB完全一致）
    for (int i = 1; i < n; i++)
    {
        double dt = t_vec[i] - t_vec[i-1];
        if (dt < 1e-6)
            continue;

        // ===================== 码相位 / 码多普勒 =====================
        double code_phase_curr = (rho_vec[i]  / CLIGHT) * 1023.0e3 + (t_vec[i] - t_vec[0]) * 1023.0e3;
        double code_phase_prev = (rho_vec[i-1]/CLIGHT) * 1023.0e3 + (t_vec[i-1] - t_vec[0]) * 1023.0e3;
        double code_dopp = (code_phase_curr - code_phase_prev) / dt - 1023.0e3;

        // ===================== 载波相位 / 载波多普勒 =====================
        double carr_dopp = (acph_vec[i] - acph_vec[i-1]) / dt - Fif * 1.0e3;

        // ===================== 差分多普勒 =====================
        double dd = carr_dopp - 1540.0 * code_dopp;
        diff_dopp.push_back(dd);
    }

    // 6. 计算平均差分多普勒
    double sum = 0.0;
    for (double v : diff_dopp)
        sum += v;
    double avg_diff = sum / (double)diff_dopp.size();
    double score = fabs(avg_diff);

    // 7. 阈值判决
    res.is_spoof[3] = (score > thr_dopp);
    res.score[3] = score;
}

// =============================================================================
// M5: 单星钟漂变化监测 (Clock Drift Monitor) —— 100% 匹配 MATLAB
// 修正点：dt_valid = 两个有效数据之间相隔的【原始历元数】
// =============================================================================
void AntiSpoofing::detectMethod05(const SatHistory& hist, SatSpoofResult& res)
{
    // 1. 从配置读取参数
    int window_size = cfg_.anti_spoofing_.chnClkDft.WindowSize;    // 窗长
    double thr_clk  = cfg_.anti_spoofing_.chnClkDft.Thr;           // 阈值
    double T_samp   = cfg_.frequency_.T_samp;  // 接收机采样周期 1s

    // ===================== 1. 收集【最新N个有效点】+ 记录它们在队列中的原始索引 =====================
    std::vector<int> valid_indices;    // 存储有效点在 SatHistory 中的索引
    std::vector<double> valid_times;   // 存储有效点的时间

    // 倒序取最新N个有效点
    for (int i = hist.size() - 1; i >= 0 && (int)valid_indices.size() < window_size; --i) {
        if (hist.effective[i] == 1) {
            valid_indices.push_back(i);
            valid_times.push_back(hist.time[i]);
        }
    }

    // 不足2个有效点 → 不计算
    if (valid_indices.size() < 2) {
        res.is_spoof[4] = false;
        res.score[4] = 0.0;
        return;
    }

    // 反转 → 时间从旧到新
    std::reverse(valid_indices.begin(), valid_indices.end());
    std::reverse(valid_times.begin(), valid_times.end());

    // ===================== 2. 按 MATLAB 公式计算钟漂 =====================
    double max_clk_dft = 0.0;
    bool is_spoof = false;

    for (int i = 1; i < (int)valid_indices.size(); ++i) {
        // 相邻两个有效点的【观测时间差】
        double dt = valid_times[i] - valid_times[i-1];

        // 相邻两个有效点在【历史队列中的索引差】→ 这就是你要的 dt_valid！
        int dt_valid = valid_indices[i] - valid_indices[i-1];

        // MATLAB 公式：ClkDft = (dt - T_sampling) / dt_valid / T_sampling
        double clk_dft = (dt - T_samp) / dt_valid / T_samp;
        double abs_val = fabs(clk_dft);

        // 记录最大值
        if (abs_val > max_clk_dft)
            max_clk_dft = abs_val;

        // 任意一个超限 → 判骗
        if (abs_val > thr_clk)
            is_spoof = true;
    }

    // 输出结果
    res.is_spoof[4] = is_spoof;
    res.score[4] = max_clk_dft;
}

// =============================================================================
// M6: 发射时间 / 伪距合理性检测 (Pseudorange Reasonableness)
// 逻辑：
// 1. 取最新 window_size 个有效数据
// 2. 不足 window_size → 有多少用多少
// 3. 计算伪距是否合理（与阈值比较）
// 4. 任意一个不合理 → 判定欺骗
// =============================================================================
void AntiSpoofing::detectMethod06(const SatHistory& hist, SatSpoofResult& res)
{
    // 1. 配置参数
    int window_size     = cfg_.anti_spoofing_.chnRho.WindowSize;    // 窗口大小
    auto thr            = cfg_.anti_spoofing_.chnRho.Thr;
    double rho_min_thr  = thr[0] * CLIGHT;   // 伪距最小值
    double rho_max_thr  = thr[1] * CLIGHT;   // 伪距最大值

    // 2. 倒序取【最新N个有效伪距】
    std::vector<double> valid_rho;
    for (int i = hist.size() - 1; i >= 0 && (int)valid_rho.size() < window_size; --i)
    {
        if (hist.effective[i] == 1)
        {
            valid_rho.push_back(hist.rho[i]);
        }
    }

    // 3. 无有效数据 → 不判决
    if (valid_rho.empty())
    {
        res.is_spoof[5] = false;
        res.score[5] = 0.0;
        return;
    }

    // 4. 检查合理性：小于最小值 或 大于最大值 → 不合理
    bool is_spoof = false;
    double max_abs_err = 0.0;

    for (double rho : valid_rho)
    {
        if (rho < rho_min_thr || rho > rho_max_thr)
        {
            is_spoof = true;
        }

        // 记录偏离程度（用于score）
        double err = 0.0;
        if (rho < rho_min_thr) err = rho_min_thr - rho;
        if (rho > rho_max_thr) err = rho - rho_max_thr;
        if (err > max_abs_err) max_abs_err = err;
    }

    // 5. 输出结果
    res.is_spoof[5] = is_spoof;
    res.score[5] = max_abs_err;
}
