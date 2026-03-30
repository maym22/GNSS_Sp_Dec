#include "antispoofing/AntiSpoofing.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_cdf.h>
#include <unordered_set>

// ============================================================================
// M9：伪距域 RAIM 欺骗检测
// 依赖：伪距PVT解算结果 rhor (m)
// 自由度：卫星数 - 4
// ============================================================================
void AntiSpoofing::detectMethod09(const UserPvtResult& userPvt, const int size, std::vector<SatSpoofResult>& sat_results)
{
    bool raimAlarm  = false;
    // 多普勒RAIM需要至少8颗卫星 + 解算有效
    if (size < 5 || !userPvt.isValid) 
        return;

    double Pfa      = cfg_.anti_spoofing_.AntRaim.Pfa;
    double Thr      = gsl_cdf_chisq_Pinv(1 - Pfa, size-4);
    double Sigma    = cfg_.anti_spoofing_.AntRaim.Sigma;

    
    double test = userPvt.rhor.squaredNorm() /(Sigma*Sigma);

    raimAlarm   = test > Thr;

    // 保存结果
    for (auto& sr : sat_results) {
        sr.is_spoof[8] = raimAlarm;
        sr.score[8] = raimAlarm ? test : 0.0;
    }

}

// 工具函数：根据选中的卫星索引，计算残差平方和
// 输入：selectedIdx -> 从 goodIndices 中挑选的卫星索引子集
// 输出：残差平方和 sum(rhor^2)
static UserPvtResult computeResidualSumSquared(UserPvtSolver& solver,
                                        const SatPVTStruct& satpvt,
                                        const std::vector<int>& selectedIdx)
{
    UserPvtResult outResult;

    // 卫星数量不足4颗，无法定位，返回极大值
    if (selectedIdx.size() < 4) {
        return outResult; // 默认残差平方和为0，isValid为false
    }

    // 调用已实现的最小二乘解算
    UserPvtResult res;
    solver.solvePesuPvt(satpvt, selectedIdx, outResult, 0);

    // 计算残差平方和
    return outResult;
}

// ==============================
// M10 生成式RAIM (Generative RAIM)
// 优化版本：全程仅使用 goodIndices 索引操作，无PRN来回转换
// 1. 枚举5颗卫星的最优核集合
// 2. 校验核有效性
// 3. 逐星加入验证，标记欺骗卫星
// ==============================
void AntiSpoofing::detectMethod10(const SatPVTStruct& satpvt, const std::vector<int>& goodIndices, std::vector<SatSpoofResult>& sat_results)
{
    // ==========================================
    // 0. 初始化与参数配置
    // ==========================================
    const double Sigma    = cfg_.anti_spoofing_.AntGRaim.Sigma;   // 噪声方差因子（与M9/M11保持一致）
    const double Pfa      = cfg_.anti_spoofing_.AntGRaim.Pfa;     // 恒虚警率
    const int    Mode     = cfg_.anti_spoofing_.AntGRaim.Mode;   // 模型选择：1=全枚举，2=随机选取1个核集合
    const int KERNEL_SIZE = 5;         // 核集合固定为5颗星

    int nAvail = goodIndices.size();
    // 可用卫星不足5颗，算法无法运行
    if (nAvail < KERNEL_SIZE) {
        std::cout << "Not enough available satellites for G-RAIM." << std::endl;
        return;
    }

    // ==========================================
    // Step 1：寻找最优5星核集合
    // 枚举所有5星组合，选择残差平方和最小的一组作为基准核
    // ==========================================
    std::vector<int> bestKernelIdx;
    UserPvtSolver solver;
    UserPvtResult kernelPvt, bestKernelPvt;
    double currRes;
    double minResidual = 1e18;

    // 生成所有 5 颗星的索引组合（基于 goodIndices）
    std::vector<bool> combMask(nAvail, false);
    std::fill(combMask.end() - KERNEL_SIZE, combMask.end(), true);
    // std::fill(combMask.begin(), combMask.begin() + KERNEL_SIZE, true);
    
    // 核集合门限：5星定位，自由度 = 5-4 = 1
    int    dofKernel = 1;
    double thrKernel = gsl_cdf_chisq_Pinv(1.0 - Pfa, dofKernel);

    do {
        std::vector<int> currIdx;
        for (int i = 0; i < nAvail; ++i) {
            if (combMask[i]) {
                currIdx.push_back(goodIndices[i]);
            }
        }

        // 计算当前组合残差平方和
        solver.solvePesuPvt(satpvt, currIdx, kernelPvt, 0);
        if(!kernelPvt.isValid) {
            currRes = 1e18; // 定位失败，视为极大残差
        } else {
            currRes = kernelPvt.rhor.squaredNorm();
        }
         
        // 更新最优核
        if (currRes < minResidual) {
            minResidual   = currRes;
            bestKernelIdx = currIdx;
            bestKernelPvt = kernelPvt;
        }

        /*
        std::cout << "Testing Kernel Indices: (";
        for (int idx : currIdx) {
            std::cout << idx << " ";
        }
        std::cout << ") Residual: " << currRes/(Sigma*Sigma) << " Threshold: " << thrKernel << std::endl;
        */
       
        if ((Mode == 2 || Mode == 4 ) && minResidual/(Sigma*Sigma) < thrKernel) {
            break; // 只要有一个满足即可
        }
        

    } while (std::next_permutation(combMask.begin(), combMask.end()));

    // 核集合残差超限，基础星座不可信，直接退出
    if (minResidual/(Sigma*Sigma) > thrKernel) {
        return;
    }

    // 打印全部的可信kenerl index
    /*
    std::cout << "Best Kernel Indices: ";
    for (int idx : bestKernelIdx) {
        std::cout << idx << " ";    
    }
    std::cout << " | Residual: " << minResidual/(Sigma*Sigma) << " | Threshold: " << thrKernel << std::endl; 
    */

    // ==========================================
    // Step 2：逐颗加入剩余卫星，验证自洽性
    // 自洽：保留；不自洽：标记为欺骗卫星
    // ==========================================
    std::vector<int> validIdx = bestKernelIdx;
    std::unordered_set<int> SpoofIdx;
    double resNew;

    if(Mode == 1 || Mode == 2){
        // 遍历所有不在核中的卫星索引
        for (int idx : goodIndices) {
            // 如果当前卫星已在可信集合中，跳过
            if (std::find(validIdx.begin(), validIdx.end(), idx) != validIdx.end()) {
                continue;
            }

            // 构造测试集合：可信集合 + 当前卫星
            std::vector<int> testIdx = validIdx;
            testIdx.push_back(idx);
            int satNumNew = testIdx.size();

            // 计算新组合残差
            //  = computeResidualSumSquared(solver, satpvt, testIdx);
            
            // 计算当前组合残差平方和
            solver.solvePesuPvt(satpvt, testIdx, kernelPvt, 0);
            if(!kernelPvt.isValid) {
                resNew = 1e18; // 定位失败，视为极大残差
            } else {
                resNew = kernelPvt.rhor.squaredNorm();
            }

            // 计算自适应门限：自由度 = 卫星数 - 4
            int dofNew     = satNumNew - 4;
            double thrNew  = gsl_cdf_chisq_Pinv(1.0 - Pfa, dofNew);

            if (resNew/(Sigma*Sigma) <= thrNew) {
                // 卫星正常，加入可信集合
                validIdx = testIdx;
            } else {
                // 卫星异常，标记为欺骗
                SpoofIdx.insert(satpvt.prn[idx]);
            }
        }
    }else{
        const double chi2_single = gsl_cdf_chisq_Pinv(1.0 - Pfa, 1);
        const double threshold   = (Sigma*Sigma) * chi2_single;

        // 遍历所有可用卫星，直接判断是否欺骗
        for (int Idx : goodIndices) {

            // 如果当前卫星已在可信集合中，跳过
            if (std::find(validIdx.begin(), validIdx.end(), Idx) != validIdx.end()) {
                continue;
            }

            // 构造测试集合：可信集合 + 当前卫星
            std::vector<int> testIdx = validIdx;
            testIdx.push_back(Idx);

            // 从卫星观测值中提取
            Eigen::Vector3d ps = satpvt.ps.row(Idx);
            double rho         = satpvt.rhos(Idx);
            double dt_s        = satpvt.dts(Idx);
            char sys           = satpvt.sys[Idx];

            // 用核PVT计算预测伪距
            Eigen::Vector3d r_es = ps - bestKernelPvt.pu;
            double pred_range    = r_es.norm();

            // 地球自转校正（和你的solvePesuPvt保持一致）
            double o_dot_e  = (sys == 'G') ? 7.2921151467e-5 : 7.292115e-5;
            double rot_corr = o_dot_e * (ps(0)*bestKernelPvt.pu(1) - ps(1)*bestKernelPvt.pu(0)) / consts::SPEED_OF_LIGHT;
            pred_range += rot_corr;

            // 加上钟差
            pred_range += consts::SPEED_OF_LIGHT * (bestKernelPvt.dtu - dt_s);

            // 单星残差
            double residual = pred_range - rho;
            double res2     = residual * residual;

            // ==========================================
            // 判决：残差平方 < 门限 → 真实卫星
            // ==========================================
            if (res2 < threshold) {
                // 卫星正常，加入可信集合
                validIdx = testIdx;
            } else {
                // 卫星异常，标记为欺骗
                SpoofIdx.insert(satpvt.prn[Idx]);
            }

            // std::cout << "Testing PRN: " << satpvt.prn[Idx] << " Residual^2: " << res2 << " Threshold: " << threshold << std::endl;
        }
    }

    // ======================
    // 3. 输出结果到 sat_results
    // ======================
    for (auto& sr : sat_results) {
        bool isSpoof = SpoofIdx.count(sr.prn) > 0;
        sr.is_spoof[9] = isSpoof;
        sr.score[9]    = isSpoof ? 1.0 : 0.0;
    }

}


// ============================================================================
// M11：多普勒域 RAIM 欺骗检测
// 参照 M9 伪距RAIM 完全同结构
// 依赖：多普勒PVT解算结果 dop_residual (Hz)
// 自由度：卫星数 - 7
// ============================================================================
void AntiSpoofing::detectMethod11(const UserPvtResult& userPvt, const int size, std::vector<SatSpoofResult>& sat_results)
{
    bool raimAlarm  = false;
    // 多普勒RAIM需要至少8颗卫星 + 解算有效
    if (size < 8 || !userPvt.dopValid) 
        return;
    
    double Pfa      = cfg_.anti_spoofing_.AntDRaim.Pfa;          // 虚警概率
    double Thr      = gsl_cdf_chisq_Pinv(1 - Pfa, size - 7);     // 自由度 = 卫星数 - 7
    double Sigma    = cfg_.anti_spoofing_.AntDRaim.Sigma;        // 多普勒观测噪声标准差（Hz）

    // 卡方检测统计量（与M9完全同格式）
    double test = userPvt.dop_residual.squaredNorm() / (Sigma * Sigma);

    std::cout << "RAIM Doppler Test: " << test << " Threshold: " << Thr << std::endl;

    raimAlarm   = test > Thr;

    // 保存结果到 M10 对应通道：is_spoof[10]、score[10]
    for (auto& sr : sat_results) {
        sr.is_spoof[10] = raimAlarm;
        sr.score[10]    = test ? test : 0.0;
    }
}


// ============================================================================
// M12：PVT 合理性检测
// 1. 高度过低报警
// 2. 速度过大报警
// 3. 钟漂过大报警
// 4. 相邻帧位置跳变报警
// ============================================================================
void AntiSpoofing::detectMethod12(const std::vector<PvtFrame>& pvt_history, std::vector<SatSpoofResult>& sat_results)
{
    bool alarm = false;

    // 阈值配置
    double MIN_HEIGHT        = cfg_.anti_spoofing_.AntPVT.MinHeight;             // 最小高度(m)
    double MAX_HEIGHT        = cfg_.anti_spoofing_.AntPVT.MaxHeight;             // 最大高度(m)
    double MAX_VELOCITY      = cfg_.anti_spoofing_.AntPVT.MaxVelocity;           // 最大速度(m/s)
    double MAX_CLOCK_DRIFT   = cfg_.anti_spoofing_.AntPVT.MaxClockDrift;         // 最大钟漂(s/s)
    double MAX_POS_JUMP      = cfg_.anti_spoofing_.AntPVT.MaxPositionChange;     // 最大位置跳变(m)

    if (pvt_history.empty()) return;

    // 当前帧 PVT
    const auto& curr     = pvt_history.back();
    const auto& curr_pvt = curr.pvt_full;

    if (!curr_pvt.isValid) return;

    // --------------------------
    // 1. 高度检测
    // --------------------------
    if (curr_pvt.pu_lla(2) < MIN_HEIGHT || curr_pvt.pu_lla(2) > MAX_HEIGHT) // 过高过低都报警
        alarm = true;

    // --------------------------
    // 2. 速度检测
    // --------------------------
    if (curr_pvt.vu.norm() > MAX_VELOCITY)
        alarm = true;

    // --------------------------
    // 3. 钟漂检测
    // --------------------------
    if (fabs(curr_pvt.ddtu) > MAX_CLOCK_DRIFT){
        alarm = true;
        // std::cout << "Clock Drift: " << curr_pvt.ddtu << " Threshold: " << MAX_CLOCK_DRIFT << std::endl;
    }
        

    // --------------------------
    // 4. 位置跳变检测（需要至少2帧）
    // --------------------------
    if (pvt_history.size() >= 2)
    {
        const auto& prev_pvt = pvt_history[pvt_history.size()-2].pvt_full;
        double posJump = (curr_pvt.pu - prev_pvt.pu).norm();
        if (posJump > MAX_POS_JUMP)
            alarm = true;
    }

    // --------------------------
    // 输出结果
    // --------------------------
    for (auto& sr : sat_results)
    {
        sr.is_spoof[11] = alarm;
        sr.score[11] = alarm ? 1.0 : 0.0;
    }
}

// ============================================================================
// M13：伪距速度 vs 多普勒速度 一致性检验（滑动窗卡方）
// 改进：只使用有效帧，无效帧跳过，动态凑数，动态自由度
// ============================================================================
void AntiSpoofing::detectMethod13(const std::vector<PvtFrame>& pvt_history, std::vector<SatSpoofResult>& sat_results)
{
    bool alarm      = false;
    double testStat = 0.0;

    // 从配置读取
    int    windowSize     = cfg_.anti_spoofing_.AntVelCon.WindowSize;               
    double Pfa            = cfg_.anti_spoofing_.AntVelCon.Pfa;
    double sigma          = cfg_.anti_spoofing_.AntVelCon.Sigma;

    // 至少需要1帧历史
    if (pvt_history.empty())
        return;

    // ==============================
    // 【鲁棒逻辑】从最新帧往前遍历，只收集有效帧，最多收集 windowSize 帧
    // ==============================
    std::vector<Eigen::Vector3d> valid_diffs;

    // 逆序遍历（最新 → 最旧）
    for (auto it = pvt_history.rbegin(); it != pvt_history.rend() && valid_diffs.size() < windowSize; ++it)
    {
        const auto& frame = *it;
        const auto& full_pvt = frame.pvt_full;

        // 只收集【同时有效】的伪距PVT + 多普勒PVT
        if (full_pvt.isValid && full_pvt.dopValid)
        {
            Eigen::Vector3d diff = full_pvt.vu - full_pvt.dop_vu;
            valid_diffs.push_back(diff);
        }
    }

    // 有效帧不足1帧 → 不检测
    int usedFrames = valid_diffs.size();
    if (usedFrames < 1)
        return;

    // ==============================
    // 构建残差向量 + 动态自由度
    // ==============================
    int dof = usedFrames * 3;
    Eigen::VectorXd residuals(dof);

    for (int i = 0; i < usedFrames; ++i)
    {
        residuals.segment<3>(i * 3) = valid_diffs[i];
    }

    // 卡方检验
    testStat = residuals.squaredNorm() / (sigma * sigma);
    double threshold = gsl_cdf_chisq_Pinv(1 - Pfa, dof);
    alarm = (testStat > threshold);

    // 输出结果
    for (auto& sr : sat_results)
    {
        sr.is_spoof[12] = alarm;
        sr.score[12] = testStat;
    }
}
