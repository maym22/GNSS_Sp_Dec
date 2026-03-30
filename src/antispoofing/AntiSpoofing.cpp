#include "antispoofing/AntiSpoofing.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_cdf.h>
#include <chrono>

AntiSpoofing::AntiSpoofing(const CfgReader& cfg,
                           const std::vector<EphemerisData>& eph_raw,
                           const std::vector<std::vector<ObsRawData>>& obs_raw,
                           const std::vector<SatPVTStruct>& pvt_results)
    : cfg_(cfg), eph_raw_(eph_raw), obs_raw_(obs_raw), satpvt_results_(pvt_results)
{
    method_enabled_[0]  = cfg.anti_spoofing_.chnEphemNormativity.flag;
    method_enabled_[1]  = cfg.anti_spoofing_.chnSatVisibility.flag;
    method_enabled_[2]  = cfg.anti_spoofing_.chnCn0.flag;
    method_enabled_[3]  = cfg.anti_spoofing_.chnDoppCon.flag;
    method_enabled_[4]  = cfg.anti_spoofing_.chnClkDft.flag;
    method_enabled_[5]  = cfg.anti_spoofing_.chnRho.flag;
    method_enabled_[6]  = cfg.anti_spoofing_.mchnCn0Corr.flag;
    method_enabled_[7]  = cfg.anti_spoofing_.mchnDoppSD.flag;
    method_enabled_[8]  = cfg.anti_spoofing_.AntRaim.flag;
    method_enabled_[9]  = cfg.anti_spoofing_.AntGRaim.flag;
    method_enabled_[10] = cfg.anti_spoofing_.AntDRaim.flag;
    method_enabled_[11] = cfg.anti_spoofing_.AntPVT.flag;
    method_enabled_[12] = cfg.anti_spoofing_.AntVelCon.flag;
    method_enabled_[13] = cfg.anti_spoofing_.IMUVal.flag;
    method_enabled_[14] = cfg.anti_spoofing_.MantPVT.flag;
    method_enabled_[15] = cfg.anti_spoofing_.MantCPDD.flag;

}

// ===================== 更新所有卫星历史（优化版：有效计数 + 自动删除无效卫星） =====================
void AntiSpoofing::updateAllSatHistory(const std::vector<ObsRawData>& obs_epoch)
{
    // 1. 当前时刻卫星签到表
    std::unordered_map<int, const ObsRawData*> current_sat;
    for (const auto& obs : obs_epoch) {
        current_sat[obs.PRN] = &obs;
    }

    // 2. 更新所有已有卫星
    for (auto it = sat_history_.begin(); it != sat_history_.end();) {
        int prn = it->first;
        SatHistory& hist = it->second;

        bool exist = current_sat.count(prn);
        if (exist) {
            // 卫星存在 → 写入真实观测
            const auto& obs = *current_sat[prn];
            hist.push(true, obs.ObsTime, obs.CNR, obs.Rho, obs.Fd, obs.AcPh);
            current_sat.erase(prn);
            ++it;
        } else {
            // 卫星不存在 → 写入0
            hist.push(false, 0, 0, 0, 0, 0);

            // 如果完全无效 → 从map删除
            if (hist.isAllInvalid()) {
                it = sat_history_.erase(it);
            } else {
                ++it;
            }
        }
    }

    // 3. 添加新卫星
    for (const auto& pair : current_sat) {
        int prn = pair.first;
        const auto& obs = *pair.second;
        SatHistory hist;
        hist.push(true, obs.ObsTime, obs.CNR, obs.Rho, obs.Fd, obs.AcPh);
        sat_history_[prn] = hist;
    }
}

// ============================================================================
// 更新PVT历史缓存：固定最大帧数，自动丢弃最旧帧
// ============================================================================
void AntiSpoofing::updatePvtHistory(double timestamp, const UserPvtResult& pvt_full, const UserPvtResult& pvt_gps, const UserPvtResult& pvt_bds)
{
    PvtFrame frame;
    frame.timestamp = timestamp;
    frame.pvt_full  = pvt_full;
    frame.pvt_gps   = pvt_gps;
    frame.pvt_bds   = pvt_bds;

    // 添加新帧
    pvt_history_.push_back(frame);

    // 超过最大帧数 → 删除最旧的
    if (pvt_history_.size() > MAX_PVT_FRAMES)
    {
        pvt_history_.erase(pvt_history_.begin());
    }
}

// ========================== 工具函数 ==========================
const EphemerisData* AntiSpoofing::findEphemeris(int prn, char sys) const {
    for (const auto& e : eph_raw_)
        if (e.PRN == prn && e.sys == sys) return &e;
    return nullptr;
}

Eigen::Vector3d AntiSpoofing::lla2ecef(const Eigen::Vector3d& lla) const {
    const double a = 6378137.0, f = 1.0/298.257223563;
    double lat = lla(0)*M_PI/180, lon=lla(1)*M_PI/180, alt=lla(2);
    double N = a / sqrt(1 - f*(2-f)*sin(lat)*sin(lat));
    double x = (N+alt)*cos(lat)*cos(lon);
    double y = (N+alt)*cos(lat)*sin(lon);
    double z = (N*(1-f*f)+alt)*sin(lat);
    return {x,y,z};
}

std::pair<double, double> AntiSpoofing::cal_azel(const Eigen::Vector3d& rec_lla, const Eigen::Vector3d& sat_pos) const {

    Eigen::Vector3d rec_ecef = lla2ecef(rec_lla);
    double lat = rec_lla(0)*M_PI/180, lon=rec_lla(1)*M_PI/180;

    Eigen::Matrix3d R;
    R << -sin(lon),           cos(lon),          0,
         -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat),
          cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat);

    Eigen::Vector3d dr = sat_pos - rec_ecef;
    Eigen::Vector3d enu = R * dr;
    double e=enu(0), n=enu(1), u=enu(2);
    double az = atan2(e,n)*180/M_PI;
    if (az<0) az+=360;
    double el = atan2(u, hypot(e,n))*180/M_PI;

    return {az, el};
}

// ============================================================================
// 封装1：筛选【未被任何方法判骗】的卫星
// ============================================================================
void AntiSpoofing::selectGoodSatellites(
    const SatPVTStruct& satpvt, const std::vector<SatSpoofResult>& sat_results,
    std::vector<int>& NormalIdx, std::vector<int>& NormalIdx_GPS, std::vector<int>& NormalIdx_BDS)
{
    // 清空输出
    NormalIdx.clear();
    NormalIdx_GPS.clear();
    NormalIdx_BDS.clear();

    for (uint i = 0; i < satpvt.prn.size(); ++i) {
        int prn = satpvt.prn[i];
        bool spoofed = false;

        for (const auto& sr : sat_results) {
            if (sr.prn == prn) {
                spoofed = (sr.is_spoof[0] && method_enabled_[0]) || (sr.is_spoof[1] && method_enabled_[1]) ||
                          (sr.is_spoof[2] && method_enabled_[2]) || (sr.is_spoof[3] && method_enabled_[3]) ||
                          (sr.is_spoof[4] && method_enabled_[4]) || (sr.is_spoof[5] && method_enabled_[5]) ||
                          (sr.is_spoof[6] && method_enabled_[6]) || (sr.is_spoof[7] && method_enabled_[7]);
                break;
            }
        }

        if (!spoofed) {
            NormalIdx.push_back(i);

            // ========================
            // 【内置自动分离 GPS/BDS】
            // ========================
            char sys = satpvt.sys[i];
            if (sys == 'G')
            {
                NormalIdx_GPS.push_back(i);
            }
            else if (sys == 'C')
            {
                NormalIdx_BDS.push_back(i);
            }
        }
    }
    return;
}


// ========================== 主检测流程 ==========================
SpoofDetectOutput AntiSpoofing::runAllDetection() {
    SpoofDetectOutput output;

    // 按照观测量的每个历元进行检测obs_raw_.size()
    for (size_t i=0; i<5; i++) {

        const auto& obs_epoch = obs_raw_[i];
        const auto& satpvt    = satpvt_results_[i];

        // ===================== 每一历元：统一更新所有卫星历史 =====================
        updateAllSatHistory(obs_epoch);

        // 初始化单时刻欺骗检测结构体
        EpochSpoofResult spoof_result;
        spoof_result.epoch_idx = i+1;
        if (!obs_epoch.empty()) { 
            spoof_result.obs_time = obs_epoch[0].ObsTime; 
            spoof_result.time = obs_epoch[0].Time; 
        }
        for (int m=0; m<16; m++) {
            spoof_result.method_enable[m] = method_enabled_[m];
        }

        // ===================== 单通道算法：Single Channel =====================
        for (uint k = 0; k < obs_epoch.size(); ++k) {

            const auto& obs = obs_epoch[k];

            // 查找对应的星历数据
            const auto* eph = findEphemeris(obs.PRN, obs.sys);
            // 若无星历数据 → 跳过该卫星（无法进行后续检测）
            if (!eph) continue;
            // 查找对应的卫星位置
            Eigen::Vector3d sp = satpvt.ps.row(k); // 取第一颗卫星位置（实际应匹配对应卫星）

            SatSpoofResult sr;
            sr.prn = eph->PRN; 
            sr.sys = obs.sys;
            const SatHistory& hist = sat_history_[obs.PRN];  // 直接读历史

            // 单一时刻即可
            if (method_enabled_[0]) detectMethod01(*eph, sr);
            if (method_enabled_[1]) detectMethod02(sp  , sr);

            // 多时刻联合
            if (method_enabled_[2]) detectMethod03(hist, sr);
            if (method_enabled_[3]) detectMethod04(hist, sr);

            if (method_enabled_[4]) detectMethod05(hist, sr);
            if (method_enabled_[5]) detectMethod06(hist, sr);

            spoof_result.sat_results.push_back(sr);
        }

        // ===================== 多通道算法：Multi Channel =====================
        // 需要多时刻处理，只处理当前历元存在，且星历存在的卫星
        if (method_enabled_[6]) detectMethod07(spoof_result.sat_results);  // M7
        if (method_enabled_[7]) detectMethod08(spoof_result.sat_results);  // M8

        // ===================== PVT域算法：RAIM, PVT校验等 =====================
        // 根据前述欺骗检测结果，判断当前时刻中的可用卫星
        std::vector<int> NormalIdx, NormalIdx_gps, NormalIdx_bds; 
        selectGoodSatellites(satpvt, spoof_result.sat_results, NormalIdx, NormalIdx_gps, NormalIdx_bds);

        // 如果卫星总数目小于4颗，无法定位，因此不再进行检测
        if (NormalIdx.size() < 4) {
            std::cout  << " 总卫星数目：" << satpvt.prn.size() << " 有效卫星数目：" << NormalIdx.size() << std::endl;
            output.push_back(spoof_result);
            continue;
        }

        // 首先进行PVT结算：包括联合/单GPS/单BDS 的伪距定位和多普勒定位，并存入历史记录中
        UserPvtSolver solver;
        UserPvtResult pvt_full, pvt_gps, pvt_bds;
        // 进行分别解算
        int PesuPvtPrintFlag = 0;
        int DoppPvtPrintFlag = 0;
        solver.solvePesuPvt(satpvt, NormalIdx, pvt_full, PesuPvtPrintFlag);
        solver.solveDoppPVT(satpvt, NormalIdx, pvt_full, DoppPvtPrintFlag);
        //solver.solvePesuPvt(satpvt, NormalIdx_gps, pvt_gps);
        //solver.solveDoppPVT(satpvt, NormalIdx_gps, pvt_gps);
        //solver.solvePesuPvt(satpvt, NormalIdx_bds, pvt_bds);
        //solver.solveDoppPVT(satpvt, NormalIdx_bds, pvt_bds);
        // 存入历史
        updatePvtHistory(obs_epoch[0].ObsTime, pvt_full, pvt_gps, pvt_bds);
        
        // ==========================
        // 先执行一次所有轻量级检测，标记卫星结果（M1-M9, M11-M13）
        // ==========================
        // RAIM 检测
        if (method_enabled_[8]) detectMethod09(pvt_full, NormalIdx.size(), spoof_result.sat_results);
        // Doppler RAIM检测
        if (method_enabled_[10]) detectMethod11(pvt_full, NormalIdx.size(), spoof_result.sat_results);
        // M12 PVT结果合理性检测（待实现）
        if (method_enabled_[11]) detectMethod12(getPvtHistory(), spoof_result.sat_results);
        // M13 DVD与IVD测速（待实现）
        if (method_enabled_[12]) detectMethod13(getPvtHistory(), spoof_result.sat_results);

        // ==========================
        // 【重量级 M10】只有前面任意一个报警时才执行
        // ==========================
        if (method_enabled_[9]){

            bool need_heavy_check = false;

            // 遍历卫星结果：只要有一颗卫星在 M9/M11/M12/M13 任意一个报警，就开启 M10
            const auto& sr = spoof_result.sat_results[0];
            if (sr.is_spoof[8] || sr.is_spoof[10] || sr.is_spoof[11] || sr.is_spoof[12]) {
                need_heavy_check = true;
            }

            // 只有需要时才执行 生成式RAIM检测
            auto start = std::chrono::high_resolution_clock::now();

            if (need_heavy_check) detectMethod10(satpvt, NormalIdx, spoof_result.sat_results);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> duration = end - start; // 毫秒

            // 打印耗时
            std::cout << " M10 Generative RAIM 已执行 | 耗时: " 
                    << duration.count() << " ms" << std::endl;
        }

        // M14 GNSS与IMU进行比较（待实现）

        // M15 双天线PVT对比

        // M16 双天线载波相位双差

        output.push_back(spoof_result);
    }
    return output;
}

// ========================== 打印 ==========================
void AntiSpoofing::printSpoofResult(const SpoofDetectOutput& out) const {
    std::cout << "\n========== Anti-Spoofing Result ==========\n";
    for (const auto& e : out) {
        std::cout << "\nEpoch " << e.epoch_idx << " | TOW=" << e.obs_time << "\n";
        for (const auto& s : e.sat_results) {
            std::cout << "  PRN" << std::setw(3) << s.prn << "(" << s.sys << "): ";
            if (e.method_enable[0])  std::cout << "M1:"  << (s.is_spoof[0]?"Sp":"Ok") << ", ";
            if (e.method_enable[1])  std::cout << "M2:"  << (s.is_spoof[1]?"Sp":"Ok") << ", ";
            if (e.method_enable[2])  std::cout << "M3:"  << (s.is_spoof[2]?"Sp":"Ok") << ", ";
            if (e.method_enable[3])  std::cout << "M4:"  << (s.is_spoof[3]?"Sp":"Ok") << ", ";
            if (e.method_enable[4])  std::cout << "M5:"  << (s.is_spoof[4]?"Sp":"Ok") << ", ";
            if (e.method_enable[5])  std::cout << "M6:"  << (s.is_spoof[5]?"Sp":"Ok") << ", ";
            if (e.method_enable[6])  std::cout << "M7:"  << (s.is_spoof[6]?"Sp":"Ok") << ", ";
            if (e.method_enable[7])  std::cout << "M8:"  << (s.is_spoof[7]?"Sp":"Ok") << ", ";
            if (e.method_enable[8])  std::cout << "M9:"  << (s.is_spoof[8]?"Sp":"Ok") << ", ";
            if (e.method_enable[9])  std::cout << "M10:" << (s.is_spoof[9]? "Sp":"Ok") << ", ";
            if (e.method_enable[10]) std::cout << "M11:" << (s.is_spoof[10]?"Sp":"Ok") << ", ";
            if (e.method_enable[11]) std::cout << "M12:" << (s.is_spoof[11]?"Sp":"Ok") << ", ";
            if (e.method_enable[12]) std::cout << "M13:" << (s.is_spoof[12]?"Sp":"Ok") << ", ";
            if (e.method_enable[13]) std::cout << "M14:" << (s.is_spoof[13]?"Sp":"Ok") << ", ";
            if (e.method_enable[14]) std::cout << "M15:" << (s.is_spoof[14]?"Sp":"Ok") << ", ";
            if (e.method_enable[15]) std::cout << "M16:" << (s.is_spoof[15]?"Sp":"Ok") << ", ";

            std::cout << "\n";
        }
    }
}