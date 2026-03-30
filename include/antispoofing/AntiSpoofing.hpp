#ifndef ANTI_SPOOFING_HPP
#define ANTI_SPOOFING_HPP

#include "config/ConfigReader.hpp"
#include "data_io/EphemerisReader.hpp"
#include "data_io/ObservationReader.hpp"
#include "pvt/SatPvtSolver.hpp"
#include "pvt/UserPvtSolver.hpp"
#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>
#include <cmath>

constexpr int DETECT_METHOD_NUM = 16;


// ===================== 卫星历史（优化版：有效计数 + 自动清理） =====================
struct SatHistory {
    static constexpr int MAX_SIZE = 50;  // 最大队列长度

    std::vector<double> time;
    std::vector<double> cn0;
    std::vector<double> rho;
    std::vector<double> dopp;
    std::vector<double> acph;
    std::vector<int> effective;

    // 入队（自动裁剪）
    void push(bool valid, double t, double c, double r, double d, double a) {
        time.push_back(valid ? t : 0.0);
        cn0.push_back(valid ? c : 0.0);
        rho.push_back(valid ? r : 0.0);
        dopp.push_back(valid ? d : 0.0);
        acph.push_back(valid ? a : 0.0);
        effective.push_back(valid ? 1 : 0);

        // 超长度 → 弹出最旧
        if ((int)time.size() > MAX_SIZE) {
            time.erase(time.begin());
            cn0.erase(cn0.begin());
            rho.erase(rho.begin());
            dopp.erase(dopp.begin());
            acph.erase(acph.begin());
            effective.erase(effective.begin());
        }
    }

    // 队列总长度（不变）
    int size() const {
        return (int)time.size();
    }

    // 真正有效数据数量（只统计effective=1）
    int validCount() const {
        int cnt = 0;
        for (int e : effective) {
            if (e == 1) cnt++;
        }
        return cnt;
    }

    // 判断是否完全无效（可删除）
    bool isAllInvalid() const {
        return validCount() == 0;
    }
};


// 单卫星检测结果
struct SatSpoofResult {
    int   prn;
    char  sys;
    bool  is_spoof[DETECT_METHOD_NUM] = {false};  // 每个方法的判决结果
    double score[DETECT_METHOD_NUM] = {0.0};      // 每个方法的得分（可用于线性判决）
};

// 单历元检测结果
struct EpochSpoofResult {
    // Time Information
    int epoch_idx;
    double obs_time;
    std::vector<int> time;

    // Method Enable Flags
    bool method_enable[DETECT_METHOD_NUM];
    std::vector<SatSpoofResult> sat_results;
};

using SpoofDetectOutput = std::vector<EpochSpoofResult>;

// 历史PVT帧：包含全卫星/GPS/BDS 三组解算结果
struct PvtFrame
{
    double          timestamp;      // 保留时间戳（不参与判断，但方便调试）
    UserPvtResult   pvt_full;       // 全部卫星解算 PVT
    UserPvtResult   pvt_gps;        // 仅GPS卫星解算 PVT
    UserPvtResult   pvt_bds;        // 仅BDS卫星解算 PVT
};


class AntiSpoofing {
private:
    const CfgReader& cfg_;
    const std::vector<EphemerisData>& eph_raw_;
    const std::vector<std::vector<ObsRawData>>& obs_raw_;
    const std::vector<SatPVTStruct>& satpvt_results_;

    bool method_enabled_[DETECT_METHOD_NUM];
    const double CLIGHT = 299792458.0;

    // ===================== 全局卫星历史 =====================
    std::unordered_map<int, SatHistory> sat_history_;

    // ===================== 全局PVT历史 =====================
    std::vector<PvtFrame> pvt_history_;
    // 最大缓存帧数（可改：10/20/30）
    const int MAX_PVT_FRAMES = 20;   
    // 新增：更新PVT历史（自动维护长度）
    void updatePvtHistory(double timestamp, const UserPvtResult& pvt_full, const UserPvtResult& pvt_gps,const UserPvtResult& pvt_bds);

    // -------------------------- 工具函数 --------------------------
    const EphemerisData* findEphemeris(int prn, char sys) const;
    Eigen::Vector3d lla2ecef(const Eigen::Vector3d& lla) const;
    std::pair<double, double> cal_azel(const Eigen::Vector3d& rec_lla, const Eigen::Vector3d& sat_pos) const;

    // -------------------------- 16 个检测算法 --------------------------
    // M1 星历规范性
    void detectMethod01(const EphemerisData& eph, SatSpoofResult& res);
    // M2 卫星可见性（仰角）
    void detectMethod02(const Eigen::Vector3d& ps, SatSpoofResult& res);
    // M3 载噪比 CNR
    void detectMethod03(const SatHistory& hist, SatSpoofResult& res);
    // M4 多普勒一致性
    void detectMethod04(const SatHistory& hist, SatSpoofResult& res);
    // M5 单星钟漂变化
    void detectMethod05(const SatHistory& hist, SatSpoofResult& res);
    // M6 发射时间合理性（伪距范围）
    void detectMethod06(const SatHistory& hist, SatSpoofResult& res);

    // M7~M8 多通道联合检测
    void detectMethod07(std::vector<SatSpoofResult>& sat_results);
    void detectMethod08(std::vector<SatSpoofResult>& sat_results);

    // PVT域检测：RAIM、PVT校验等
    void detectMethod09(const UserPvtResult& userPvt, const int size, std::vector<SatSpoofResult>& sat_results);
    void detectMethod10(const SatPVTStruct& satpvt, const std::vector<int>& goodIndices, std::vector<SatSpoofResult>& sat_results);
    void detectMethod11(const UserPvtResult& userPvt, const int size, std::vector<SatSpoofResult>& sat_results); 

    void detectMethod12(const std::vector<PvtFrame>& pvt_history, std::vector<SatSpoofResult>& sat_results);
    void detectMethod13(const std::vector<PvtFrame>& pvt_history, std::vector<SatSpoofResult>& sat_results);

    void detectMethod14(...) {}
    void detectMethod15(...) {}
    void detectMethod16(...) {}

public:
    AntiSpoofing(const CfgReader& cfg,
                 const std::vector<EphemerisData>& eph_raw,
                 const std::vector<std::vector<ObsRawData>>& obs_raw,
                 const std::vector<SatPVTStruct>& pvt_results);
    void selectGoodSatellites(const SatPVTStruct& satpvt, const std::vector<SatSpoofResult>& sat_results,
                              std::vector<int>& NormalIdx, std::vector<int>& NormalIdx_gps, std::vector<int>& NormalIdx_bds);
    void extractNormalObservations(const SatPVTStruct& satpvt, const std::vector<int>& goodIndices,
                                    Eigen::MatrixXd& goodPos, Eigen::MatrixXd& goodVel, Eigen::VectorXd& goodClk, 
                                    Eigen::VectorXd& goodRho, Eigen::VectorXd& goodDrho, 
                                    std::vector<int>& goodPrns, std::vector<char>& goodSys);

    // 【关键】每一历元统一更新所有卫星历史
    void updateAllSatHistory(const std::vector<ObsRawData>& obs_epoch);
    // 获取PVT历史
    const std::vector<PvtFrame>& getPvtHistory() const { return pvt_history_; }

    SpoofDetectOutput runAllDetection();
    void printSpoofResult(const SpoofDetectOutput& output) const;
};

#endif