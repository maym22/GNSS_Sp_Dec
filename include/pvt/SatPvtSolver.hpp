#ifndef PVT_SOLVER_HPP
#define PVT_SOLVER_HPP

#include "config/ConfigReader.hpp"
#include "data_io/EphemerisReader.hpp"
#include "data_io/ObservationReader.hpp"
#include <vector>
#include <eigen3/Eigen/Dense>  // 推荐使用Eigen处理矩阵运算（需安装：sudo apt install libeigen3-dev）

// 对应Matlab pos_struct：PVT解算结果结构体
struct SatPVTStruct {
    // 卫星相关解算结果 
    Eigen::MatrixXd ps;       // 卫星位置 (m)，每行对应一颗卫星 [x,y,z]
    Eigen::MatrixXd vs;       // 卫星速度 (m/s)，每行对应一颗卫星 [vx,vy,vz]
    Eigen::VectorXd dts;      // 卫星钟差 (s)
    Eigen::VectorXd ddts;     // 卫星钟速 (s/s)
    //
    Eigen::VectorXd rhos;     // 伪距观测值 (m)
    Eigen::VectorXd drhos;    // 伪距率（多普勒）(m/s)
    std::vector<char>   sys;  // 系统
    std::vector<int>    prn;  // 卫星号

};

// PVT解算核心类（复用统一配置、星历/观测数据）
class SatPVTSolver {
private:
    const CfgReader& cfg_;               // 统一配置引用
    const std::vector<EphemerisData>& eph_raw_;  // 星历数据
    const double CLIGHT = 299792458.0;   // 光速 (m/s)

    // 私有辅助函数：单颗卫星轨道计算（对应ephPos核心逻辑）
    bool calcSatOrbit(const EphemerisData& eph, const ObsRawData& obs, 
                      Eigen::Vector3d& sat_pos, Eigen::Vector3d& sat_vel, 
                      double& dt, double& ddt);

public:
    // 构造函数：传入配置、星历数据
    SatPVTSolver(const CfgReader& cfg, const std::vector<EphemerisData>& eph_raw) 
        : cfg_(cfg), eph_raw_(eph_raw) {}

    // 核心接口：解算单个历元的PVT（对应Matlab ephPos）
    SatPVTStruct solveEpochPVT(const std::vector<ObsRawData>& obs_epoch);

    // 核心接口：批量解算所有历元的PVT（对应Matlab Step3循环）
    std::vector<SatPVTStruct> solveAllPVT(const std::vector<std::vector<ObsRawData>>& obs_raw);

    // 打印函数：验证PVT解算结果
    void printPVTResult(const std::vector<SatPVTStruct>& pvt_results) const;
};

#endif // SAT_PVT_SOLVER_HPP