# pragma once
#ifndef USER_PVT_SOLVER_HPP
#define USER_PVT_SOLVER_HPP

#include <eigen3/Eigen/Dense>
#include <vector>
#include "pvt/Consts.hpp"   // 你工程里的常量定义：c、地球自转、椭球参数等
#include "pvt/SatPvtSolver.hpp" // 卫星PVT解算结果结构体 

struct UserPvtResult
{
    Eigen::Vector3d pu;             // ECEF 位置
    Eigen::Vector3d vu;             // ECEF 速度
    Eigen::Vector3d pu_lla;         // 纬度/经度/高度
    Eigen::Vector3d vu_enu;         // ENU 速度
    double dtu = 0.0;               // 钟差 (s)
    double ddtu = 0.0;              // 钟漂 (s/s)
    Eigen::VectorXd rhor;           // 伪距残差
    Eigen::VectorXd drhor;          // 伪距率残差
    bool isValid = false;
    std::vector<int>    prns;       // 解算卫星信息      
    std::vector<char>   sys;

    // ======================
    // 【新增：多普勒PVT（完整位置+速度+钟漂）】
    // ======================
    Eigen::Vector3d dop_pu;         // 多普勒解算 ECEF 位置
    Eigen::Vector3d dop_vu;         // 多普勒解算 ECEF 速度
    Eigen::Vector3d dop_pu_lla;     // 多普勒 LLA
    Eigen::Vector3d dop_vu_enu;     // 多普勒 ENU 速度
    double          dop_ddtu;       // 多普勒钟漂 (s/s)
    Eigen::VectorXd dop_residual;   // 多普勒残差（用于Doppler RAIM）
    Eigen::MatrixXd dop_H;          // 多普勒观测矩阵 H
    bool dopValid = false;          // 多普勒解算有效标志

    std::vector<int>    dop_prns;       // 解算卫星信息      
    std::vector<char>   dop_sys;
    
};

class UserPvtSolver
{
public:
    void solvePesuPvt(const SatPVTStruct& satpvt, const std::vector<int>& goodIndices, UserPvtResult& result, int printflag);
    void solveDoppPVT(const SatPVTStruct& satpvt, const std::vector<int>& goodIndices, UserPvtResult& result, int printflag);

    void printPesuPvtDebug(const UserPvtResult& res);
    void printDoppPvtDebug(const UserPvtResult& res);
private:
    Eigen::MatrixXd calGeoMat(const Eigen::MatrixXd& ps, const Eigen::Vector3d& pu);
    Eigen::VectorXd rotCorrection(const Eigen::MatrixXd& ps, const Eigen::Vector3d& pu);
    Eigen::MatrixXd calLoS(const Eigen::MatrixXd& ps, const Eigen::Vector3d& pu);
    Eigen::MatrixXd vec3d_crossProduct(const Eigen::MatrixXd& v1, const Eigen::MatrixXd& v2);
    Eigen::Vector3d ecef2lla(const Eigen::Vector3d& ecef);
    Eigen::Vector3d ecefVel2EnuVel(const Eigen::Vector3d& vu_ecef, const Eigen::Vector3d& lla); 
    
};

#endif