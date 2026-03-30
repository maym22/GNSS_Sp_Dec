#include "pvt/SatPvtSolver.hpp"
#include "pvt/Consts.hpp"   // 你工程里的常量定义：c、地球自转、椭球参数等
#include <stdexcept>
#include <cmath>

// 单颗卫星轨道计算（星历→位置/速度/钟差，对应ephPos核心逻辑）
bool SatPVTSolver::calcSatOrbit(const EphemerisData& eph, const ObsRawData& obs,
                             Eigen::Vector3d& sat_pos, Eigen::Vector3d& sat_vel,
                             double& dt, double& ddt)
{
    // ===================== 常量定义（与 MATLAB consts 完全一致）=====================
    const double mu_gps     = 3.9860050e14;
    const double mu_bds      = 3.986004418e14;
    const double o_dot_e_gps = 7.2921151467e-5;
    const double o_dot_e_bds = 7.292115e-5;
    const double F          = -4.442807633e-10; // 相对论修正系数
    const double CLIGHT     = 299792458.0;

    // ===================== 钟差 eph2clk（完全 MATLAB 迭代）=====================
    auto eph2clk = [&](double tsv_clk, const EphemerisData& eph_clk) -> std::pair<double, double> {
        double ts = tsv_clk - eph_clk.Toc_TOW;
        double t_ = ts;
        double t = ts - (eph_clk.af0 + eph_clk.af1 * t_ + eph_clk.af2 * t_ * t_);

        while (fabs(t - t_) > 1e-15) {
            t_ = t;
            t = ts - (eph_clk.af0 + eph_clk.af1 * t_ + eph_clk.af2 * t_ * t_);
        }

        double dt_out, ddt_out;
        if (eph_clk.sys == 'G') {
            dt_out = eph_clk.af0 + eph_clk.af1 * t + eph_clk.af2 * t*t - eph_clk.TGD;
        } else if (eph_clk.sys == 'C') {
            dt_out = eph_clk.af0 + eph_clk.af1 * t + eph_clk.af2 * t*t - 2 * eph_clk.TGD;
        } else {
            dt_out = 0.0;
        }
        ddt_out = eph_clk.af1 + 2 * eph_clk.af2 * t;
        return {dt_out, ddt_out};
    };

    // 发射时间 tsv（与 MATLAB 完全一致：观测时间 - 伪距/c）
    double tsv = obs.ObsTime - obs.Rho / CLIGHT - 14*(obs.sys == 'C');
    auto [tsv_1, tsv_2] = eph2clk(tsv, eph);
    tsv  = tsv - tsv_1;

    // ===================== 1. 调用 eph2pt 计算卫星位置 + 钟差 =====================
    auto eph2pt = [&](double tsv_in) -> Eigen::Vector4d {
        double mu, o_dot_e;
        if (eph.sys == 'G') {
            mu      = mu_gps;
            o_dot_e = o_dot_e_gps;
        } else if (eph.sys == 'C') {
            mu      = mu_bds;
            o_dot_e = o_dot_e_bds;
        } else {
            mu      = mu_gps;
            o_dot_e = o_dot_e_gps;
        }

        // 星历参数（完全与 MATLAB 对应）
        double Crs     = eph.Crs;
        double Delta_n = eph.Delta_n;
        double M_0     = eph.M_0;
        double Cuc     = eph.Cuc;
        double e       = eph.e;
        double Cus     = eph.Cus;
        double sqrt_a  = eph.sqrt_a;
        double Toe     = eph.Toe;
        double Cic     = eph.Cic;
        double Omega_0 = eph.Omega_0;
        double Cis     = eph.Cis;
        double i_0     = eph.i_0;
        double Crc     = eph.Crc;
        double omega   = eph.omega;
        double Omega_dot = eph.Omega_dot;
        double I_dot   = eph.I_dot;

        // 时间 tk（完全 MATLAB 逻辑）
        double tk = tsv_in - Toe;
        if (tk > 302400)  tk -= 604800;
        if (tk < -302400) tk += 604800;

        // 平均角速度 n
        double n = sqrt(mu / pow(sqrt_a, 6)) + Delta_n;

        // 平近点角 Mk
        double Mk = M_0 + n * tk;

        // 偏近点角迭代（完全 MATLAB 格式）
        double Ek_ = Mk;
        double Ek = Ek_ - (Ek_ - e * sin(Ek_) - Mk) / (1 - e * cos(Ek_));
        while (fabs(Ek - Ek_) > 1e-10) {
            Ek_ = Ek;
            Ek = Ek_ - (Ek_ - e * sin(Ek_) - Mk) / (1 - e * cos(Ek_));
        }

        // 真近点角 vk
        double vk = atan2(sqrt(1 - e*e) * sin(Ek), cos(Ek) - e);

        // 升交角距 phik
        double phik = vk + omega;

        // 摄动修正
        double delta_uk = Cus * sin(2 * phik) + Cuc * cos(2 * phik);
        double delta_rk = Crs * sin(2 * phik) + Crc * cos(2 * phik);
        double delta_ik = Cis * sin(2 * phik) + Cic * cos(2 * phik);

        // 修正后轨道参数
        double uk  = phik + delta_uk;
        double rk  = sqrt_a*sqrt_a * (1 - e*cos(Ek)) + delta_rk;
        double ik  = i_0 + I_dot * tk + delta_ik;

        // 轨道面坐标
        double xk = rk * cos(uk);
        double yk = rk * sin(uk);

        // 升交点赤经 Omegak（MATLAB 原版）
        double Omegak = Omega_0 + (Omega_dot - o_dot_e) * tk
                            - o_dot_e * (Toe - 0 * (eph.sys == 'C'));

        // ECEF 位置，未修正地球自转
        double px = xk * cos(Omegak) - yk * cos(ik) * sin(Omegak);
        double py = xk * sin(Omegak) + yk * cos(ik) * cos(Omegak);
        double pz = yk * sin(ik);

        // 卫星钟差
        auto [dtsv, ddtsv] = eph2clk(tsv_in, eph);
        double dts = dtsv + F * e * sqrt_a * sin(Ek); // 相对论修正

        return Eigen::Vector4d(px, py, pz, dts);
    };

    // ===================== 2. 调用 eph2pvt 速度：差分 1ms =====================
    // 其实也可以直接计算速度，但差分更简答
    double dt_vel = 1e-3;
    Eigen::Vector4d pt0 = eph2pt(tsv);
    Eigen::Vector4d pt1 = eph2pt(tsv + dt_vel);

    sat_pos[0] = pt0(0);
    sat_pos[1] = pt0(1);
    sat_pos[2] = pt0(2);

    sat_vel[0] = (pt1(0) - pt0(0)) / dt_vel;
    sat_vel[1] = (pt1(1) - pt0(1)) / dt_vel;
    sat_vel[2] = (pt1(2) - pt0(2)) / dt_vel;

    // 最终钟差 & 钟漂
    dt  = pt0(3);
    ddt = (eph2pt(tsv + dt_vel)(3) - pt0(3)) / dt_vel;

    return true;
}

// 解算单个历元的卫星PVT（对应Matlab ephPos）
SatPVTStruct SatPVTSolver::solveEpochPVT(const std::vector<ObsRawData>& obs_epoch) {
    SatPVTStruct pvt;
    int valid_sat = 0;

    // 1. 初始化矩阵/向量（预留空间）
    pvt.sys.reserve(obs_epoch.size());
    pvt.prn.reserve(obs_epoch.size());
    pvt.rhos.resize(obs_epoch.size());
    pvt.drhos.resize(obs_epoch.size());
    // 卫星相关
    pvt.ps.resize(obs_epoch.size(), 3);
    pvt.vs.resize(obs_epoch.size(), 3);
    pvt.dts.resize(obs_epoch.size());
    pvt.ddts.resize(obs_epoch.size());


    // 2. 遍历当前历元的所有观测卫星
    for (size_t i = 0; i < obs_epoch.size(); ++i) {
        const ObsRawData& obs = obs_epoch[i];
        // 匹配对应卫星的星历
        auto it = std::find_if(eph_raw_.begin(), eph_raw_.end(), 
            [&](const EphemerisData& eph) {
                return eph.PRN == obs.PRN && eph.sys == obs.sys;
            });
        
        if (it == eph_raw_.end()) continue; // 无星历数据跳过

        // 计算卫星位置/速度/钟差
        Eigen::Vector3d sat_pos, sat_vel;
        double dt, ddt;
        if (!calcSatOrbit(*it, obs, sat_pos, sat_vel, dt, ddt)) continue;

        // 存储结果
        pvt.ps.row(valid_sat) = sat_pos;
        pvt.vs.row(valid_sat) = sat_vel;
        pvt.dts(valid_sat)    = dt;
        pvt.ddts(valid_sat)   = ddt;
        pvt.rhos(valid_sat)   = obs.Rho; // 伪距修正钟差  - CLIGHT * dt
        pvt.drhos(valid_sat)  = - obs.Fd * CLIGHT / obs.Fc; // 多普勒转伪距率

        pvt.prn.push_back(obs.PRN);
        pvt.sys.push_back(obs.sys);

        valid_sat++;
    }

    // 3. 截断到有效卫星数量
    pvt.ps.conservativeResize(valid_sat, 3);
    pvt.vs.conservativeResize(valid_sat, 3);
    pvt.dts.conservativeResize(valid_sat);
    pvt.ddts.conservativeResize(valid_sat);
    pvt.rhos.conservativeResize(valid_sat);
    pvt.drhos.conservativeResize(valid_sat);
    pvt.prn.shrink_to_fit();
    pvt.sys.shrink_to_fit();    

    return pvt;
}

// 批量解算所有历元的PVT（对应Matlab Step3循环）
std::vector<SatPVTStruct> SatPVTSolver::solveAllPVT(const std::vector<std::vector<ObsRawData>>& obs_data) {
    std::vector<SatPVTStruct> pvt_results;
    pvt_results.reserve(obs_data.size());

    // 遍历每个历元解算PVT
    for (size_t o = 0; o < obs_data.size(); ++o) {
        const auto& obs_epoch = obs_data[o];
        if (obs_epoch.empty()) {
            pvt_results.emplace_back();
            continue;
        }
        pvt_results.push_back(solveEpochPVT(obs_epoch));
        
    }

    // std::cout << "[PVT解算] 完成 " << pvt_results.size() << " 个历元的SatPVT解算" << std::endl;
    return pvt_results;
}

// 打印PVT解算结果（验证正确性）
void SatPVTSolver::printPVTResult(const std::vector<SatPVTStruct>& pvt_results) const {
    if (pvt_results.empty()) {
        std::cout << "\n[PVT打印] 无有效PVT解算结果！" << std::endl;
        return;
    }

    std::cout << "\n======================================= SatPVT解算结果 =======================================" << std::endl;
    for (size_t epoch = 0; epoch < 1; ++epoch) { // pvt_results.size()
        const auto& pvt = pvt_results[epoch];

        // 打印卫星解算结果
        std::cout << "卫星数量：" << pvt.ps.rows() << std::endl;
        std::cout << std::setw(6) << "PRN" 
                  << std::setw(15) << "X (m)" 
                  << std::setw(15) << "Y (m)" 
                  << std::setw(15) << "Z (m)" 
                  << std::setw(15) << "vX (m/s)" 
                  << std::setw(15) << "vY (m/s)" 
                  << std::setw(15) << "vZ (m/s)" 
                  << std::setw(12) << "钟差 (ns)" << std::endl;
        std::cout << "----------------------------------------------------------------------------------------" << std::endl;

        for (int i = 0; i < pvt.ps.rows(); ++i) {
            std::cout << std::setw(6)  << pvt.prn[i]
                      << std::setw(15) << std::fixed << std::setprecision(3) << pvt.ps(i,0)
                      << std::setw(15) << std::fixed << std::setprecision(3) << pvt.ps(i,1)
                      << std::setw(15) << std::fixed << std::setprecision(3) << pvt.ps(i,2)
                      << std::setw(15) << std::fixed << std::setprecision(3) << pvt.vs(i,0)
                      << std::setw(15) << std::fixed << std::setprecision(3) << pvt.vs(i,1)
                      << std::setw(15) << std::fixed << std::setprecision(3) << pvt.vs(i,2)
                      << std::setw(15) << std::fixed << std::setprecision(3) << pvt.rhos(i) << std::endl;
        }
    }
    std::cout << "\n=========================================================================================" << std::endl;
}