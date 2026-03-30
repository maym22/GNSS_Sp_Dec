#include "pvt/UserPvtSolver.hpp"
#include <iostream>
#include <iomanip>

void UserPvtSolver::solvePesuPvt(const SatPVTStruct& satpvt, const std::vector<int>& goodIndices, UserPvtResult& result, int printflag)
{
    int M = goodIndices.size();

    if (M < 4) {
        result.isValid = false;
        if (printflag == 1) {
            printPesuPvtDebug(result);
        }
        return;
    }
    // ==========================
    // 内部自动提取优质观测（完全封装）
    // ==========================
    Eigen::MatrixXd ps(M,3), vs(M,3);
    Eigen::VectorXd dts(M), rhos(M), drhos(M);
    std::vector<int> Prns;
    std::vector<char> Sys;

    for (int k=0; k<M; ++k) {
        int i     = goodIndices[k];
        ps.row(k) = satpvt.ps.row(i);
        vs.row(k) = satpvt.vs.row(i);
        dts(k)    = satpvt.dts(i);
        rhos(k)   = satpvt.rhos(i);
        drhos(k)  = satpvt.drhos(i);
        Prns.push_back(satpvt.prn[i]);
        Sys.push_back(satpvt.sys[i]);
    }

    // 初始化
    Eigen::Vector3d pu = Eigen::Vector3d::Zero();
    double dtu         = 0.0;
    double error       = 1e5;
    const int iter_max = 100;
    int iter_cnt       = 0;

    // 迭代最小二乘解位置
    while (error > 1e-5 && iter_cnt < iter_max) {
        Eigen::MatrixXd H = calGeoMat(ps, pu);
        Eigen::VectorXd r = rotCorrection(ps, pu);
        Eigen::VectorXd b_pos = (r.array() + consts::SPEED_OF_LIGHT * (dtu - dts.array())) - rhos.array();

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_H(H);
        Eigen::VectorXd dx_pos = -qr_H.solve(b_pos);

        pu  += dx_pos.head(3);
        dtu += dx_pos(3) / consts::SPEED_OF_LIGHT;

        error = dx_pos.cwiseAbs().maxCoeff();
        iter_cnt++;
    }

    // 最终残差
    Eigen::MatrixXd H_final = calGeoMat(ps, pu);
    Eigen::VectorXd r_final = rotCorrection(ps, pu);
    Eigen::VectorXd rhor = (r_final.array() + consts::SPEED_OF_LIGHT * (dtu - dts.array())) - rhos.array();

    // 解速度
    Eigen::MatrixXd LoS = calLoS(ps, pu);
    Eigen::MatrixXd G(LoS.rows(), 4);
    G << -LoS, Eigen::VectorXd::Ones(LoS.rows());

    Eigen::VectorXd b_vel = (vs.array() * (-LoS).array()).rowwise().sum() + drhos.array();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_G(G);
    Eigen::VectorXd dx_vel = qr_G.solve(b_vel);

    Eigen::Vector3d vu = dx_vel.head(3);
    double ddtu = dx_vel(3);
    Eigen::VectorXd drhor = b_vel;

    // 输出结果
    result.pu      = pu;
    result.vu      = vu;
    result.pu_lla  = ecef2lla(pu);
    result.vu_enu  = ecefVel2EnuVel(vu, result.pu_lla);
    result.dtu     = dtu;
    result.ddtu    = ddtu/consts::SPEED_OF_LIGHT; // 转换回秒/秒
    result.rhor    = rhor;
    result.drhor   = drhor;
    result.isValid = true;
    result.prns    = Prns;      // 自动保存
    result.sys     = Sys;       // 自动保存

    if (printflag == 1) {
        printPesuPvtDebug(result);
    }
    return;
}


// ============================================================================
// 【统一架构】多普勒PVT解算（位置+速度+钟漂）
// 结构与 solvePvtIterative 完全一致
// 输入：drho (m/s)
// 输出：dop_residual (Hz)
// 状态量：7维 = [pos3, vel3, ddtu1]
// 最小卫星数：7
// ============================================================================
void UserPvtSolver::solveDoppPVT(const SatPVTStruct& satpvt, const std::vector<int>& goodIndices, UserPvtResult& result, int printflag)
{
    int M = goodIndices.size();
    const int N_STATE = 7;

    // 至少7颗卫星
    if (M < N_STATE) {
        result.dopValid = false;
        if (printflag == 1) {
            printDoppPvtDebug(result);
        }
        return;
    }

    // ===========================
    // 1. 提取观测（统一格式）
    // ===========================
    Eigen::MatrixXd goodPos(M, 3);
    Eigen::MatrixXd goodVel(M, 3);
    Eigen::VectorXd goodDrho(M);   // m/s
    Eigen::VectorXd goodDdts(M);
    std::vector<int> goodPrns;
    std::vector<char> goodSys;

    for (int k = 0; k < M; ++k) {
        int i = goodIndices[k];
        goodPos.row(k) = satpvt.ps.row(i);
        goodVel.row(k) = satpvt.vs.row(i);
        goodDrho(k)    = satpvt.drhos(i);
        goodDdts(k)    = satpvt.ddts(i);
        goodPrns.push_back(satpvt.prn[i]);
        goodSys.push_back(satpvt.sys[i]);
    }

    // ===========================
    // 单位转换：m/s → Hz
    // dopp = - drho * f0 / c
    // ===========================
    double c = consts::SPEED_OF_LIGHT;
    Eigen::VectorXd f0(M);
    for (int i = 0; i < M; ++i) {
        f0(i) = (goodSys[i] == 'G') ? consts::GPS_FREQ : consts::BDS_FREQ;
    }
    Eigen::VectorXd dopp = -goodDrho.array() * f0.array() / c;

    // ===========================
    // 2. 迭代初始化（与伪距PVT完全相同）
    // ===========================
    Eigen::Vector3d pu_k = Eigen::Vector3d::Zero();
    Eigen::Vector3d vu_k = Eigen::Vector3d::Zero();
    double ddtu_k = 0.0;

    double error_k = 1e5;
    const int iter_max = 100;
    int iter_cnt = 0;

    Eigen::MatrixXd H(M, N_STATE);
    Eigen::VectorXd b(M);

    // ===========================
    // 3. 迭代最小二乘（结构完全统一）
    // ===========================
    while (error_k > 1e-5 && iter_cnt < iter_max)
    {
        // --------------------------
        // 几何量计算
        // --------------------------
        Eigen::MatrixXd delta_p = goodPos.rowwise() - pu_k.transpose();
        Eigen::MatrixXd delta_v = goodVel.rowwise() - vu_k.transpose();
        Eigen::VectorXd r_norm = delta_p.rowwise().norm();
        Eigen::MatrixXd los = delta_p.array().colwise() / r_norm.array();

        // --------------------------
        // 多普勒 H 矩阵（MATLAB原版）
        // --------------------------
        Eigen::MatrixXd vec2 = delta_v.array().colwise() / r_norm.array();
        Eigen::MatrixXd vec3 = vec3d_crossProduct(los, vec2);
        Eigen::MatrixXd h_pos = vec3d_crossProduct(los, vec3).array().colwise() * f0.array() / c;
        Eigen::MatrixXd h_vel = los.array().colwise() * (- f0.array()) / c;
        Eigen::VectorXd h_ddt = f0;

        H.block(0, 0, M, 3) = h_pos;
        H.block(0, 3, M, 3) = h_vel;
        H.col(6) = h_ddt;

        // --------------------------
        // 观测残差 b (Hz)
        // --------------------------
        Eigen::VectorXd proj_v = (delta_v.array() * los.array()).rowwise().sum();
        b = proj_v.array() * f0.array() / c
            + (ddtu_k - goodDdts.array()) * f0.array()
            + dopp.array();

        // --------------------------
        // 最小二乘（统一）
        // --------------------------
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(H);
        Eigen::VectorXd dx = -qr.solve(b);

        pu_k   += dx.segment<3>(0);
        vu_k   += dx.segment<3>(3);
        ddtu_k += dx(6);

        error_k = dx.cwiseAbs().maxCoeff();
        iter_cnt++;
    }

    // ===========================
    // 4. 结果输出
    // ===========================
    result.dop_pu       = pu_k;
    result.dop_vu       = vu_k;
    result.dop_ddtu     = ddtu_k;
    result.dop_residual = b;  // 单位：Hz
    result.dop_H        = H;
    result.dopValid     = true;

    result.dop_prns    = goodPrns;      // 自动保存
    result.dop_sys     = goodSys;       // 自动保存
    // 坐标转换
    result.dop_pu_lla = ecef2lla(result.dop_pu);
    result.dop_vu_enu = ecefVel2EnuVel(result.dop_vu, result.dop_pu_lla);

    if (printflag == 1) {
        printDoppPvtDebug(result);
    }

    return;
}


// 对应MATLAB vec3d_crossProduct
Eigen::MatrixXd UserPvtSolver::vec3d_crossProduct(const Eigen::MatrixXd& v1, const Eigen::MatrixXd& v2)
{
    int M = v1.rows();
    Eigen::MatrixXd res(M, 3);

    for (int i = 0; i < M; ++i) {
        double x1 = v1(i,0), y1 = v1(i,1), z1 = v1(i,2);
        double x2 = v2(i,0), y2 = v2(i,1), z2 = v2(i,2);
        res(i,0) = y1*z2 - z1*y2;
        res(i,1) = z1*x2 - x1*z2;
        res(i,2) = x1*y2 - y1*x2;
    }
    return res;
}

// 几何矩阵 H
Eigen::MatrixXd UserPvtSolver::calGeoMat(const Eigen::MatrixXd& ps, const Eigen::Vector3d& pu)
{
    int M = ps.rows();
    Eigen::MatrixXd H(M, 4);
    for (int i = 0; i < M; ++i) {
        Eigen::Vector3d diff = ps.row(i) - pu.transpose();
        double norm = diff.norm();
        H(i, 0) = -diff(0) / norm;
        H(i, 1) = -diff(1) / norm;
        H(i, 2) = -diff(2) / norm;
        H(i, 3) = 1.0;
    }
    return H;
}

// 地球自转改正
Eigen::VectorXd UserPvtSolver::rotCorrection(const Eigen::MatrixXd& ps, const Eigen::Vector3d& pu)
{
    int M = ps.rows();
    Eigen::VectorXd r(M);
    for (int i = 0; i < M; ++i) {
        Eigen::Vector3d diff = ps.row(i) - pu.transpose();
        double geom_range = diff.norm();
        double sagnac = consts::EARTH_ROTATION_RATE * (ps(i, 0) * pu(1) - ps(i, 1) * pu(0)) / consts::SPEED_OF_LIGHT;
        r(i) = geom_range + sagnac;
    }
    return r;
}

// 视线向量 LoS
Eigen::MatrixXd UserPvtSolver::calLoS(const Eigen::MatrixXd& ps, const Eigen::Vector3d& pu)
{
    int M = ps.rows();
    Eigen::MatrixXd LoS(M, 3);
    for (int i = 0; i < M; ++i) {
        Eigen::Vector3d diff = ps.row(i) - pu.transpose();
        LoS.row(i) = diff.normalized();
    }
    return LoS;
}

// ECEF -> LLA 坐标转换
Eigen::Vector3d UserPvtSolver::ecef2lla(const Eigen::Vector3d& ecef)
{
    double x = ecef(0), y = ecef(1), z = ecef(2);
    double lon = std::atan2(y, x);
    double p = std::sqrt(x*x + y*y);

    double lat = std::atan2(z * consts::SEMI_MAJOR_AXIS, p * consts::SEMI_MINOR_AXIS);
    const int iter_max = 10;
    double tol = 1e-12;

    for (int i = 0; i < iter_max; ++i) {
        double sin_lat = std::sin(lat);
        double N       = consts::SEMI_MAJOR_AXIS / std::sqrt(1.0 - consts::FIRST_ECCENTRICITY_SQUARED * sin_lat * sin_lat);
        double lat_new = std::atan2(z + consts::FIRST_ECCENTRICITY_SQUARED * N * sin_lat, p);
        if (std::abs(lat_new - lat) < tol) break;
        lat = lat_new;
    }

    double sin_lat = std::sin(lat);
    double N       = consts::SEMI_MAJOR_AXIS / std::sqrt(1.0 - consts::FIRST_ECCENTRICITY_SQUARED * sin_lat * sin_lat);
    double alt     = p / std::cos(lat) - N;

    return Eigen::Vector3d(lat * 180.0 / consts::PI, lon * 180.0 / consts::PI, alt);
}

// ============================================================================
// ECEF 速度 转 ENU 速度
// ============================================================================
Eigen::Vector3d UserPvtSolver::ecefVel2EnuVel(
    const Eigen::Vector3d& vu_ecef, 
    const Eigen::Vector3d& lla)
{
    double lat = lla(0) * consts::PI / 180.0;
    double lon = lla(1) * consts::PI / 180.0;

    double sLat = sin(lat);
    double cLat = cos(lat);
    double sLon = sin(lon);
    double cLon = cos(lon);

    // 标准 ECEF->ENU 旋转矩阵（对速度同样适用）
    Eigen::Matrix3d R;
    R << -sLon,        cLon,        0,
         -sLat*cLon,  -sLat*sLon,   cLat,
          cLat*cLon,   cLat*sLon,   sLat;

    return R * vu_ecef;
}

// ==========================
// 调试打印函数（新增）
// ==========================
void UserPvtSolver::printPesuPvtDebug(const UserPvtResult& res)
{
    if (!res.isValid) {
        std::cout << "\n[PVT Solver] 解算无效！卫星数量不足4颗" << std::endl;
        return;
    }

    std::cout << "\n=====================================" << std::endl;
    std::cout << "          PVT 解算结果（调试）" << std::endl;
    std::cout << "=====================================" << std::endl;

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "位置 ECEF (m):   " << res.pu.transpose() << std::endl;
    std::cout << "位置 LLA (deg,m): " << res.pu_lla(0) << "  " << res.pu_lla(1) << "  " << res.pu_lla(2) << std::endl;
    std::cout << "速度 ECEF (m/s): " << res.vu.transpose() << std::endl;
    std::cout << "速度 ENU  (m/s): " << res.vu_enu.transpose() << std::endl;
    std::cout << "钟差 (m):        " << res.dtu * 299792458.0 << "  (s): " << res.dtu << std::endl;
    std::cout << "钟漂 (m/s):      " << res.ddtu * 299792458.0 << "  (s/s): " << res.ddtu << std::endl;
    std::cout << "参与卫星数:      " << res.rhor.size() << std::endl;

    // 残差统计
    double maxRes = res.rhor.cwiseAbs().maxCoeff();
    double meanRes = res.rhor.cwiseAbs().mean();
    std::cout << "伪距残差 最大值: " << maxRes << " m" << std::endl;
    std::cout << "伪距残差 平均值: " << meanRes << " m" << std::endl;

    std::cout << "=====================================\n" << std::endl;
}

void UserPvtSolver::printDoppPvtDebug(const UserPvtResult& res)
{
    std::cout << "\n=====================================\n";
    std::cout << "        多普勒PVT解算结果\n";
    std::cout << "=====================================\n";

    if (res.dopValid) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "多普勒位置 ECEF: " << res.dop_pu.transpose() << "\n";
        std::cout << "多普勒位置 LLA: "  << res.dop_pu_lla.transpose() << "\n";
        std::cout << "多普勒速度 ECEF: " << res.dop_vu.transpose() << "\n";
        std::cout << "多普勒速度 ENU:  " << res.dop_vu_enu.transpose() << "\n";
        std::cout << "多普勒钟漂:      " << res.dop_ddtu << "\n";
        std::cout << "多普勒残差平方:  " << res.dop_residual.squaredNorm() << "\n";
    } else {
        std::cout << "多普勒PVT无效\n";
    }
}