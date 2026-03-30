#include <iostream>
#include <config/ConfigReader.hpp>
#include <data_io/EphemerisReader.hpp>
#include <data_io/ObservationReader.hpp>
#include <antispoofing/AntiSpoofing.hpp>
#include <pvt/SatPvtSolver.hpp>

int main(int argc, char* argv[]) {

    try {
        // 1.1 初始化配置读取器（替换为你的yaml文件路径）
        CfgReader cfg_reader("/home/mym/ws_mym/Sf_ver0/config/config.yaml");
        std::cout << "[√] 配置文件读取成功！" << std::endl;
        // std::cout << "\n===== 配置读取结果 =====" << std::endl;
        // cfg_reader.printAllConfig();


        // 1.2 读取星历数据
        EphemerisReader eph_reader(cfg_reader);
        std::vector<EphemerisData> eph_data = eph_reader.readEph();
        std::cout << "[√] 星历数据读取成功！" << std::endl;
        // 验证读取结果
        // std::cout << "\n===== 星历读取结果 =====" << std::endl;
        // eph_reader.printEphemerisData(eph_data);

        // 1.3 读取观测量数据
        ObsReader obs_reader(cfg_reader);
        std::vector<std::vector<ObsRawData>> obs_data = obs_reader.readObs();
        std::cout << "[√] 观测量数据读取成功！" << std::endl;
        // 打印验证
        // std::cout << "\n===== 观测量读取结果 =====" << std::endl;
        // obs_reader.printObsData(obs_data);

        // 1.4: SatPVT解算（核心新增逻辑）
        SatPVTSolver pvt_solver(cfg_reader, eph_data);
        auto pvt_results = pvt_solver.solveAllPVT(obs_data);
        std::cout << "[√] SatPVT解算完成！" << std::endl;
        // std::cout << "\n===== SatPVT解算结果 =====" << std::endl;
        // pvt_solver.printPVTResult(pvt_results);

        // 2. 欺骗防御主流程
        AntiSpoofing as_detector(cfg_reader, eph_data, obs_data, pvt_results);
        auto spoof_output = as_detector.runAllDetection();
        as_detector.printSpoofResult(spoof_output);

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "程序失败" << e.what() << std::endl;
        return 1;
    }
}