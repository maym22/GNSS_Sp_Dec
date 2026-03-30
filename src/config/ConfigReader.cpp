#include <config/ConfigReader.hpp>

// 构造函数：调用加载配置方法
CfgReader::CfgReader(const std::string& config_path) {
    loadConfig(config_path);
}

// 核心：解析yaml配置文件
void CfgReader::loadConfig(const std::string& config_path) {
    try {

        YAML::Node config = YAML::LoadFile(config_path);
        if (!config) {
            throw std::runtime_error("配置文件为空或格式错误：" + config_path);
        }

        // ========== 解析卫星星座配置 ==========
        if (config["satellite"]) {
            satellite_.useGPS = config["satellite"]["useGPS"].as<int>();
            satellite_.useBDS = config["satellite"]["useBDS"].as<int>();
            satellite_.sat_max_num = satellite_.useGPS*32 + satellite_.useBDS*30;
        } else {
            throw std::runtime_error("缺失 'satellite' 配置节点");
        }

        // ========== 解析频率配置 ==========
        if (config["frequency"]) {
            frequency_.Fif = config["frequency"]["Fif"].as<double>();
            frequency_.Fs = config["frequency"]["Fs"].as<double>();
            frequency_.T_samp = config["frequency"]["T_samp"].as<double>();
        } else {
            throw std::runtime_error("缺失 'frequency' 配置节点");
        }

        // ========== 解析时间配置 ==========
        if (config["datetime"]) {
            datetime_.AprxTime.start_gps_week = config["datetime"]["AprxTime"]["start_gps_week"].as<std::string>();
            datetime_.AprxTime.ephemeris_read = config["datetime"]["AprxTime"]["ephemeris_read"].as<std::string>();
            datetime_.AprxTime.ref_observation = config["datetime"]["AprxTime"]["ref_observation"].as<std::string>();
            datetime_.DurMax = config["datetime"]["DurMax"].as<std::string>();
        } else {
            throw std::runtime_error("缺失 'datetime' 配置节点");
        }


        // ========== 解析观测数据配置 ==========
        if (config["observation"]) {
            observation_.data_path = config["observation"]["data_path"].as<std::string>();
            observation_.ObsMax = config["observation"]["ObsMax"].as<int>();
            observation_.ObsPer = config["observation"]["ObsPer"].as<int>();
        } else {
            throw std::runtime_error("缺失 'observation' 配置节点");
        }

        // ========== 解析星历数据配置 ==========
        if (config["ephemeris"]) {
            ephemeris_.data_path = config["ephemeris"]["data_path"].as<std::string>();
        } else {
            throw std::runtime_error("缺失 'ephemeris' 配置节点");
        }


        // ========== 解析反欺骗配置 ==========
        if (config["anti_spoofing"]) {
            YAML::Node as_node = config["anti_spoofing"];
            
            // M1: 星历规范性检查
            anti_spoofing_.chnEphemNormativity.flag   = as_node["chnEphemNormativity"]["flag"].as<bool>();
            anti_spoofing_.chnEphemNormativity.useGPS = satellite_.useGPS;
            anti_spoofing_.chnEphemNormativity.useBDS = satellite_.useBDS;


            // M2: 卫星可见性检查
            anti_spoofing_.chnSatVisibility.flag = as_node["chnSatVisibility"]["flag"].as<bool>();
            anti_spoofing_.chnSatVisibility.Thr = as_node["chnSatVisibility"]["Thr"].as<int>();
            anti_spoofing_.chnSatVisibility.rec_pos_lla = as_node["chnSatVisibility"]["rec_pos_lla"].as<std::vector<double>>();

            // M3: C/N0异常检测
            anti_spoofing_.chnCn0.flag = as_node["chnCn0"]["flag"].as<bool>();
            anti_spoofing_.chnCn0.WindowSize =  as_node["chnCn0"]["WindowSize"].as<int>();
            anti_spoofing_.chnCn0.Thr = as_node["chnCn0"]["Thr"].as<int>();

            // M4: 多普勒一致性检查
            anti_spoofing_.chnDoppCon.flag = as_node["chnDoppCon"]["flag"].as<bool>();
            anti_spoofing_.chnDoppCon.Thr = as_node["chnDoppCon"]["Thr"].as<int>();
            anti_spoofing_.chnDoppCon.Fif = as_node["chnDoppCon"]["Fif"].as<double>();

            // M5: 单星钟漂变化检查
            anti_spoofing_.chnClkDft.flag = as_node["chnClkDft"]["flag"].as<bool>();
            anti_spoofing_.chnClkDft.Thr = as_node["chnClkDft"]["Thr"].as<int>();

            // M6: 信号传播时间有效性检查
            anti_spoofing_.chnRho.flag = as_node["chnRho"]["flag"].as<bool>();
            anti_spoofing_.chnRho.Thr = as_node["chnRho"]["Thr"].as<std::vector<double>>();
 

            // M7: C/N0相关性检查
            anti_spoofing_.mchnCn0Corr.flag          = as_node["mchnCn0Corr"]["flag"].as<bool>();
            anti_spoofing_.mchnCn0Corr.WindowSize    = as_node["mchnCn0Corr"]["WindowSize"].as<int>();
            anti_spoofing_.mchnCn0Corr.MinWindowSize = as_node["mchnCn0Corr"]["MinWindowSize"].as<int>();
            anti_spoofing_.mchnCn0Corr.Pfa           = as_node["mchnCn0Corr"]["Pfa"].as<double>();

            // M8: 多普勒单差检查
            anti_spoofing_.mchnDoppSD.flag          = as_node["mchnDoppSD"]["flag"].as<bool>();
            anti_spoofing_.mchnDoppSD.MinWindowSize = as_node["mchnDoppSD"]["MinWindowSize"].as<int>();
            anti_spoofing_.mchnDoppSD.WindowSize    = as_node["mchnDoppSD"]["WindowSize"].as<int>();
            anti_spoofing_.mchnDoppSD.Thr           = as_node["mchnDoppSD"]["Thr"].as<double>();

            // M9: 标准RAIM
            anti_spoofing_.AntRaim.flag   = as_node["AntRaim"]["flag"].as<bool>();
            anti_spoofing_.AntRaim.Pfa    = as_node["AntRaim"]["Pfa"].as<double>();
            anti_spoofing_.AntRaim.Sigma  = as_node["AntRaim"]["Sigma"].as<double>();

            // M10: 生成式RAIM
            anti_spoofing_.AntGRaim.flag  = as_node["AntGRaim"]["flag"].as<bool>();
            anti_spoofing_.AntGRaim.Pfa   = as_node["AntGRaim"]["Pfa"].as<double>();
            anti_spoofing_.AntGRaim.Sigma = as_node["AntGRaim"]["Sigma"].as<double>();
            anti_spoofing_.AntGRaim.Mode  = as_node["AntGRaim"]["Mode"].as<int>();

            // M11: 多普勒RAIM
            anti_spoofing_.AntDRaim.flag  = as_node["AntDRaim"]["flag"].as<bool>();
            anti_spoofing_.AntDRaim.Pfa   = as_node["AntDRaim"]["Pfa"].as<double>();
            anti_spoofing_.AntDRaim.Sigma = as_node["AntDRaim"]["Sigma"].as<double>();

            // M12: PVT验证
            anti_spoofing_.AntPVT.flag = as_node["AntPVT"]["flag"].as<bool>();
            anti_spoofing_.AntPVT.MinHeight = as_node["AntPVT"]["MinHeight"].as<double>();
            anti_spoofing_.AntPVT.MaxHeight = as_node["AntPVT"]["MaxHeight"].as<double>();
            anti_spoofing_.AntPVT.MaxVelocity = as_node["AntPVT"]["MaxVelocity"].as<double>();
            anti_spoofing_.AntPVT.MaxClockDrift = as_node["AntPVT"]["MaxClockDrift"].as<double>();
            anti_spoofing_.AntPVT.MaxPositionChange = as_node["AntPVT"]["MaxPositionChange"].as<double>();

            // M13: 速度一致性
            anti_spoofing_.AntVelCon.flag       = as_node["AntVelCon"]["flag"].as<bool>();
            anti_spoofing_.AntVelCon.WindowSize = as_node["AntVelCon"]["WindowSize"].as<int>();
            anti_spoofing_.AntVelCon.Pfa        = as_node["AntVelCon"]["Pfa"].as<double>();
            anti_spoofing_.AntVelCon.Sigma      = as_node["AntVelCon"]["Sigma"].as<double>();
            
            // M14: IMU验证
            anti_spoofing_.IMUVal.flag = as_node["IMUVal"]["flag"].as<bool>();
            anti_spoofing_.IMUVal.DeltaVelThr = as_node["IMUVal"]["DeltaVelThr"].as<int>();
            anti_spoofing_.IMUVal.ObsCnt = as_node["IMUVal"]["ObsCnt"].as<int>();

            // M15: PVT交叉验证
            anti_spoofing_.MantPVT.flag = as_node["MantPVT"]["flag"].as<bool>();
            anti_spoofing_.MantPVT.PosThr = as_node["MantPVT"]["PosThr"].as<int>();
            anti_spoofing_.MantPVT.TimThr = as_node["MantPVT"]["TimThr"].as<int>();
            anti_spoofing_.MantPVT.VelThr = as_node["MantPVT"]["VelThr"].as<int>();
            anti_spoofing_.MantPVT.ObsCnt = as_node["MantPVT"]["ObsCnt"].as<int>();

            // M16: CPDD
            anti_spoofing_.MantCPDD.flag = as_node["MantCPDD"]["flag"].as<bool>();
            anti_spoofing_.MantCPDD.Thr = as_node["MantCPDD"]["Thr"].as<double>();
            anti_spoofing_.MantCPDD.ObsCnt = as_node["MantCPDD"]["ObsCnt"].as<int>();
        
                      
            std::cout << "===== 卫星星座配置 =====" << std::endl;

        } else {
            throw std::runtime_error("缺失 'anti_spoofing' 配置节点");
        }

    } catch (const YAML::BadFile& e) {
        throw std::runtime_error("无法打开配置文件：" + config_path + "，错误：" + e.what());
    } catch (const YAML::RepresentationException& e) {
        throw std::runtime_error("配置参数缺失/格式错误：" + std::string(e.what()));
    } catch (const std::exception& e) {
        throw std::runtime_error("配置解析失败：" + std::string(e.what()));
    }
}

// 调试用：打印所有配置（方便验证）
void CfgReader::printAllConfig() const {
    std::cout << "===== 卫星星座配置 =====" << std::endl;
    std::cout << "useGPS: " << satellite_.useGPS << std::endl;
    std::cout << "useBDS: " << satellite_.useBDS << std::endl;
    std::cout << "sat_max_num: " << satellite_.sat_max_num << std::endl;

    std::cout << "\n===== 频率配置 (kHz) =====" << std::endl;
    std::cout << "Fif: " << frequency_.Fif << std::endl;
    std::cout << "Fs: " << frequency_.Fs << std::endl;

    std::cout << "\n===== 时间配置 =====" << std::endl;
    std::cout << "start_gps_week: " << datetime_.AprxTime.start_gps_week << std::endl;
    std::cout << "ephemeris_read: " << datetime_.AprxTime.ephemeris_read << std::endl;
    std::cout << "ref_observation: " << datetime_.AprxTime.ref_observation << std::endl;
    std::cout << "DurMax: " << datetime_.DurMax << std::endl;

    std::cout << "\n===== 观测数据配置 =====" << std::endl;
    std::cout << "data_path: " << observation_.data_path << std::endl;
    std::cout << "ObsMax: " << observation_.ObsMax << std::endl;
    std::cout << "ObsPer: " << observation_.ObsPer << " ms" << std::endl;

    std::cout << "\n===== 星历数据配置 =====" << std::endl;
    std::cout << "data_path: " << ephemeris_.data_path << std::endl;
 
    std::cout << "\n===== 反欺骗配置 =====" << std::endl;
    // M1
    std::cout << "M1-星历规范性: useGPS=" << anti_spoofing_.chnEphemNormativity.useGPS 
              << ", useBDS=" << anti_spoofing_.chnEphemNormativity.useBDS << std::endl;
    // M2
    std::cout << "M2-卫星可见性: Thr=" << anti_spoofing_.chnSatVisibility.Thr << " deg, rec_pos_lla=[";
    for (size_t i=0; i<anti_spoofing_.chnSatVisibility.rec_pos_lla.size(); ++i) {
        std::cout << anti_spoofing_.chnSatVisibility.rec_pos_lla[i];
        if (i < anti_spoofing_.chnSatVisibility.rec_pos_lla.size()-1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    // M3
    std::cout << "M3-C/N0异常: Thr=" << anti_spoofing_.chnCn0.Thr << std::endl;
    // M6
    std::cout << "M6-信号传播时间: Thr=[" << anti_spoofing_.chnRho.Thr[0] << ", " 
              << anti_spoofing_.chnRho.Thr[1] << "] s" << std::endl;
    // M13
    std::cout << "M13-PVT验证: PosThr=" << anti_spoofing_.AntPVT.MaxPositionChange << " m, DftThr=" 
              << anti_spoofing_.AntPVT.MaxClockDrift << ", Height=" << anti_spoofing_.AntPVT.MinHeight << " m" << std::endl;
    // 其他子模块可按需补充打印
}