#pragma once
#ifndef CFG_READ_HPP
#define CFG_READ_HPP

#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>

// 配置读取类：完整解析你的yaml配置文件
class CfgReader {

public:
    // ========== 配置参数结构体（按yaml层级封装）==========
    // 卫星星座配置
    struct SatelliteConfig {
        int useGPS;
        int useBDS;
        int sat_max_num;
    } satellite_;

    // 频率配置
    struct FrequencyConfig {
        double Fif;
        double Fs;
        double T_samp; 
    } frequency_;

    // 时间配置
    struct AprxTimeConfig {
        std::string start_gps_week;
        std::string ephemeris_read;
        std::string ref_observation;
    };
    struct DatetimeConfig {
        AprxTimeConfig AprxTime;
        std::string DurMax;
    } datetime_;

    // 观测数据配置
    struct ObservationConfig {
        std::string data_path;
        int ObsMax;
        int ObsPer;
    } observation_;

    // 星历数据配置
    struct EphemerisConfig {
        std::string data_path;
    } ephemeris_;

    // 反欺骗配置 - 子模块
    struct ChnEphemNormativityConfig {
        bool flag;
        int useGPS;
        int useBDS;
    };
    struct ChnSatVisibilityConfig {
        bool flag;
        int Thr;
        std::vector<double> rec_pos_lla;
    };
    struct ChnCn0Config {
        bool flag;
        int WindowSize;
        int Thr;
    };
    struct ChnDoppConConfig {
        bool flag;
        int WindowSize;
        int Thr;
        double Fif;
    };
    struct ChnClkDftConfig {
        bool flag;
        int WindowSize;
        int Thr;
    };
    struct ChnRhoConfig {
        bool flag;
        int WindowSize;
        std::vector<double> Thr;
    };
    struct MchnCn0CorrConfig {
        bool   flag;
        double Pfa;
        int    WindowSize;
        int    MinWindowSize;
    };
    struct MchnDoppSDConfig {
        bool flag;
        int MinWindowSize;
        int WindowSize;
        double Thr;
    };
    struct AntRaimConfig {
        bool   flag;
        double Pfa;
        double Sigma;
    };
    struct AntGRaimConfig {
        bool   flag;
        int    Mode;
        double Pfa;
        double Sigma;
    };
    struct AntDRaimConfig {
        bool   flag;
        double Pfa;
        double Sigma;
    };
    struct AntPVTConfig {
        bool flag;
        double MinHeight;
        double MaxHeight;
        double MaxVelocity;
        double MaxClockDrift;
        double MaxPositionChange;
    };
    struct AntVelConConfig {
        bool flag;
        int WindowSize;
        double Pfa;
        double Sigma;
    };
    
    struct IMUValConfig {
        bool flag;
        int DeltaVelThr;
        int ObsCnt;
    };
    struct MantPVTConfig {
        bool flag;
        int PosThr;
        int TimThr;
        int VelThr;
        int ObsCnt;
    };
    struct MantCPDDConfig {
        bool flag;
        double Thr;
        int ObsCnt;
    };

    // 反欺骗总配置
    struct AntiSpoofingConfig {
        ChnEphemNormativityConfig chnEphemNormativity;
        ChnSatVisibilityConfig chnSatVisibility;
        ChnCn0Config chnCn0;
        ChnDoppConConfig chnDoppCon;
        ChnClkDftConfig chnClkDft;
        ChnRhoConfig chnRho;
        MchnCn0CorrConfig mchnCn0Corr;
        MchnDoppSDConfig mchnDoppSD;
        AntRaimConfig AntRaim;
        AntGRaimConfig AntGRaim;
        AntDRaimConfig AntDRaim;
        AntPVTConfig AntPVT;
        AntVelConConfig AntVelCon;
        IMUValConfig IMUVal;
        MantPVTConfig MantPVT;
        MantCPDDConfig MantCPDD;
    } anti_spoofing_;

    // 私有方法：加载并解析yaml配置
    void loadConfig(const std::string& config_path);

    // 构造函数：传入yaml文件路径，自动解析配置
    CfgReader(const std::string& config_path);

    // ========== 获取配置的接口（GETTER）==========
    // 卫星配置
    const SatelliteConfig& getSatelliteConfig() const { return satellite_; }
    // 频率配置
    const FrequencyConfig& getFrequencyConfig() const { return frequency_; }
    // 时间配置
    const DatetimeConfig& getDatetimeConfig() const { return datetime_; }
    // 观测数据配置
    const ObservationConfig& getObservationConfig() const { return observation_; }
    // 星历数据配置
    const EphemerisConfig& getEphemerisConfig() const { return ephemeris_; }
    // 反欺骗配置
    const AntiSpoofingConfig& getAntiSpoofingConfig() const { return anti_spoofing_; }

    // 调试用：打印所有配置（验证读取是否正确）
    void printAllConfig() const;
};

#endif // CFG_READ_HPP