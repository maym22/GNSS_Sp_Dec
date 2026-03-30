#pragma once
#ifndef OBS_READER_HPP
#define OBS_READER_HPP

#include "config/ConfigReader.hpp"  // 统一配置类
#include "pvt/Consts.hpp"   // 你工程里的常量定义：c、地球自转、椭球参数等
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <regex>
#include <stdexcept>

// 匹配Matlab的obs_gstruct：观测量核心结构体
struct ObsRawData {
    char sys;                // 系统类型（G=GPS, C=BDS）
    int PRN;                 // 卫星PRN号（BDS偏移32后）
    std::vector<int> Time;   // 原始时间数组 [年,月,日,时,分,秒]（对应Matlab rec_d）
    double ObsTime;          // TOW时间（秒，对应Matlab rec_t）
    double Fc;               // 载波频率（Hz）：GPS=1575.42e6, BDS=1561.098e6
    double Rho;              // 伪距（m，对应Matlab Rho）
    double AcPh;             // 载波相位（cycle，对应Matlab AcPh）
    double Fd;               // 多普勒频移（Hz，对应Matlab Fd）
    double CNR;              // 载噪比（dB-Hz，对应Matlab CNR）
};

// 观测量读取类（复用CfgReader统一配置）
class ObsReader {
private:
    const CfgReader& cfg_;               // 统一配置引用
    const double GPS_FREQ = 1575.42e6;   // GPS载波频率(Hz)
    const double BDS_FREQ = 1561.098e6;  // BDS载波频率(Hz)

    // 私有辅助函数
    int skipObsHeader(const std::vector<std::string>& flines);  // 跳过文件头
    std::vector<int> parseRecD(const std::string& line);        
    std::vector<double> parseRecD_t(const std::string& line);     // 解析时间数组rec_d
    double dt2WNTOW_Obs(const std::vector<int>& rec_d, char sys);   // 时间转TOW（匹配Matlab dateTool）
    bool parseObsLine(const std::string& line, const std::vector<int>& rec_d, const double ObsTimeRes, ObsRawData& obs_data); // 解析观测行

public:
    // 构造函数：传入统一配置
    ObsReader(const CfgReader& cfg) : cfg_(cfg) {}

    // 核心接口：读取观测量（完全匹配Matlab readObs）
    std::vector<std::vector<ObsRawData>> readObs();

    // 打印函数：格式化输出（便于和Matlab核对）
    void printObsData(const std::vector<std::vector<ObsRawData>>& obs_output) const;
};

#endif // OBS_READER_HPP