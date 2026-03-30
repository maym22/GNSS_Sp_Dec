#pragma once
#ifndef EPHEMERIS_READER_HPP
#define EPHEMERIS_READER_HPP

#include "config/ConfigReader.hpp"  // 引入已定义的CfgReader类
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <regex>
#include <cstring>

// 匹配Matlab的eph_struct结构
struct EphemerisData {
    // 基础信息
    char sys;                // 系统(G/C)
    int  PRN;                 // 卫星PRN号
    // Toc相关
    std::tm Toc_cal;         // Toc日历时间
    int Toc_WN;              // Toc周数
    double Toc_TOW;          // Toc时间Of周
    // 钟差参数
    double af0;
    double af1;
    double af2;
    // 轨道参数
    double IODE;
    double Crs;
    double Delta_n;
    double M_0;
    double Cuc;
    double e;
    double Cus;
    double sqrt_a;
    double Toes;
    double Toe;
    double Cic;
    double Omega_0;
    double Cis;
    double i_0;
    double Crc;
    double omega;
    double Omega_dot;
    double I_dot;
    int WN;
    double Health;
    double TGD;

    // 构造函数初始化
    EphemerisData() : sys('\0'), PRN(0), Toc_WN(0), Toc_TOW(0.0),
                      af0(0.0), af1(0.0), af2(0.0), IODE(0.0), Crs(0.0),
                      Delta_n(0.0), M_0(0.0), Cuc(0.0), e(0.0), Cus(0.0),
                      sqrt_a(0.0), Toes(0.0), Toe(0.0), Cic(0.0), Omega_0(0.0),
                      Cis(0.0), i_0(0.0), Crc(0.0), omega(0.0), Omega_dot(0.0),
                      I_dot(0.0), WN(0), Health(0.0), TGD(0.0) {
        memset(&Toc_cal, 0, sizeof(std::tm));
    }
};

// 星历读取类（复用CfgReader配置）
class EphemerisReader {
private:
    // 持有CfgReader的常量引用（避免拷贝，保证配置统一）
    const CfgReader& cfg_;
    // 辅助参数（GPST转BDT偏移秒数，可后续放入yaml配置）
    const double GPS_to_BDT_offset = 14.0;

    // 辅助函数：跳过文件头直到END OF HEADER
    int skipHeader(std::vector<std::string>& flines);
    // 辅助函数：计算时间差（秒）
    double calculateTimeDiff(const std::tm& t1, const std::tm& t2);
    // 辅助函数：字符串替换（D->e）
    std::string replaceDWithE(const std::string& str);
    // 辅助函数：解析单行数值
    std::vector<double> parseLine(const std::string& line);
    // 辅助函数：转换时间到WN和TOW
    void dt2WNTOW(const std::tm& dt, char sys, int& WN, double& TOW);
    // 辅助函数：将yaml中的时间字符串（如"2026-01-13T18:00:00"）转换为std::tm
    bool parseDatetimeStr(const std::string& datetime_str, std::tm& out_tm);

public:
    // 构造函数：传入CfgReader引用（核心！复用已有配置）
    EphemerisReader(const CfgReader& cfg) : cfg_(cfg) {}
    // 读取星历数据核心接口
    std::vector<EphemerisData> readEph();
    void printEphemerisData(const std::vector<EphemerisData>& eph_data) const;
};

#endif // EPHEMERIS_READER_HPP