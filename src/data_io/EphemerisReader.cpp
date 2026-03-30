#include "data_io/EphemerisReader.hpp"
#include <sstream>
#include <regex>
#include <iomanip>
#include <stdexcept>
#include <iostream>

// 跳过文件头
int EphemerisReader::skipHeader(std::vector<std::string>& flines) {
    int line_idx = 0;
    for (; line_idx < flines.size(); ++line_idx) {
        if (flines[line_idx].find("END OF HEADER") != std::string::npos) {
            break;
        }
    }
    return line_idx + 1; // 跳过END OF HEADER行
}

// 替换字符串中的D为e（科学计数法兼容）
std::string EphemerisReader::replaceDWithE(const std::string& str) {
    std::string res = str;
    std::replace(res.begin(), res.end(), 'D', 'e');
    return res;
}

// 计算两个tm时间的差值（秒）
double EphemerisReader::calculateTimeDiff(const std::tm& t1, const std::tm& t2) {
    time_t time1 = mktime(const_cast<std::tm*>(&t1));
    time_t time2 = mktime(const_cast<std::tm*>(&t2));
    if (time1 == -1 || time2 == -1) {
        throw std::runtime_error("时间转换失败");
    }
    return std::difftime(time1, time2);
}

// 解析单行数值
std::vector<double> EphemerisReader::parseLine(const std::string& line) {
    std::vector<double> res;
    std::istringstream iss(line);
    double val;
    while (iss >> val) {
        res.push_back(val);
    }
    return res;
}

// 时间字符串转std::tm（解析"2026-01-13 18:00:00"格式）
bool EphemerisReader::parseDatetimeStr(const std::string& datetime_str, std::tm& out_tm) {
    memset(&out_tm, 0, sizeof(std::tm));
    std::istringstream ss(datetime_str);
    ss >> std::get_time(&out_tm, "%Y-%m-%d %H:%M:%S");
    return !ss.fail();
}

// 时间转换为WN和TOW（简化实现，可根据GPST/BDT规则完善）
void EphemerisReader::dt2WNTOW(const std::tm& dt, char sys, int& WN, double& TOW) {
    time_t t = mktime(const_cast<std::tm*>(&dt));
    if (t == -1) {
        WN = 0;
        TOW = 0.0;
        return;
    }

    // 考虑闰秒的影响
    t = t + 18; 

    // GPST周数计算（起始时间：1980-01-06 00:00:00）
    tm epoch_tm = {0};
    if (sys == 'G') { // GPS起始：1980-01-06 00:00:00 GPST
        epoch_tm.tm_year = 80;
        epoch_tm.tm_mon  = 0;
        epoch_tm.tm_mday = 6;
    } else if (sys == 'C') { // BDS起始：2006-01-01 00:00:00 BDT
        epoch_tm.tm_year = 106;
        epoch_tm.tm_mon  = 0;
        epoch_tm.tm_mday = 1;
    }
    const time_t gpst_epoch = mktime(&epoch_tm);

    double sec_since_epoch = std::difftime(t, gpst_epoch);
    WN = static_cast<int>(sec_since_epoch / (7 * 24 * 3600));
    TOW = fmod(sec_since_epoch, 7 * 24 * 3600);
    // BDT偏移修正
    if (sys == 'C') {
        TOW -= GPS_to_BDT_offset;
    }
}

// 核心：读取星历数据（复用CfgReader配置）
std::vector<EphemerisData> EphemerisReader::readEph() {
    // 1. 从CfgReader获取配置参数
    const auto& sat_cfg      = cfg_.getSatelliteConfig();       // 卫星配置（useGPS/useBDS/sat_max_num）
    const auto& datetime_cfg = cfg_.getDatetimeConfig();        // 时间配置（ephemeris_read参考时间）
    const auto& eph_cfg      = cfg_.getEphemerisConfig();       // 星历路径配置（data_path）
    
    // 星历文件完整路径（拼接data_dir和默认文件名，可根据suffix扩展）
    std::string eph_file_path = eph_cfg.data_path;

    // DurMax转换（从"03:59:59"字符串转秒数）
    std::tm dur_tm = {0};
    std::istringstream dur_ss(datetime_cfg.DurMax);
    dur_ss >> std::get_time(&dur_tm, "%H:%M:%S");
    double DurMax = dur_tm.tm_hour * 3600 + dur_tm.tm_min * 60 + dur_tm.tm_sec;

    // 2. 初始化星历数组和时间差数组
    std::vector<EphemerisData> eph_output(sat_cfg.sat_max_num);
    std::vector<double> delta_t(sat_cfg.sat_max_num, DurMax); // 初始化为最大容差

    // 3. 解析参考时间（ephemeris_read）
    std::tm ref_time;
    if (!parseDatetimeStr(datetime_cfg.AprxTime.ephemeris_read, ref_time)) {
        throw std::runtime_error("参考时间解析失败：" + datetime_cfg.AprxTime.ephemeris_read);
    }


    // 4. 读取星历文件
    std::vector<std::string> flines;
    std::ifstream file(eph_file_path);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开星历文件：" + eph_file_path);
    }
    std::string line;
    while (std::getline(file, line)) {
        flines.push_back(line);
    }
    file.close();

    // 5. 跳过文件头
    uint line_idx = skipHeader(flines);

    // 6. 遍历解析星历数据
    std::string use_sys;
    if (sat_cfg.useGPS) use_sys += "G";
    if (sat_cfg.useBDS) use_sys += "C";

    while (line_idx < flines.size()) {
        std::string fline = flines[line_idx];
        line_idx++;

        // 过滤非目标系统的星历
        if (fline.empty() || use_sys.find(fline[0]) == std::string::npos) {
            continue;
        }

        EphemerisData eph_tmp;
        eph_tmp.sys = fline[0];
        fline = replaceDWithE(fline);
        
        // 正则替换：数字和正负号间加空格（匹配Matlab regexprep逻辑）
        std::regex reg(R"((\d)([-+]))");
        fline = std::regex_replace(fline, reg, "$1 $2");

        // 解析第一行核心数据
        std::vector<double> str = parseLine(fline.substr(1));
        if (str.size() < 10) {
            continue;
        }

        // 计算PRN（BDS需偏移GPS的32）
        eph_tmp.PRN = static_cast<int>(str[0]);
        if (eph_tmp.sys == 'C' && sat_cfg.useGPS) {
            eph_tmp.PRN += 32;
        }
        // 边界检查
        if (eph_tmp.PRN <= 0 || eph_tmp.PRN > sat_cfg.sat_max_num) {
            continue;
        }

        if (str.size() < 7) { // 确保有PRN+6位时间字段
            continue;
        }

        int year      = static_cast<int>(str[1]);
        int month     = static_cast<int>(str[2]);
        int day       = static_cast<int>(str[3]);
        int hour      = static_cast<int>(str[4]);
        int minute    = static_cast<int>(str[5]);
        int second    = static_cast<int>(str[6]); // 秒可能含小数，取整或保留均可

        // 2. 拼接为 "YYYY-MM-DD HH:MM:SS" 格式的字符串（匹配parseDatetimeStr的解析格式）
        char toc_str_buf[32];
        // 若秒是整数：直接拼接
        sprintf(toc_str_buf, "%04d-%02d-%02d %02d:%02d:%02d", 
                year, month, day, hour, minute, second);
        
        std::string toc_datetime_str = toc_str_buf;
        // 解析Toc时间并计算与参考时间的差值
        if (!parseDatetimeStr(toc_datetime_str, eph_tmp.Toc_cal)) {
            continue;
        }

        // 转换Toc到WN和TOW
        dt2WNTOW(eph_tmp.Toc_cal, eph_tmp.sys, eph_tmp.Toc_WN, eph_tmp.Toc_TOW);

        // BDT时间偏移修正
        if (eph_tmp.sys == 'C') {
            time_t toc_time = mktime(&eph_tmp.Toc_cal);
            toc_time -= static_cast<time_t>(GPS_to_BDT_offset);
            eph_tmp.Toc_cal = *localtime(&toc_time);
        }

        // 筛选最近星历（时间差更小则更新）
        double diff = fabs(calculateTimeDiff(eph_tmp.Toc_cal, ref_time));
        if (diff >= delta_t[eph_tmp.PRN - 1]) {
            continue;
        }
        delta_t[eph_tmp.PRN - 1] = diff;

        // 填充钟差参数
        eph_tmp.af0 = str[7];
        eph_tmp.af1 = str[8];
        eph_tmp.af2 = str[9];

        // 解析第二行（轨道参数1）
        if (line_idx >= flines.size()) break;
        std::string line2 = replaceDWithE(flines[line_idx++]);
        std::vector<double> str2 = parseLine(line2);
        if (str2.size() >= 4) {
            eph_tmp.IODE = str2[0];
            eph_tmp.Crs = str2[1];
            eph_tmp.Delta_n = str2[2];
            eph_tmp.M_0 = str2[3];
        }

        // 解析第三行（轨道参数2）
        if (line_idx >= flines.size()) break;
        std::string line3 = replaceDWithE(flines[line_idx++]);
        std::vector<double> str3 = parseLine(line3);
        if (str3.size() >= 4) {
            eph_tmp.Cuc = str3[0];
            eph_tmp.e = str3[1];
            eph_tmp.Cus = str3[2];
            eph_tmp.sqrt_a = str3[3];
        }

        // 解析第四行（轨道参数3）
        if (line_idx >= flines.size()) break;
        std::string line4 = replaceDWithE(flines[line_idx++]);
        std::vector<double> str4 = parseLine(line4);
        if (str4.size() >= 4) {
            eph_tmp.Toes = str4[0] - (eph_tmp.sys == 'C' ? GPS_to_BDT_offset : 0);
            eph_tmp.Toe = eph_tmp.Toes;
            eph_tmp.Cic = str4[1];
            eph_tmp.Omega_0 = str4[2];
            eph_tmp.Cis = str4[3];
        }

        // 解析第五行（轨道参数4）
        if (line_idx >= flines.size()) break;
        std::string line5 = replaceDWithE(flines[line_idx++]);
        std::vector<double> str5 = parseLine(line5);
        if (str5.size() >= 4) {
            eph_tmp.i_0 = str5[0];
            eph_tmp.Crc = str5[1];
            eph_tmp.omega = str5[2];
            eph_tmp.Omega_dot = str5[3];
        }

        // 解析第六行（轨道参数5）
        if (line_idx >= flines.size()) break;
        std::string line6 = replaceDWithE(flines[line_idx++]);
        std::vector<double> str6 = parseLine(line6);
        if (str6.size() >= 3) {
            eph_tmp.I_dot = str6[0];
            eph_tmp.WN = static_cast<int>(str6[2]);
        }

        // 解析第七行（健康状态+TGD）
        if (line_idx >= flines.size()) break;
        std::string line7 = replaceDWithE(flines[line_idx++]);
        std::vector<double> str7 = parseLine(line7);
        if (str7.size() >= 3) {
            eph_tmp.Health = str7[1];
            eph_tmp.TGD = str7[2];
        }

        // 保存有效星历
        eph_output[eph_tmp.PRN - 1] = eph_tmp;
    }

    // 7. 统计有效星历数量
    int valid_count = 0;
    for (const auto& eph : eph_output) {
        if (eph.sys != '\0') {
            valid_count++;
        }
    }

    std::cout << "[星历读取] 有效星历数量：" << valid_count << "/" << sat_cfg.sat_max_num << std::endl;

    return eph_output;
}

// 完整打印星历数据（按字段顺序，适配Matlab核对）
void EphemerisReader::printEphemerisData(const std::vector<EphemerisData>& eph_data) const {
    std::cout << "\n========================================" << std::endl;
    std::cout << "              卫星星历解析结果              " << std::endl;
    std::cout << "========================================\n" << std::endl;

    int valid_count = 0;
    for (size_t i = 0; i < eph_data.size(); ++i) {
        const auto& eph = eph_data[i];
        // 仅打印有效星历
        if (eph.sys == '\0') {
            continue;
        }
        valid_count++;

        // 打印星历头部信息
        std::cout << "---------- 卫星 " << valid_count << " (PRN: " << eph.PRN << ") ----------" << std::endl;
        std::cout << "基础信息：" << std::endl;
        std::cout << "  系统类型: " << eph.sys << " (G=GPS, C=BDS)" << std::endl;
        std::cout << "  PRN编号: " << eph.PRN << std::endl;
        std::cout << "  健康状态: " << eph.Health << std::endl;

        // 钟差参数
        std::cout << "钟差参数：" << std::endl;
        std::cout << "  af0 (s): " << std::fixed << std::setprecision(15) << eph.af0 << std::endl;
        std::cout << "  af1 (s/s): " << std::fixed << std::setprecision(15) << eph.af1 << std::endl;
        std::cout << "  af2 (s/s²): " << std::fixed << std::setprecision(15) << eph.af2 << std::endl;

        // Toc相关
        std::cout << "Toc时间参数：" << std::endl;
        std::cout << "  Toc周数 (WN): " << eph.Toc_WN << std::endl;
        std::cout << "  Toc时间Of周 (TOW, s): " << std::fixed << std::setprecision(3) << eph.Toc_TOW << std::endl;
        // 格式化打印Toc日历时间
        char toc_time_buf[64];
        strftime(toc_time_buf, sizeof(toc_time_buf), "%Y-%m-%d %H:%M:%S", &eph.Toc_cal);
        std::cout << "  Toc日历时间: " << toc_time_buf << std::endl;

        // 轨道参数
        std::cout << "轨道参数：" << std::endl;
        std::cout << "  IODE: " << std::fixed << std::setprecision(6) << eph.IODE << std::endl;
        std::cout << "  Crs (m): " << std::fixed << std::setprecision(3) << eph.Crs << std::endl;
        std::cout << "  Delta_n (rad/s): " << std::fixed << std::setprecision(15) << eph.Delta_n << std::endl;
        std::cout << "  M_0 (rad): " << std::fixed << std::setprecision(15) << eph.M_0 << std::endl;
        std::cout << "  Cuc (rad): " << std::fixed << std::setprecision(15) << eph.Cuc << std::endl;
        std::cout << "  偏心率 e: " << std::fixed << std::setprecision(15) << eph.e << std::endl;
        std::cout << "  Cus (rad): " << std::fixed << std::setprecision(15) << eph.Cus << std::endl;
        std::cout << "  sqrt(a) (m^0.5): " << std::fixed << std::setprecision(3) << eph.sqrt_a << std::endl;
        std::cout << "  半长轴 a (m): " << std::fixed << std::setprecision(3) << pow(eph.sqrt_a, 2) << std::endl; // 计算a
        std::cout << "  Toes (s): " << std::fixed << std::setprecision(3) << eph.Toes << std::endl;
        std::cout << "  Toe (s): " << std::fixed << std::setprecision(3) << eph.Toe << std::endl;
        std::cout << "  Cic (rad): " << std::fixed << std::setprecision(15) << eph.Cic << std::endl;
        std::cout << "  Omega_0 (rad): " << std::fixed << std::setprecision(15) << eph.Omega_0 << std::endl;
        std::cout << "  Cis (rad): " << std::fixed << std::setprecision(15) << eph.Cis << std::endl;
        std::cout << "  i_0 (rad): " << std::fixed << std::setprecision(15) << eph.i_0 << std::endl;
        std::cout << "  Crc (m): " << std::fixed << std::setprecision(3) << eph.Crc << std::endl;
        std::cout << "  omega (rad): " << std::fixed << std::setprecision(15) << eph.omega << std::endl;
        std::cout << "  Omega_dot (rad/s): " << std::fixed << std::setprecision(15) << eph.Omega_dot << std::endl;
        std::cout << "  I_dot (rad/s): " << std::fixed << std::setprecision(15) << eph.I_dot << std::endl;
        std::cout << "  WN: " << eph.WN << std::endl;
        std::cout << "  TGD (s): " << std::fixed << std::setprecision(15) << eph.TGD << std::endl;

        std::cout << std::endl; // 空行分隔不同卫星
    }

    // 打印统计信息
    std::cout << "========================================" << std::endl;
    std::cout << "星历解析完成：共解析 " << valid_count << " 颗有效卫星星历" << std::endl;
    std::cout << "========================================\n" << std::endl;
}