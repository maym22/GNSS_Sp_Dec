#include "data_io/ObservationReader.hpp"

// 跳过观测文件头（直到END OF HEADER）
int ObsReader::skipObsHeader(const std::vector<std::string>& flines) {
    for (int i = 0; i < flines.size(); ++i) {
        if (flines[i].find("END OF HEADER") != std::string::npos) {
            return i + 1; // 跳过header行
        }
    }
    throw std::runtime_error("观测文件未找到END OF HEADER标记");
}

// 解析时间数组rec_d（对应Matlab sscanf(fline(2:end), "%d%d%d%d%d%f",6)'）
std::vector<int> ObsReader::parseRecD(const std::string& line) {
    std::vector<int> rec_d(6, 0);
    std::istringstream ss(line.substr(1)); // 跳过'>'号
    ss >> rec_d[0] >> rec_d[1] >> rec_d[2] >> rec_d[3] >> rec_d[4] >> rec_d[5];
    if (ss.fail()) {
        throw std::runtime_error("时间行解析失败：" + line);
    }
    return rec_d;
}

// 解析时间数组rec_d（对应Matlab sscanf(fline(2:end), "%d%d%d%d%d%f",6)'）
std::vector<double> ObsReader::parseRecD_t(const std::string& line) {
    std::vector<double> rec_t(6, 0);
    std::istringstream ss(line.substr(1)); // 跳过'>'号
    ss >> rec_t[0] >> rec_t[1] >> rec_t[2] >> rec_t[3] >> rec_t[4] >> rec_t[5];
    if (ss.fail()) {
        throw std::runtime_error("时间行解析失败：" + line);
    }
    return rec_t;
}

// 时间数组转TOW（匹配Matlab dateTool dt2WNTOW_Obs)
double ObsReader::dt2WNTOW_Obs(const std::vector<int>& rec_d, char sys) {
    // rec_d = [年,月,日,时,分,秒]
    std::tm tm_time = {0};
    tm_time.tm_year = rec_d[0] - 1900;
    tm_time.tm_mon = rec_d[1] - 1;
    tm_time.tm_mday = rec_d[2];
    tm_time.tm_hour = rec_d[3];
    tm_time.tm_min = rec_d[4];
    tm_time.tm_sec = rec_d[5];

    time_t t = mktime(&tm_time);
    if (t == -1) return 0.0;

    // 定义GPS/BDS起始epoch
    time_t epoch;
    std::tm epoch_tm = {0};
    if (sys == 'G') { // GPS epoch: 1980-01-06 00:00:00
        epoch_tm.tm_year = 80;
        epoch_tm.tm_mon = 0;
        epoch_tm.tm_mday = 6;
    } else if (sys == 'C') { // BDS epoch: 2006-01-01 00:00:00
        epoch_tm.tm_year = 106;
        epoch_tm.tm_mon = 0;
        epoch_tm.tm_mday = 1;
    } else {
        return 0.0;
    }
    epoch = mktime(&epoch_tm);

    // 闰秒修正（GPST=UTC+18s, BDS=UTC+4s）星历已经转过，这里就不需要了
    const int LEAP_SEC_GPS = 0;
    const int LEAP_SEC_BDS = 0;
    if (sys == 'G') t += LEAP_SEC_GPS;
    if (sys == 'C') t += LEAP_SEC_BDS;

    // 计算TOW（周内秒数）
    double sec_since_epoch = std::difftime(t, epoch);
    return fmod(sec_since_epoch, 7 * 24 * 3600);
}

// 解析单行观测数据（匹配Matlab fline解析逻辑）
bool ObsReader::parseObsLine(const std::string& line, const std::vector<int>& rec_d, const double ObsTimeRes, ObsRawData& obs_data) {
    // 1. 解析系统类型
    obs_data.sys = line[0];

    // 2. 过滤未启用的系统（匹配Matlab use_sys）
    const auto& sat_cfg = cfg_.getSatelliteConfig();
    bool is_sys_enabled = false;
    if ((obs_data.sys == 'G' && sat_cfg.useGPS) || (obs_data.sys == 'C' && sat_cfg.useBDS)) {
        is_sys_enabled = true;
    }
    if (!is_sys_enabled) {
        return false;
    }

    // 3. 解析PRN（BDS偏移32，匹配Matlab逻辑）
    int prn_raw = std::stoi(line.substr(1, 2)); // fline(2:3)
    if (obs_data.sys == 'C' && sat_cfg.useGPS) {
        prn_raw += 32;
    }
    obs_data.PRN = prn_raw;

    // 4. 时间相关（匹配Matlab Time/ObsTime）
    obs_data.Time = rec_d;
    obs_data.ObsTime = dt2WNTOW_Obs(rec_d, obs_data.sys) + ObsTimeRes;

    // 5. 载波频率（匹配Matlab Fc）
    if (obs_data.sys == 'G') {
        obs_data.Fc = consts::GPS_FREQ;
    } else if (obs_data.sys == 'C') {
        obs_data.Fc = consts::BDS_FREQ;
    }

    // 6. 解析观测值（Rho/AcPh/Fd/CNR，对应Matlab fline((-9:3)+o*16)）
    // 格式说明：每16个字符为一个字段，依次是Rho(1-16), AcPh(17-32), Fd(33-48), CNR(49-64)
    auto getField = [&] -> std::vector<std::string> {
        // 步骤1：按空格分割整行到临时数组（核心修改）
        std::vector<std::string> fields;
        std::istringstream iss(line);
        std::string field;
        while (iss >> field) { // 自动跳过任意数量空格，提取非空字段
            fields.push_back(field);
        }
        
        return fields;

    };

    std::vector<std::string> data = getField();
    obs_data.Rho  = std::stod(data[1]);  // Rho: 第2个字段
    obs_data.AcPh = std::stod(data[2]);  // AcPh: 第3个字段
    obs_data.Fd   = std::stod(data[3]);  // Fd: 第4个字段
    obs_data.CNR  = std::stod(data[4]);  // CNR: 第5个字段

    return true;
}

// 核心接口：读取观测量（完全匹配Matlab readObs逻辑）
std::vector<std::vector<ObsRawData>> ObsReader::readObs() {
    // 1. 从统一配置获取参数
    const auto& obs_cfg  = cfg_.getObservationConfig();
    const auto& sat_cfg  = cfg_.getSatelliteConfig();
    std::string obs_path = obs_cfg.data_path;
    int obs_max          = obs_cfg.ObsMax;

    // 2. 读取文件所有行
    std::vector<std::string> flines;
    std::ifstream obs_file(obs_path);
    if (!obs_file.is_open()) {
        throw std::runtime_error("无法打开观测文件：" + obs_path);
    }
    std::string line;
    while (std::getline(obs_file, line)) {
        flines.push_back(line);
    }
    obs_file.close();

    // 3. 跳过文件头
    int line_idx = skipObsHeader(flines);

    // 4. 按块解析观测数据（匹配Matlab逻辑）
    std::vector<std::vector<ObsRawData>> obs_output;
    std::vector<int> rec_d;
    std::vector<double> rec_t;
    int obs_idx = 0;
    double ObsTimeRes = 0;

    while (line_idx < flines.size() && obs_idx < obs_max) {
        line = flines[line_idx];
        line_idx++;

        // 解析时间块（以>开头）
        if (line[0] == '>') {
            obs_idx++;
            rec_d = parseRecD(line);
            rec_t = parseRecD_t(line);
            ObsTimeRes = rec_t[5] - (double)rec_d[5];
            /*
            std::string& line1 = line;
            
            std::vector<int>    rec_d(6, 0);
            std::vector<double> rec_t(6, 0);
            std::istringstream ss(line1.substr(1)); // 跳过'>'号
            ss >> rec_d[0] >> rec_d[1] >> rec_d[2] >> rec_d[3] >> rec_d[4] >> rec_d[5];
            // ss >> rec_t[0] >> rec_t[1] >> rec_t[2] >> rec_t[3] >> rec_t[4] >> rec_t[5];
            if (ss.fail()) {
                    throw std::runtime_error("时间行解析失败：" + line);
            }
            // ObsTimeRes = rec_t[5] - (double)rec_d[5];
            */
            obs_output.emplace_back(); // 新增一个数据块

            continue;
        }

        // 解析卫星观测行
        if (obs_idx == 0 || rec_d.empty()) continue; // 未初始化时间块跳过
        ObsRawData obs_data;
        if (parseObsLine(line, rec_d, ObsTimeRes, obs_data)) {
            obs_output[obs_idx - 1].push_back(obs_data);
        }
    }

    // 5. 截断到实际读取数量
    obs_output.resize(obs_idx);

    // 6. 日志输出（匹配Matlab logger）
    std::cout << "[观测量读取] 成功读取 " << obs_idx << " 个数据块" << std::endl;
    if (!obs_output.empty()) {
        // 输出时间范围
        auto& first_time = obs_output[0][0].Time;
        auto& last_time = obs_output.back()[0].Time;
        std::cout << "[观测量读取] 时间范围：" 
                  << first_time[0] << "-" << first_time[1] << "-" << first_time[2] << " "
                  << first_time[3] << ":" << first_time[4] << ":" << first_time[5] 
                  << " 至 "
                  << last_time[0] << "-" << last_time[1] << "-" << last_time[2] << " "
                  << last_time[3] << ":" << last_time[4] << ":" << last_time[5] << std::endl;

        // 统计最大/最小观测数
        int max_obs = 0, min_obs = 0;
        for (const auto& block : obs_output) {
            max_obs = std::max(max_obs, (int)block.size());
            min_obs = std::min(min_obs, (int)block.size());
        }
        std::cout << "[观测量读取] 最大/最小观测数：" << max_obs << "/" << min_obs << std::endl;
    }

    return obs_output;
}

// 打印函数：完全匹配Matlab输出格式，便于核对
void ObsReader::printObsData(const std::vector<std::vector<ObsRawData>>& obs_output) const {
    if (obs_output.empty()) {
        std::cout << "\n[观测量打印] 无有效观测数据！" << std::endl;
        return;
    }

    std::cout << "\n======================================= 观测数据 =======================================" << std::endl;
    for (int block_idx = 0; block_idx < obs_output.size(); ++block_idx) {
        const auto& block = obs_output[block_idx];
        std::cout << "\n---------- 数据块 " << block_idx + 1 << " ----------" << std::endl;
        
        // 打印块时间
        const auto& rec_d = block[0].Time;
        std::cout << "时间：" << rec_d[0] << "-" << rec_d[1] << "-" << rec_d[2] << " "
                  << rec_d[3] << ":" << rec_d[4] << ":" << rec_d[5] 
                  << " (TOW: " << std::fixed << std::setprecision(3) << block[0].ObsTime << "s)" << std::endl;

        // 打印表头
        std::cout << std::setw(4) << "Sys" 
                  << std::setw(6) << "PRN" 
                  << std::setw(18) << "Fc(Hz)" 
                  << std::setw(15) << "Rho(m)" 
                  << std::setw(18) << "AcPh(cycle)" 
                  << std::setw(12) << "Fd(Hz)" 
                  << std::setw(10) << "CNR(dB-Hz)" << std::endl;
        std::cout << "----------------------------------------------------------------------------------------" << std::endl;

        // 打印块内卫星数据
        for (const auto& obs : block) {
            std::cout << std::setw(4) << obs.sys 
                      << std::setw(6) << obs.PRN 
                      << std::setw(18) << std::fixed << std::setprecision(3) << obs.Fc/1e6 << "e6" 
                      << std::setw(15) << std::fixed << std::setprecision(3) << obs.Rho 
                      << std::setw(18) << std::fixed << std::setprecision(6) << obs.AcPh 
                      << std::setw(12) << std::fixed << std::setprecision(3) << obs.Fd 
                      << std::setw(10) << std::fixed << std::setprecision(1) << obs.CNR << std::endl;
        }
    }

    // 统计信息
    std::cout << "\n========================================================================================" << std::endl;
    std::cout << "总计数据块数：" << obs_output.size() << std::endl;
    int total_sat = 0;
    for (const auto& block : obs_output) total_sat += block.size();
    std::cout << "总计卫星观测数：" << total_sat << std::endl;
    std::cout << "========================================================================================" << std::endl;
}