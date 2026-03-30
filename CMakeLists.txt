cmake_minimum_required(VERSION 3.10)
project(SignalDetection)

# C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall")


# 查找yaml-cpp库
find_library(YAML_CPP_LIBRARIES NAMES yaml-cpp)
if(NOT YAML_CPP_LIBRARIES)
    message(FATAL_ERROR "yaml-cpp library not found!")
endif()

# 查找GSL库
find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
    message(FATAL_ERROR "GSL library not found!")
endif()

# 头文件路径
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/usr/include)  # 确保yaml-cpp头文件可找到
include_directories(${GSL_INCLUDE_DIRS})

# 源文件列表
set(SOURCES
    src/main.cpp
    src/data_io/EphemerisReader.cpp
    src/data_io/ObservationReader.cpp
    src/antispoofing/AntiSpoofing.cpp
    src/antispoofing/AsFunSchannel.cpp
    src/antispoofing/AsFunMchannel.cpp
    src/antispoofing/AsFunPVT.cpp
    src/pvt/SatPvtSolver.cpp
    src/pvt/UserPvtSolver.cpp
    src/config/ConfigReader.cpp
)

# 生成可执行文件
add_executable(SignalDetection ${SOURCES})

# 链接yaml-cpp库
target_link_libraries(SignalDetection ${YAML_CPP_LIBRARIES})
target_link_libraries(SignalDetection ${GSL_LIBRARIES})

# 复制配置文件到build目录（方便运行）
# configure_file(${PROJECT_SOURCE_DIR}/config/config.yaml ${PROJECT_BINARY_DIR}/config/config.yaml COPYONLY)