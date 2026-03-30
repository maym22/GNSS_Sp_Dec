#pragma once
#ifndef CONSTS_HPP
#define CONSTS_HPP

namespace consts {
    constexpr double SEQ_LIM = 3;

    constexpr double PI = 3.141592653589793;

    constexpr double GPS_FREQ = 1575.42e6;   // GPS载波频率(Hz)
    constexpr double BDS_FREQ = 1561.098e6;  // BDS载波频率(Hz)

    constexpr double SPEED_OF_LIGHT = 299792458.0;
    constexpr double INTERMEDIATE_FREQ = 2556000.0;
    constexpr double RADIO_FREQ = 1575420000.0;
    constexpr double CODE_FREQ = 1023000.0;
    constexpr double SAMPLE_FREQ = 15360000.0;
    constexpr double EARTH_ROTATION_RATE = 7.2921151467e-5;

    constexpr double SEMI_MAJOR_AXIS = 6378137.0;
    constexpr double SEMI_MINOR_AXIS = 6356752.314245;
    constexpr double FLATTENING = 1.0 / 298.257223563;
    constexpr double FIRST_ECCENTRICITY_SQUARED = 6.69437999014e-3;
    constexpr double SECOND_ECCENTRICITY_SQUARED = 6.73949674228e-3;
    constexpr double EARTH_ANGULAR_VELOCITY = 7.2921151467e-5;
    constexpr double EQUATORIAL_GRAVITY = 9.7803253359;
}

#endif // SAT_PVT_SOLVER_HPP