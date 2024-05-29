#include "mperf.h"

#include <chrono>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

// performance measurement
std::vector<std::string> desc;
std::vector<double> times;

void new_measurement(const std::string &d, double t)
{
    desc.push_back(d);
    times.push_back(t);
}

void clear_measurement()
{
    times.clear();
    desc.clear();
}

std::string get_mesaurement_info()
{
    std::stringstream ss;
    // total time
    double total_time = 0.0;
    for (int i = 0; i < times.size(); i++)
    {
        total_time += times[i];
    }
    // .3f means 3 decimal places
    ss << "Total time: " << std::setprecision(3) << total_time << " ms\n";

    for (int i = 0; i < desc.size(); i++)
    {
        ss << desc[i] << ": " << times[i] << " ms\n";
    }

    return ss.str();
}
