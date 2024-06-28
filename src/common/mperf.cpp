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

void update_measurement(const std::string &d, double t){
    for (int i = 0; i < desc.size(); i++)
    {
        if (desc[i] == d)
        {
            times[i] = t;
            return;
        }
    }
    new_measurement(d, t);
}

void clear_measurement()
{
    times.clear();
    desc.clear();
}

std::string get_mesaurement_info()
{
    std::stringstream ss;

    for (int i = 0; i < desc.size(); i++)
    {
        ss << desc[i] << ": " << times[i] << " ms\n";
    }

    return ss.str();
}
