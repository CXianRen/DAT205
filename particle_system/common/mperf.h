#ifndef __M_PERF_H__
#define __M_PERF_H__

#include <vector>
#include <string>

#include <chrono>

void new_measurement(const std::string &d, double t);

void clear_measurement();

std::string get_mesaurement_info();

#define T_START \
    {           \
        auto start = std::chrono::high_resolution_clock::now();

#define T_END(x)                                                      \
    auto end = std::chrono::high_resolution_clock::now();             \
    std::chrono::duration<double, std::milli> elapsed = end - start;  \
    new_measurement(x, elapsed.count());                              \
    }

#endif // __M_PERF_H__
