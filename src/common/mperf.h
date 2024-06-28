#ifndef __M_PERF_H__
#define __M_PERF_H__

#include <vector>
#include <string>

#include <chrono>

void new_measurement(const std::string &d, double t);

void update_measurement(const std::string &d, double t);

void clear_measurement();

std::string get_mesaurement_info();

#define T_START(perf_str)                                       \
    {                                                           \
        auto start = std::chrono::high_resolution_clock::now(); \
        auto t_perf_str = perf_str;                            \
        new_measurement(t_perf_str, 0.0);

#define T_END                                                        \
    auto end = std::chrono::high_resolution_clock::now();            \
    std::chrono::duration<double, std::milli> elapsed = end - start; \
    update_measurement(t_perf_str, elapsed.count());                 \
    }

#endif // __M_PERF_H__
