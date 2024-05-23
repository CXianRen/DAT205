#ifndef __TEST_COMMON_H__

#define __TEST_COMMON_H__

#include <iostream>

#define TEST_LOG(msg) \
    std::cout << "[TEST]:" << msg << std::endl

// measure time
#include <chrono>

#define TIME_START \
    auto time_check_start = std::chrono::high_resolution_clock::now();

#define TIME_END(msg)                                                   \
    auto time_check_end = std::chrono::high_resolution_clock::now();    \
    std::cout << "[TIME]:" << msg << " "                                \
              << std::chrono::duration_cast<std::chrono::milliseconds>( \
                     time_check_end - time_check_start)                 \
                     .count()                                           \
              << "ms" << std::endl;

#endif