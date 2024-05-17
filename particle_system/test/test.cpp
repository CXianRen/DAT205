
#include <iostream>
#include <vector>
#define TEST_LOG(msg) \
    std::cout << "[TEST]:"<< msg << std::endl;

#include "common/mmath.h"
// Function to compare two arrays
template <typename T>
bool compareArrays(const T* arr1, const T* arr2, int size) {
    for (int i = 0; i < size; ++i) {
        if (arr1[i] != arr2[i]) {
            TEST_LOG("Arrays are not equal at idx:" << i << " arr1:" << arr1[i] << " arr2:" << arr2[i]);
            return false;
        }
    }
    return true;
}

int main() {

    return 0;
}