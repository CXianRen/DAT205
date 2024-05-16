
#include <iostream>

#define gN 4 
#include "particle/const.h"

#define TEST_LOG(msg) \
    std::cout << "[TEST]:"<< msg << std::endl;

#include "particle/simulator.h"


int main()
{
    TEST_LOG("Test!");
    TEST_LOG("gX: " << gX);
    TEST_LOG("gY: " << gY);
    TEST_LOG("gZ: " << gZ);
    TEST_LOG("gSIZE: " << gSIZE);

    auto sim = new Simulator();
    sim->init();
}