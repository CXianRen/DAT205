#include <test_common.h>
#include <common/mmath.h>
#include "test.h"

void test_acc3d()
{
    TEST_LOG("Test ACC3D");

    // access 3d in pos: 10,10,10, size: 20,20,20 data, ACC3D returns 10+10*20+10*20*20
    if (ACC3D(10, 10, 10, 20, 20) != 4210)
    {
        TEST_LOG("ACC3D failed, " << "expected: " << 4010 << ", got: " << ACC3D(10, 10, 10, 20, 20));
        return;
    }
    // access 0,0,0
    if (ACC3D(0, 0, 0, 20, 20) != 0)
    {
        TEST_LOG("ACC3D failed, " << "expected: " << 0 << ", got: " << ACC3D(0, 0, 0, 20, 20));
        return;
    }

    // access 10,0,0
    if (ACC3D(10, 0, 0, 20, 20) != 10)
    {
        TEST_LOG("ACC3D failed, " << "expected: " << 10 << ", got: " << ACC3D(10, 0, 0, 20, 20));
        return;
    }

    TEST_LOG("\t\tACC3D passed");
}

void test_mmath()
{
    test_acc3d();
}