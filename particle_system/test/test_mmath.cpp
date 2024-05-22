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

void test_acc2d(){
    TEST_LOG("Test ACC2D");

    // access 2d in pos: 10,10, size: 20,20 data, ACC2D returns 10+10*20
    if (ACC2D(10, 10, 20) != 210)
    {
        TEST_LOG("ACC2D failed, " << "expected: " << 210 << ", got: " << ACC2D(10, 10, 20));
        return;
    }
    // access 0,0
    if (ACC2D(0, 0, 20) != 0)
    {
        TEST_LOG("ACC2D failed, " << "expected: " << 0 << ", got: " << ACC2D(0, 0, 20));
        return;
    }

    // access 10,0
    if (ACC2D(10, 0, 20) != 10)
    {
        TEST_LOG("ACC2D failed, " << "expected: " << 10 << ", got: " << ACC2D(10, 0, 20));
        return;
    }

    TEST_LOG("\t\tACC2D passed");

}

void test_get_gradiant_3d_x(){
    TEST_LOG("Test GET_GRADIANT_3D_X");

    float* data = new float[27];
    data[ACC3D(0,0,0,3,3)] = 1.0f;
    data[ACC3D(1,0,0,3,3)] = 2.0f;
    data[ACC3D(2,0,0,3,3)] = 3.0f;

    auto g  = GET_GRADIANT_3D_X(1,0,0, data, 3, 3, 1.0f);
    if (g != 1.0f)
    {
        TEST_LOG("GET_GRADIANT_3D_X failed, " << "expected: " << 2.0f << ", got: " << g
            << "\n\t\tdata[ACC3D(1+1,0,0,3,3)]: " << data[ACC3D(1+1,0,0,3,3)]
            << "\n\t\tdata[ACC3D(1-1,0,0,3,3)]: " << data[ACC3D(1-1,0,0,3,3)]
        );

        return;
    }
    TEST_LOG("\t\tGET_GRADIANT_3D_X passed");
}

void test_get_gradiant_3d_y(){
    TEST_LOG("Test GET_GRADIANT_3D_Y");

    float* data = new float[27];
    data[ACC3D(0,0,0,3,3)] = 1.0f;
    data[ACC3D(0,1,0,3,3)] = 2.0f;
    data[ACC3D(0,2,0,3,3)] = 3.0f;

    auto g  = GET_GRADIANT_3D_Y(0,1,0, data, 3, 3, 1.0f);
    if (g != 1.0f)
    {
        TEST_LOG("GET_GRADIANT_3D_Y failed, " << "expected: " << 2.0f << ", got: " << g
            << "\n\t\tdata[ACC3D(0,1+1,0,3,3)]: " << data[ACC3D(0,1+1,0,3,3)]
            << "\n\t\tdata[ACC3D(0,1-1,0,3,3)]: " << data[ACC3D(0,1-1,0,3,3)]
        );

        return;
    }
    TEST_LOG("\t\tGET_GRADIANT_3D_Y passed");
}

void test_get_gradiant_3d_z(){
    TEST_LOG("Test GET_GRADIANT_3D_Z");

    float* data = new float[27];
    data[ACC3D(0,0,0,3,3)] = 1.0f;
    data[ACC3D(0,0,1,3,3)] = 2.0f;
    data[ACC3D(0,0,2,3,3)] = 3.0f;

    auto g  = GET_GRADIANT_3D_Z(0,0,1, data, 3, 3, 1.0f);
    if (g != 1.0f)
    {
        TEST_LOG("GET_GRADIANT_3D_Z failed, " << "expected: " << 2.0f << ", got: " << g
            << "\n\t\tdata[ACC3D(0,0,1+1,3,3)]: " << data[ACC3D(0,0,1+1,3,3)]
            << "\n\t\tdata[ACC3D(0,0,1-1,3,3)]: " << data[ACC3D(0,0,1-1,3,3)]
        );

        return;
    }
    TEST_LOG("\t\tGET_GRADIANT_3D_Z passed");
}


void test_get_gradiant_2d_x(){
    TEST_LOG("Test GET_GRADIANT_2D_X");

    float* data = new float[9];
    data[ACC2D(0,0,3)] = 1.0f;
    data[ACC2D(1,0,3)] = 2.0f;
    data[ACC2D(2,0,3)] = 3.0f;

    auto g  = GET_GRADIANT_2D_X(1,0, data, 3, 1.0f);
    if (g != 1.0f)
    {
        TEST_LOG("GET_GRADIANT_2D_X failed, " << "expected: " << 2.0f << ", got: " << g
            << "\n\t\tdata[ACC2D(1+1,0,3)]: " << data[ACC2D(1+1,0,3)]
            << "\n\t\tdata[ACC2D(1-1,0,3)]: " << data[ACC2D(1-1,0,3)]
        );

        return;
    }
    TEST_LOG("\t\tGET_GRADIANT_2D_X passed");

}

void test_get_gradiant_2d_y(){
    TEST_LOG("Test GET_GRADIANT_2D_Y");

    float* data = new float[9];
    data[ACC2D(0,0,3)] = 1.0f;
    data[ACC2D(0,1,3)] = 2.0f;
    data[ACC2D(0,2,3)] = 3.0f;

    auto g  = GET_GRADIANT_2D_Y(0,1, data, 3, 1.0f);
    if (g != 1.0f)
    {
        TEST_LOG("GET_GRADIANT_2D_Y failed, " << "expected: " << 2.0f << ", got: " << g
            << "\n\t\tdata[ACC2D(0,1+1,3)]: " << data[ACC2D(0,1+1,3)]
            << "\n\t\tdata[ACC2D(0,1-1,3)]: " << data[ACC2D(0,1-1,3)]
        );

        return;
    }
    TEST_LOG("\t\tGET_GRADIANT_2D_Y passed");
}

void test_mmath()
{
    test_acc3d();
    test_acc2d();
    test_get_gradiant_3d_x();
    test_get_gradiant_3d_y();
    test_get_gradiant_3d_z();
    test_get_gradiant_2d_y();
}