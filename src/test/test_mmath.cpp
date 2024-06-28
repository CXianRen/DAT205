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

// void test_get_gradiant_3d_x(){
//     TEST_LOG("Test GET_GRADIANT_3D_X");

//     float* data = new float[27];
//     data[ACC3D(0,0,0,3,3)] = 1.0f;
//     data[ACC3D(1,0,0,3,3)] = 2.0f;
//     data[ACC3D(2,0,0,3,3)] = 3.0f;

//     auto g  = GET_GRADIANT_3D_X(1,0,0, data, 3, 3, 1.0f);
//     if (g != 1.0f)
//     {
//         TEST_LOG("GET_GRADIANT_3D_X failed, " << "expected: " << 2.0f << ", got: " << g
//             << "\n\t\tdata[ACC3D(1+1,0,0,3,3)]: " << data[ACC3D(1+1,0,0,3,3)]
//             << "\n\t\tdata[ACC3D(1-1,0,0,3,3)]: " << data[ACC3D(1-1,0,0,3,3)]
//         );

//         return;
//     }
//     TEST_LOG("\t\tGET_GRADIANT_3D_X passed");
// }

// void test_get_gradiant_3d_y(){
//     TEST_LOG("Test GET_GRADIANT_3D_Y");

//     float* data = new float[27];
//     data[ACC3D(0,0,0,3,3)] = 1.0f;
//     data[ACC3D(0,1,0,3,3)] = 2.0f;
//     data[ACC3D(0,2,0,3,3)] = 3.0f;

//     auto g  = GET_GRADIANT_3D_Y(0,1,0, data, 3, 3, 1.0f);
//     if (g != 1.0f)
//     {
//         TEST_LOG("GET_GRADIANT_3D_Y failed, " << "expected: " << 2.0f << ", got: " << g
//             << "\n\t\tdata[ACC3D(0,1+1,0,3,3)]: " << data[ACC3D(0,1+1,0,3,3)]
//             << "\n\t\tdata[ACC3D(0,1-1,0,3,3)]: " << data[ACC3D(0,1-1,0,3,3)]
//         );

//         return;
//     }
//     TEST_LOG("\t\tGET_GRADIANT_3D_Y passed");
// }

// void test_get_gradiant_3d_z(){
//     TEST_LOG("Test GET_GRADIANT_3D_Z");

//     float* data = new float[27];
//     data[ACC3D(0,0,0,3,3)] = 1.0f;
//     data[ACC3D(0,0,1,3,3)] = 2.0f;
//     data[ACC3D(0,0,2,3,3)] = 3.0f;

//     auto g  = GET_GRADIANT_3D_Z(0,0,1, data, 3, 3, 1.0f);
//     if (g != 1.0f)
//     {
//         TEST_LOG("GET_GRADIANT_3D_Z failed, " << "expected: " << 2.0f << ", got: " << g
//             << "\n\t\tdata[ACC3D(0,0,1+1,3,3)]: " << data[ACC3D(0,0,1+1,3,3)]
//             << "\n\t\tdata[ACC3D(0,0,1-1,3,3)]: " << data[ACC3D(0,0,1-1,3,3)]
//         );

//         return;
//     }
//     TEST_LOG("\t\tGET_GRADIANT_3D_Z passed");
// }

// void test_get_gradiant_2d_x(){
//     TEST_LOG("Test GET_GRADIANT_2D_X");

//     float* data = new float[9];
//     data[ACC2D(0,0,3)] = 1.0f;
//     data[ACC2D(1,0,3)] = 2.0f;
//     data[ACC2D(2,0,3)] = 3.0f;

//     auto g  = GET_GRADIANT_2D_X(1,0, data, 3, 1.0f);
//     if (g != 1.0f)
//     {
//         TEST_LOG("GET_GRADIANT_2D_X failed, " << "expected: " << 2.0f << ", got: " << g
//             << "\n\t\tdata[ACC2D(1+1,0,3)]: " << data[ACC2D(1+1,0,3)]
//             << "\n\t\tdata[ACC2D(1-1,0,3)]: " << data[ACC2D(1-1,0,3)]
//         );

//         return;
//     }
//     TEST_LOG("\t\tGET_GRADIANT_2D_X passed");

// }

// void test_get_gradiant_2d_y(){
//     TEST_LOG("Test GET_GRADIANT_2D_Y");

//     float* data = new float[9];
//     data[ACC2D(0,0,3)] = 1.0f;
//     data[ACC2D(0,1,3)] = 2.0f;
//     data[ACC2D(0,2,3)] = 3.0f;

//     auto g  = GET_GRADIANT_2D_Y(0,1, data, 3, 1.0f);
//     if (g != 1.0f)
//     {
//         TEST_LOG("GET_GRADIANT_2D_Y failed, " << "expected: " << 2.0f << ", got: " << g
//             << "\n\t\tdata[ACC2D(0,1+1,3)]: " << data[ACC2D(0,1+1,3)]
//             << "\n\t\tdata[ACC2D(0,1-1,3)]: " << data[ACC2D(0,1-1,3)]
//         );

//         return;
//     }
//     TEST_LOG("\t\tGET_GRADIANT_2D_Y passed");
// }

// void test_2d_lapalace(){
//     TEST_LOG("Test build_2d_laplace");
//     int diag[9] = {2,3,2,3,4,3,2,3,2};
//     auto L = build_2d_laplace<float>(3,3);
//     for (int i = 0; i < 9; i++)
//     {
//         if (L.coeff(i,i) != diag[i])
//         {
//             TEST_LOG("build_2d_laplace failed, " << "expected: " << diag[i] << ", got: " << L.coeff(i,i));
//             return;
//         }
//     }
//     TEST_LOG("\n" << L.toDense());
//     TEST_LOG("\t\tbuild_2d_laplace passed");
// }

// void test_2d_divergence(){
//     TEST_LOG("Test build_2d_divergence");
//     float u[9] = {1,1,1,
//                   1,2,1,
//                   1,1,1
//                 };
//     float v[9] = {1,1,1,
//                   1,2,1,
//                   1,1,1
//                 };

//     auto d = GET_DIVERGENCE_2D(1,1,u,v,3,1.0,1.0);
//     // d = u(i,j+1)- u (i,j) + v(i+1,j) - v(i,j)
//     if (d != -2.f)
//     {
//         TEST_LOG("GET_DIVERGENCE_2D failed, " << "expected: " << -2.f << ", got: " << d
//             << "\n u(2,1): " << u[ACC2D(2,1,3)]
//             << "\n u(1,1): " << u[ACC2D(1,1,3)]
//             << "\n v(1,1): " << v[ACC2D(1,1,3)]
//             << "\n v(1,2): " << v[ACC2D(1,2,3)]);
//         return;
//     }
//     TEST_LOG("\t\tGET_DIVERGENCE_2D passed");

// }

void test_mmath()
{
    test_acc3d();
    // test_acc2d();
    // test_get_gradiant_3d_x();
    // test_get_gradiant_3d_y();
    // test_get_gradiant_3d_z();
    // test_get_gradiant_2d_y();
    // test_2d_lapalace();
    // test_2d_divergence();
}