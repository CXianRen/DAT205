#include "Solver.h"

#include <stdio.h>
#include <stdlib.h>

#include <perf.h>
#include "constants.h"
#include "common/debug.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#define checkCudaErrors(err) err;

void genLaplace(int *row_ptr, int *col_ind, float *val, int M, int N, int nz,
                float *rhs)
{
    assert(M == N);
    int n = (int)sqrt((double)N);
    assert(n * n == N);
    printf("laplace dimension = %d\n", n);
    int idx = 0;

    // loop over degrees of freedom
    for (int i = 0; i < N; i++)
    {
        int ix = i % n;
        int iy = i / n;

        row_ptr[i] = idx;

        // up
        if (iy > 0)
        {
            val[idx] = 1.0;
            col_ind[idx] = i - n;
            idx++;
        }
        else
        {
            rhs[i] -= 1.0;
        }

        // left
        if (ix > 0)
        {
            val[idx] = 1.0;
            col_ind[idx] = i - 1;
            idx++;
        }
        else
        {
            rhs[i] -= 0.0;
        }

        // center
        val[idx] = -4.0;
        col_ind[idx] = i;
        idx++;

        // right
        if (ix < n - 1)
        {
            val[idx] = 1.0;
            col_ind[idx] = i + 1;
            idx++;
        }
        else
        {
            rhs[i] -= 0.0;
        }

        // down
        if (iy < n - 1)
        {
            val[idx] = 1.0;
            col_ind[idx] = i + n;
            idx++;
        }
        else
        {
            rhs[i] -= 0.0;
        }
    }

    row_ptr[N] = idx;
}

void cudaSolve(
    Eigen::SparseMatrix<float, Eigen::RowMajor> &At,
    Eigen::VectorXf &bt,
    Eigen::VectorXf &xt)
{
    //--------------------------------------------------------------------------

    const int max_iter = 1000;
    int k, M = 0, N = 0, nz = 0, *I = NULL, *J = NULL;
    int *d_col, *d_row;
    int qatest = 0;
    const float tol = 1e-5f;
    float *x, *rhs;
    float r0, r1, alpha, beta;
    float *d_val, *d_x;
    float *d_zm1, *d_zm2, *d_rm2;
    float *d_r, *d_p, *d_omega, *d_y;
    float *val = NULL;
    float *d_valsILU0;
    float rsum, diff, err = 0.0;
    float qaerr1, qaerr2 = 0.0;
    float dot, numerator, denominator, nalpha;
    const float floatone = 1.0;
    const float floatzero = 0.0;

    int nErrors = 0;

    printf("Nz = %d\n", At.nonZeros());
    
    //  I : row pointer, J : column pointer, val : value pointer,
    //  N : number of rows, nz : number of non-zero elements

    // N = At.rows();
    // nz = At.nonZeros();
    // I = At.outerIndexPtr();
    // J = At.innerIndexPtr();
    // val = At.valuePtr();

    // x = xt.data();
    // rhs = bt.data();
    // printf("rhs size = %d\n", bt.size());
    // printf("x size = %d\n", xt.size());

    /* Generate a Laplace matrix in CSR (Compressed Sparse Row) format */

    M = N = 32*64*32;
    nz = 5 * N - 4 * (int)sqrt((double)N);
    I = (int *)malloc(sizeof(int) * (N + 1));   // csr row pointers for matrix A
    J = (int *)malloc(sizeof(int) * nz);       // csr column indices for matrix A
    val = (float *)malloc(sizeof(float) * nz); // csr values for matrix A
    x = (float *)malloc(sizeof(float) * N);
    rhs = (float *)malloc(sizeof(float) * N);    


    for (int i = 0; i < N; i++)
    {
        rhs[i] = 0.0;  // Initialize RHS
        x[i] = 0.0;    // Initial solution approximation
    }

    genLaplace(I, J, val, M, N, nz, rhs);

    printf("nz = %d\n", nz);

    /* Create CUBLAS context */
    cublasHandle_t cublasHandle = NULL;
    checkCudaErrors(cublasCreate(&cublasHandle));

    /* Create CUSPARSE context */
    cusparseHandle_t cusparseHandle = NULL;
    checkCudaErrors(cusparseCreate(&cusparseHandle));

    /* Description of the A matrix */
    cusparseMatDescr_t descr = 0;
    checkCudaErrors(cusparseCreateMatDescr(&descr));
    checkCudaErrors(cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL));
    checkCudaErrors(cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO));

    /* Allocate required memory */
    checkCudaErrors(cudaMalloc((void **)&d_col, nz * sizeof(int)));
    checkCudaErrors(cudaMalloc((void **)&d_row, (N + 1) * sizeof(int)));
    checkCudaErrors(cudaMalloc((void **)&d_val, nz * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_x, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_y, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_r, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_p, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_omega, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_valsILU0, nz * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_zm1, (N) * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_zm2, (N) * sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_rm2, (N) * sizeof(float)));

    /* Wrap raw data into cuSPARSE generic API objects */
    cusparseDnVecDescr_t vecp = NULL, vecX = NULL, vecY = NULL, vecR = NULL, vecZM1 = NULL;
    checkCudaErrors(cusparseCreateDnVec(&vecp, N, d_p, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecX, N, d_x, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecY, N, d_y, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecR, N, d_r, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecZM1, N, d_zm1, CUDA_R_32F));
    cusparseDnVecDescr_t vecomega = NULL;
    checkCudaErrors(cusparseCreateDnVec(&vecomega, N, d_omega, CUDA_R_32F));

    /* Initialize problem data */
    checkCudaErrors(cudaMemcpy(
        d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(
        d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(
        d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(
        d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(
        d_x, x, N * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(
        d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice));

    cusparseSpMatDescr_t matA = NULL;
    cusparseSpMatDescr_t matM_lower, matM_upper;
    cusparseFillMode_t fill_lower = CUSPARSE_FILL_MODE_LOWER;
    cusparseDiagType_t diag_unit = CUSPARSE_DIAG_TYPE_UNIT;
    cusparseFillMode_t fill_upper = CUSPARSE_FILL_MODE_UPPER;
    cusparseDiagType_t diag_non_unit = CUSPARSE_DIAG_TYPE_NON_UNIT;

    checkCudaErrors(cusparseCreateCsr(
        &matA, N, N, nz, d_row, d_col, d_val, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));

    /* Copy A data to ILU(0) vals as input*/
    checkCudaErrors(cudaMemcpy(
        d_valsILU0, d_val, nz * sizeof(float), cudaMemcpyDeviceToDevice));

    // Lower Part
    checkCudaErrors(cusparseCreateCsr(&matM_lower, N, N, nz, d_row, d_col, d_valsILU0,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));

    checkCudaErrors(cusparseSpMatSetAttribute(matM_lower,
                                              CUSPARSE_SPMAT_FILL_MODE,
                                              &fill_lower, sizeof(fill_lower)));
    checkCudaErrors(cusparseSpMatSetAttribute(matM_lower,
                                              CUSPARSE_SPMAT_DIAG_TYPE,
                                              &diag_unit, sizeof(diag_unit)));
    // M_upper
    checkCudaErrors(cusparseCreateCsr(&matM_upper, N, N, nz, d_row, d_col, d_valsILU0,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
    checkCudaErrors(cusparseSpMatSetAttribute(matM_upper,
                                              CUSPARSE_SPMAT_FILL_MODE,
                                              &fill_upper, sizeof(fill_upper)));
    checkCudaErrors(cusparseSpMatSetAttribute(matM_upper,
                                              CUSPARSE_SPMAT_DIAG_TYPE,
                                              &diag_non_unit,
                                              sizeof(diag_non_unit)));

    /* Create ILU(0) info object */
    int bufferSizeLU = 0;
    size_t bufferSizeMV, bufferSizeL, bufferSizeU;
    void *d_bufferLU, *d_bufferMV, *d_bufferL, *d_bufferU;
    cusparseSpSVDescr_t spsvDescrL, spsvDescrU;
    cusparseMatDescr_t matLU;
    csrilu02Info_t infoILU = NULL;

    checkCudaErrors(cusparseCreateCsrilu02Info(&infoILU));
    checkCudaErrors(cusparseCreateMatDescr(&matLU));
    checkCudaErrors(cusparseSetMatType(matLU, CUSPARSE_MATRIX_TYPE_GENERAL));
    checkCudaErrors(cusparseSetMatIndexBase(matLU, CUSPARSE_INDEX_BASE_ZERO));

    /* Allocate workspace for cuSPARSE */
    checkCudaErrors(cusparseSpMV_bufferSize(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matA,
        vecp, &floatzero, vecomega, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT,
        &bufferSizeMV));
    checkCudaErrors(cudaMalloc(&d_bufferMV, bufferSizeMV));

    checkCudaErrors(cusparseScsrilu02_bufferSize(
        cusparseHandle, N, nz, matLU, d_val, d_row, d_col, infoILU, &bufferSizeLU));
    checkCudaErrors(cudaMalloc(&d_bufferLU, bufferSizeLU));

    checkCudaErrors(cusparseSpSV_createDescr(&spsvDescrL));
    checkCudaErrors(cusparseSpSV_bufferSize(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matM_lower, vecR, vecX, CUDA_R_32F,
        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, &bufferSizeL));
    checkCudaErrors(cudaMalloc(&d_bufferL, bufferSizeL));

    checkCudaErrors(cusparseSpSV_createDescr(&spsvDescrU));
    checkCudaErrors(cusparseSpSV_bufferSize(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matM_upper, vecR, vecX, CUDA_R_32F,
        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrU, &bufferSizeU));
    checkCudaErrors(cudaMalloc(&d_bufferU, bufferSizeU));

    /* Preconditioned Conjugate Gradient using ILU.
       --------------------------------------------
       Follows the description by Golub & Van Loan,
       "Matrix Computations 3rd ed.", Algorithm 10.3.1  */

    // printf("\nConvergence of CG using ILU(0) preconditioning: \n");

    /* Perform analysis for ILU(0) */
    checkCudaErrors(cusparseScsrilu02_analysis(
        cusparseHandle, N, nz, descr, d_valsILU0, d_row, d_col, infoILU,
        CUSPARSE_SOLVE_POLICY_USE_LEVEL, d_bufferLU));

    /* generate the ILU(0) factors */
    checkCudaErrors(cusparseScsrilu02(
        cusparseHandle, N, nz, matLU, d_valsILU0, d_row, d_col, infoILU,
        CUSPARSE_SOLVE_POLICY_USE_LEVEL, d_bufferLU));

    /* perform triangular solve analysis */
    checkCudaErrors(cusparseSpSV_analysis(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone,
        matM_lower, vecR, vecX, CUDA_R_32F,
        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, d_bufferL));

    checkCudaErrors(cusparseSpSV_analysis(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone,
        matM_upper, vecR, vecX, CUDA_R_32F,
        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrU, d_bufferU));

    /* reset the initial guess of the solution to zero */
    for (int i = 0; i < N; i++)
    {
        x[i] = 0.0;
    }
    checkCudaErrors(cudaMemcpy(
        d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(
        d_x, x, N * sizeof(float), cudaMemcpyHostToDevice));

    k = 0;
    checkCudaErrors(cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1));

    while (r1 > tol * tol && k <= max_iter)
    {
        // preconditioner application: d_zm1 = U^-1 L^-1 d_r
        checkCudaErrors(cusparseSpSV_solve(cusparseHandle,
                                           CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone,
                                           matM_lower, vecR, vecY, CUDA_R_32F,
                                           CUSPARSE_SPSV_ALG_DEFAULT,
                                           spsvDescrL));

        checkCudaErrors(cusparseSpSV_solve(cusparseHandle,
                                           CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matM_upper,
                                           vecY, vecZM1,
                                           CUDA_R_32F,
                                           CUSPARSE_SPSV_ALG_DEFAULT,
                                           spsvDescrU));
        k++;

        if (k == 1)
        {
            checkCudaErrors(cublasScopy(cublasHandle, N, d_zm1, 1, d_p, 1));
        }
        else
        {
            checkCudaErrors(cublasSdot(
                cublasHandle, N, d_r, 1, d_zm1, 1, &numerator));
            checkCudaErrors(cublasSdot(
                cublasHandle, N, d_rm2, 1, d_zm2, 1, &denominator));
            beta = numerator / denominator;
            checkCudaErrors(cublasSscal(cublasHandle, N, &beta, d_p, 1));
            checkCudaErrors(cublasSaxpy(
                cublasHandle, N, &floatone, d_zm1, 1, d_p, 1));
        }

        checkCudaErrors(cusparseSpMV(
            cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matA,
            vecp, &floatzero, vecomega, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT,
            d_bufferMV));
        checkCudaErrors(cublasSdot(
            cublasHandle, N, d_r, 1, d_zm1, 1, &numerator));
        checkCudaErrors(cublasSdot(
            cublasHandle, N, d_p, 1, d_omega, 1, &denominator));
        alpha = numerator / denominator;
        checkCudaErrors(cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1));
        checkCudaErrors(cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1));
        checkCudaErrors(cublasScopy(cublasHandle, N, d_zm1, 1, d_zm2, 1));
        nalpha = -alpha;
        checkCudaErrors(cublasSaxpy(
            cublasHandle, N, &nalpha, d_omega, 1, d_r, 1));
        checkCudaErrors(cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1));
    }

    printf("  iteration = %3d, residual = %e \n", k, sqrt(r1));

    checkCudaErrors(cudaMemcpy(
        x, d_x, N * sizeof(float), cudaMemcpyDeviceToHost));

    /* Destroy descriptors */
    checkCudaErrors(cusparseDestroyCsrilu02Info(infoILU));
    checkCudaErrors(cusparseDestroyMatDescr(matLU));
    checkCudaErrors(cusparseSpSV_destroyDescr(spsvDescrL));
    checkCudaErrors(cusparseSpSV_destroyDescr(spsvDescrU));
    checkCudaErrors(cusparseDestroySpMat(matM_lower));
    checkCudaErrors(cusparseDestroySpMat(matM_upper));
    checkCudaErrors(cusparseDestroySpMat(matA));
    checkCudaErrors(cusparseDestroyDnVec(vecp));
    checkCudaErrors(cusparseDestroyDnVec(vecomega));
    checkCudaErrors(cusparseDestroyDnVec(vecR));
    checkCudaErrors(cusparseDestroyDnVec(vecX));
    checkCudaErrors(cusparseDestroyDnVec(vecY));
    checkCudaErrors(cusparseDestroyDnVec(vecZM1));

    /* Destroy contexts */
    checkCudaErrors(cusparseDestroy(cusparseHandle));
    checkCudaErrors(cublasDestroy(cublasHandle));

    /* Free device memory */
    free(I);
    free(J);
    free(val);
    free(x);
    free(rhs);
    checkCudaErrors(cudaFree(d_bufferMV));
    checkCudaErrors(cudaFree(d_bufferLU));
    checkCudaErrors(cudaFree(d_bufferL));
    checkCudaErrors(cudaFree(d_bufferU));
    checkCudaErrors(cudaFree(d_col));
    checkCudaErrors(cudaFree(d_row));
    checkCudaErrors(cudaFree(d_val));
    checkCudaErrors(cudaFree(d_x));
    checkCudaErrors(cudaFree(d_y));
    checkCudaErrors(cudaFree(d_r));
    checkCudaErrors(cudaFree(d_p));
    checkCudaErrors(cudaFree(d_omega));
    checkCudaErrors(cudaFree(d_valsILU0));
    checkCudaErrors(cudaFree(d_zm1));
    checkCudaErrors(cudaFree(d_zm2));
    checkCudaErrors(cudaFree(d_rm2));

    // cudaDeviceReset causes the driver to clean up all state. While
    // not mandatory in normal operation, it is good practice.  It is also
    // needed to ensure correct operation when the application is being
    // profiled. Calling cudaDeviceReset causes all profile data to be
    // flushed before the application exits
    cudaDeviceReset();
}