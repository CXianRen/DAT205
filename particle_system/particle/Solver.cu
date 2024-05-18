#include "Solver.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <cusparse.h>


void cudaSolve()
{
    // problem definition
    const int A_num_rows = 40;
    const int A_num_cols = 40;
    const int A_nnz = 80; // assuming a dense matrix
    int hA_csrOffsets[A_num_rows + 1];
    int hA_columns[A_nnz];
    float hA_values[A_nnz];
    float hX[A_num_cols];
    float hY[A_num_rows];
    float hY_result[A_num_rows];
    float alpha = 1.0f;
    float beta = 0.0f;

    // Initialize your data here. For example:
    for (int i = 0; i < A_num_rows + 1; i++)
    {
        hA_csrOffsets[i] = i * A_num_cols; // assuming a dense matrix
    }
    for (int i = 0; i < A_nnz; i++)
    {
        hA_columns[i] = i % A_num_cols;                // assuming a dense matrix
        hA_values[i] = static_cast<float>(i % 10 + 1); // some arbitrary values
    }
    for (int i = 0; i < A_num_cols; i++)
    {
        hX[i] = static_cast<float>(i % 10 + 1); // some arbitrary values
    }
    for (int i = 0; i < A_num_rows; i++)
    {
        hY[i] = 0.0f;
        hY_result[i] = 0.0f; // you'll need to compute the correct results
    }

    //--------------------------------------------------------------------------
    // Device memory management
    int *dA_csrOffsets, *dA_columns;
    float *dA_values, *dX, *dY;
    cudaMalloc((void **)&dA_csrOffsets, (A_num_rows + 1) * sizeof(int));
    cudaMalloc((void **)&dA_columns, A_nnz * sizeof(int));
    cudaMalloc((void **)&dA_values, A_nnz * sizeof(float));
    cudaMalloc((void **)&dX, A_num_cols * sizeof(float));
    cudaMalloc((void **)&dY, A_num_rows * sizeof(float));

    cudaMemcpy(dA_csrOffsets, hA_csrOffsets, (A_num_rows + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dA_columns, hA_columns, A_nnz * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dA_values, hA_values, A_nnz * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dX, hX, A_num_cols * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dY, hY, A_num_rows * sizeof(float), cudaMemcpyHostToDevice);
    //--------------------------------------------------------------------------
    // CUSPARSE APIs
    cusparseHandle_t handle = NULL;
    cusparseSpMatDescr_t matA;
    cusparseDnVecDescr_t vecX, vecY;
    void *dBuffer = NULL;
    size_t bufferSize = 0;
    cusparseCreate(&handle);
    cusparseCreateCsr(&matA, A_num_rows, A_num_cols, A_nnz,
                      dA_csrOffsets, dA_columns, dA_values,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                      CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
    cusparseCreateDnVec(&vecX, A_num_cols, dX, CUDA_R_32F);
    cusparseCreateDnVec(&vecY, A_num_rows, dY, CUDA_R_32F);
    cusparseSpMV_bufferSize(
        handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
        CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
    cudaMalloc(&dBuffer, bufferSize);

    // execute SpMV
    cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                 &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
                 CUSPARSE_SPMV_ALG_DEFAULT, dBuffer);

    // destroy matrix/vector descriptors
    cusparseDestroySpMat(matA);
    cusparseDestroyDnVec(vecX);
    cusparseDestroyDnVec(vecY);
    cusparseDestroy(handle);
    //--------------------------------------------------------------------------
    // device result check
    cudaMemcpy(hY, dY, A_num_rows * sizeof(float), cudaMemcpyDeviceToHost);
    int correct = 1;
    for (int i = 0; i < A_num_rows; i++)
    {
        if (hY[i] != hY_result[i])
        {                // direct floating point comparison is not
            correct = 0; // reliable
            break;
        }
    }
    if (correct)
        printf("spmv_csr_example test PASSED\n");
    else
        printf("spmv_csr_example test FAILED: wrong result\n");
    //--------------------------------------------------------------------------
    // device memory deallocation
    cudaFree(dBuffer);
    cudaFree(dA_csrOffsets);
    cudaFree(dA_columns);
    cudaFree(dA_values);
    cudaFree(dX);
    cudaFree(dY);
}