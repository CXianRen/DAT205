# how does Cuda arrange threads

1. SM 
   1. the actually worker (processor to handle a group of threads (Block))
2. Block & Warp
   1. threads in a Block will be grouped into Warps, for example, a 128 block, will contains 4 warp, with 32 threads to each one.
   2. max blocks: 

# a cuda kernel example
```cpp

__global__ void addVec(int* A, int*B, int* C, int N){
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(idx < N){
        // do you job
    }
}

int main(){
    // ...
    // task size 
    int N;

    // number of thread per block
    int threadsPerBlock = 256;
    
    // number of block per grid, 
    // make sure N%threadsPerBlock still allocated into a block.
    int blocksPerGrid = (N + threadsPerBlock -1 )/threadsPerBlock;

    // start the kernel
    addVec<<<blocksPerGrid, threadsPerBlock>>>();

}
```
