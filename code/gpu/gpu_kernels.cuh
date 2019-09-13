#ifndef GPU_KERNELS_CUH
#define GPU_KERNELS_CUH

#include "split_pointer.hpp"

#include <cstdio>

template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void rwTest(IndexType n, SplitPointer<NumberSegments, OriginalPointerType> rw, OriginalPointerType factor)
{
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < n) {
        OriginalPointerType val = rw.read<SegmentsToUse>(idx);
        rw.write<SegmentsToUse>(idx, val * factor);
    }
}


template<MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void setVector(IndexType n, SplitPointer<NumberSegments, OriginalPointerType> vec, OriginalPointerType initData)
{
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < n) {
        vec.writeAll(idx, initData);
    }
}

template<typename ValueType, typename IndexType>
__global__ void setVector(IndexType n, ValueType *__restrict__ vec, ValueType initData)
{
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < n) {
        vec[idx] = initData;
    }
}



template<unsigned int blockSize, typename IndexType, typename WorkingType>
__device__ __forceinline__ void blockAddReduce(WorkingType *shared)
{
    IndexType tIdx = threadIdx.x;
    __syncthreads();
    
    #pragma unroll
    for (IndexType i = blockSize >> 1; i > 0; i = i >> 1) {
        if (tIdx < i) {
            shared[tIdx] += shared[tIdx + i];
        }
        __syncthreads();
    }
}


/*
Requirements:
- BlockSize MUST be a power of 2 (and equal to blockDim.x)

Reduction from <size> elements to <gridDim.x>
*/
template<unsigned int blockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void addReduction(IndexType size, const SplitPointer<NumberSegments, OriginalPointerType> in, OriginalPointerType *__restrict__ out)
{
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType tIdx = threadIdx.x;
    
    __shared__ OriginalPointerType shared[blockSize];
    OriginalPointerType sum = static_cast<OriginalPointerType>(0.0);
    for (IndexType i = idx; i < size; i += blockSize * gridDim.x) {
        const OriginalPointerType val = in.read<SegmentsToUse>(i);
        sum += val;
    }
    shared[tIdx] = sum;
    __syncthreads();
    
    blockAddReduce<blockSize, IndexType>(shared);
    
    if (tIdx == 0) {
        out[blockIdx.x] = shared[0];
    }
}


/*
Requirements:
- BlockSize MUST be a power of 2 (and equal to blockDim.x)

Reduction from <size> elements to <gridDim.x>
*/
template<unsigned int blockSize, typename OriginalPointerType, typename IndexType>
__global__ void addReduction(IndexType size, const OriginalPointerType *__restrict__ in, OriginalPointerType *__restrict__ out)
{
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType tIdx = threadIdx.x;
    
    __shared__ OriginalPointerType shared[blockSize];
    OriginalPointerType sum = static_cast<OriginalPointerType>(0.0);
    for (IndexType i = idx; i < size; i += blockSize * gridDim.x) {
        const OriginalPointerType val = in[i];
        sum += val;
    }
    shared[tIdx] = sum;
    __syncthreads();
    
    blockAddReduce<blockSize, IndexType>(shared);
    
    if (tIdx == 0) {
        out[blockIdx.x] = shared[0];
    }
}


template<unsigned int blockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void norm1Step1(IndexType n, const SplitPointer<NumberSegments, OriginalPointerType> in, SplitPointer<NumberSegments, OriginalPointerType> out)
{
    using std::abs;
    __shared__ OriginalPointerType shared[blockSize];
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType tIdx = threadIdx.x;
    
    OriginalPointerType sum = static_cast<OriginalPointerType>(0.0);
    for (IndexType i = idx; i < n; i += blockSize * gridDim.x) {
        sum += abs(in.read<SegmentsToUse>(i));
    }
    
    shared[tIdx] = sum;
    __syncthreads();
    blockAddReduce<blockSize, IndexType>(shared);    //res[blockIdx.x] will be written
    
    if (tIdx == 0) {
        out.write<SegmentsToUse>(blockIdx.x, shared[0]);
    }
}

template<unsigned int blockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void norm1Step1(IndexType n, const SplitPointer<NumberSegments, OriginalPointerType> in, OriginalPointerType *__restrict__ out)
{
    using std::abs;
    __shared__ OriginalPointerType shared[blockSize];
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType tIdx = threadIdx.x;
    
    OriginalPointerType sum = static_cast<OriginalPointerType>(0.0);
    for (IndexType i = idx; i < n; i += blockSize * gridDim.x) {
        sum += abs(in.read<SegmentsToUse>(i));
    }
    
    shared[tIdx] = sum;
    __syncthreads();
    blockAddReduce<blockSize, IndexType>(shared);    //res[blockIdx.x] will be written
    
    if (tIdx == 0) {
        out[blockIdx.x] = shared[0];
    }
}


template<unsigned int blockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void norm1DiffStep1(IndexType n, const SplitPointer<NumberSegments, OriginalPointerType> in1,
    const SplitPointer<NumberSegments, OriginalPointerType> in2,
    SplitPointer<NumberSegments, OriginalPointerType> out)
{
    using std::abs;
    __shared__ OriginalPointerType shared[blockSize];
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType tIdx = threadIdx.x;
    
    OriginalPointerType sum = static_cast<OriginalPointerType>(0.0);
    for (IndexType i = idx; i < n; i += blockSize * gridDim.x) {
        sum += abs(in1.read<SegmentsToUse>(i) - in2.read<SegmentsToUse>(i));
    }
    
    shared[tIdx] = sum;
    __syncthreads();
    blockAddReduce<blockSize, IndexType>(shared);    //res[blockIdx.x] will be written
    
    if (tIdx == 0) {
        out.write<SegmentsToUse>(blockIdx.x, shared[0]);
    }
}

template<unsigned int blockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void norm1DiffStep1(IndexType n, const SplitPointer<NumberSegments, OriginalPointerType> in1,
    const SplitPointer<NumberSegments, OriginalPointerType> in2,
    OriginalPointerType *__restrict__ out)
{
    using std::abs;
    __shared__ OriginalPointerType shared[blockSize];
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType tIdx = threadIdx.x;
    
    OriginalPointerType sum = static_cast<OriginalPointerType>(0.0);
    for (IndexType i = idx; i < n; i += blockSize * gridDim.x) {
        sum += abs(in1.read<SegmentsToUse>(i) - in2.read<SegmentsToUse>(i));
    }
    
    shared[tIdx] = sum;
    __syncthreads();
    blockAddReduce<blockSize, IndexType>(shared);    //res[blockIdx.x] will be written
    
    if (tIdx == 0) {
        out[blockIdx.x] = shared[0];
    }
}


template<unsigned int blockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void selectiveNorm1Step1(
    IndexType indexVectorSize,
    const IndexType * __restrict__ vecIndices,
    const SplitPointer<NumberSegments, OriginalPointerType> vecSplit,
   SplitPointer<NumberSegments, OriginalPointerType> resSplit)
{
    using std::abs;
    __shared__ OriginalPointerType shared[blockSize];
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType tIdx = threadIdx.x;
    
    OriginalPointerType sum = static_cast<OriginalPointerType>(0.0);
    for (IndexType i = idx; i < indexVectorSize; i += blockSize * gridDim.x) {
        const IndexType vecIdx = vecIndices[i];
        sum += abs(vecSplit.read<SegmentsToUse>(vecIdx));
    }
    
    shared[tIdx] = sum;
    __syncthreads();
    blockAddReduce<blockSize, IndexType>(shared);    //res[blockIdx.x] will be written
    
    if (tIdx == 0) {
        resSplit.write<SegmentsToUse>(blockIdx.x, shared[0]);
    }
}

template<unsigned int blockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void selectiveNorm1Step1(
    IndexType indexVectorSize,
    const IndexType * __restrict__ vecIndices,
    const SplitPointer<NumberSegments, OriginalPointerType> vecSplit,
   OriginalPointerType *__restrict__ result)
{
    using std::abs;
    __shared__ OriginalPointerType shared[blockSize];
    const IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    const IndexType tIdx = threadIdx.x;
    
    OriginalPointerType sum = static_cast<OriginalPointerType>(0.0);
    for (IndexType i = idx; i < indexVectorSize; i += blockSize * gridDim.x) {
        const IndexType vecIdx = vecIndices[i];
        sum += abs(vecSplit.read<SegmentsToUse>(vecIdx));
    }
    
    shared[tIdx] = sum;
    __syncthreads();
    blockAddReduce<blockSize, IndexType>(shared);    //res[blockIdx.x] will be written
    
    if (tIdx == 0) {
        result[blockIdx.x] = shared[0];
    }
}


template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void ellPageRank_wki(IndexType n, IndexType k, const IndexType * __restrict__ indices,
    const SplitPointer<NumberSegments, OriginalPointerType> valsSplit,
    OriginalPointerType dampingFactor,
        const SplitPointer<NumberSegments, OriginalPointerType> oldRankSplit,
       SplitPointer<NumberSegments, OriginalPointerType> newRankSplit)
{
    const IndexType row = blockDim.x * blockIdx.x + threadIdx.x;
    if (row < n) {
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        for (IndexType i = 0; i < k; ++i) {
            const IndexType matIdx = i*n + row;
            const IndexType rankIdx = indices[matIdx];
            const OriginalPointerType matValue = valsSplit.read<SegmentsToUse>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.read<SegmentsToUse>(rankIdx);
            localRes += matValue * rankValue;
        }

        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * localRes;
        newRankSplit.write<SegmentsToUse>(row, newRank);
    }
}

// PageRank as Described by matlab: https://www.mathworks.com/help/matlab/examples/use-page-rank-algorithm-to-rank-websites.html
template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void ellPageRank(IndexType n, IndexType k, const IndexType * __restrict__ indices,
    const SplitPointer<NumberSegments, OriginalPointerType> valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPointer<NumberSegments, OriginalPointerType> oldRankSplit,
   SplitPointer<NumberSegments, OriginalPointerType> newRankSplit)
{
    const IndexType row = blockDim.x * blockIdx.x + threadIdx.x;
    if (row < n) {
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        for (IndexType i = 0; i < k; ++i) {
            const IndexType matIdx = i*n + row;
            const IndexType rankIdx = indices[matIdx];
            const OriginalPointerType matValue = valsSplit.read<partPointersToRead>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.read<partPointersToRead>(rankIdx);
            localRes += matValue * rankValue;
        }

        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * (localRes + correction);
        newRankSplit.write<partPointersToWrite>(row, newRank);
    }
}

// PageRank as Described by matlab: https://www.mathworks.com/help/matlab/examples/use-page-rank-algorithm-to-rank-websites.html
template<int ThreadsPerRow, MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void csrPageRank(IndexType m, IndexType n, const IndexType * __restrict__ colIndices,
    const IndexType *__restrict__ rowStartIndices,
    const SplitPointer<NumberSegments, OriginalPointerType> valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPointer<NumberSegments, OriginalPointerType> oldRankSplit,
   SplitPointer<NumberSegments, OriginalPointerType> newRankSplit)
{
    const IndexType idx = blockDim.x * blockIdx.x + threadIdx.x;
    const IndexType row = idx / ThreadsPerRow;
    const IndexType idxShift = idx % ThreadsPerRow;
    if (row < m) {
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        for (IndexType i = rowStartIndices[row] + idxShift; i < rowStartIndices[row + 1]; i+= ThreadsPerRow) {
            const IndexType matIdx = i;
            const IndexType rankIdx = colIndices[i];
            const OriginalPointerType matValue = valsSplit.read<partPointersToRead>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.read<partPointersToRead>(rankIdx);
            localRes += matValue * rankValue;
        }
        
        #pragma unroll
        for (int i = ThreadsPerRow >> 1; i > 0; i = i >> 1) {
#if __CUDACC_VER_MAJOR__ < 9
            localRes += __shfl_xor(localRes, i, ThreadsPerRow);
#else
            localRes += __shfl_xor_sync(std::uint32_t{0xffffffff}, localRes, i, ThreadsPerRow);
#endif
        }
        if (idxShift == 0) {
            const OriginalPointerType newRank =
                (static_cast<OriginalPointerType>(1) - dampingFactor)
                 / static_cast<OriginalPointerType>(n)
                + dampingFactor * (localRes + correction);
            newRankSplit.write<partPointersToWrite>(row, newRank);
        }
    }
}

template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void normalizeVector(IndexType n,
       SplitPointer<NumberSegments, OriginalPointerType> vecValSplit,
       OriginalPointerType norm)
{
    const IndexType idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n) {
        const OriginalPointerType vecValue = vecValSplit.read<SegmentsToUse>(idx);
        vecValSplit.write<SegmentsToUse>(idx, vecValue / norm);
    }
}


/* Thread assignment should be in a way so one thread block handles an integer amount of memory blocks
   After execution, the segment pointer of vec should only be used as a pointer of the type OriginalPointerType.
*/
template<MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void convertToOriginal(IndexType n, SplitPointer<NumberSegments, OriginalPointerType> vec)
{
    const IndexType idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (NumberSegments <= 1) {
        return;
    }
    if (idx < n) {
        OriginalPointerType val = vec.readAll(idx);
        OriginalPointerType *convertedPtr = reinterpret_cast<OriginalPointerType *>(vec.segmentPtr);
        __syncthreads();
        convertedPtr[idx] = val;
    }
}


// Thread assignment should be in a way so one thread block handles an integer amount of memory blocks
template<MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void convertToSplitPointer(IndexType n, SplitPointer<NumberSegments, OriginalPointerType> vec)
{
    const IndexType idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (NumberSegments <= 1) {
        return;
    }
    if (idx < n) {
        OriginalPointerType *convertedPtr = reinterpret_cast<OriginalPointerType *>(vec.segmentPtr);
        OriginalPointerType val = convertedPtr[idx];
        __syncthreads();
        vec.writeAll(idx, val);
    }
}

template<typename FloatType>
__device__ __forceinline__ FloatType myAbs(FloatType value) {
    return abs(value);
}

template<>
__device__ __forceinline__ double myAbs<double>(double value) {
    return fabs(value);
}

template<>
__device__ __forceinline__ float myAbs<float>(float value) {
    return fabsf(value);
}

template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
__global__ void testForPrecision(
        IndexType n,
        OriginalPointerType relPrecisionEpsilon,
        OriginalPointerType relConvergenceEpsilon,
        SplitPointer<NumberSegments, OriginalPointerType> oldVec,
        SplitPointer<NumberSegments, OriginalPointerType> newVec,
        bool *changePrecision,
        bool *isConverged)
{
    const IndexType idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n) {
        const OriginalPointerType oldVal = oldVec.read<SegmentsToUse>(idx);
        const OriginalPointerType newVal = newVec.read<SegmentsToUse>(idx);
        
        const OriginalPointerType refVal = (oldVal < newVal) ? newVal : oldVal;
        const OriginalPointerType relDiff = myAbs(oldVal - newVal) / refVal;
        
        if (relDiff <= relPrecisionEpsilon) {
        //if (relDiff == relPrecisionEpsilon) {
            *changePrecision = true;
        }
        if (relConvergenceEpsilon <= relDiff) {
            *isConverged = false;
        }
    }
}



#endif //GPU_KERNELS_CUH
