#ifndef SPLIT_POINTER_HELPER_CUH
#define SPLIT_POINTER_HELPER_CUH

#include "split_pointer.hpp"
#include "gpu_kernels.cuh"
#include "macros.hpp"

#include <cstddef>  // for std::size_t
#include <type_traits>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <algorithm>

template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
cudaError_t mallocSplitPtr(SplitPointer<NumberSegments, OriginalPointerType> *splitPtr, std::size_t size) {
    return cudaMalloc(&(splitPtr->segmentPtr), splitPtr->getMinimumBytesFor(size));
}

template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
cudaError_t freeSplitPtr(SplitPointer<NumberSegments, OriginalPointerType> &splitPtr) {
    auto ret = cudaFree(splitPtr.segmentPtr);
    splitPtr.segmentPtr = nullptr;
    return ret;
}

template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
cudaError_t copy(const SplitPointer<NumberSegments, OriginalPointerType> &dest, const SplitPointer<NumberSegments, OriginalPointerType> &src, std::size_t size, cudaMemcpyKind cudaKind) {
    return cudaMemcpy(dest.segmentPtr, src.segmentPtr, src.getMinimumBytesFor(size), cudaKind);
}


template<unsigned int BlockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType cudaNorm1(
    const SplitPointer<NumberSegments, OriginalPointerType> &d_vec, IndexType size,
    //const SplitPointer<NumberSegments, OriginalPointerType> &d_tempSplit,
    OriginalPointerType *d_temp)
{
    OriginalPointerType *d_result = d_temp;
    dim3 blockNorm1(BlockSize);
    dim3 gridNorm1( std::min((size - 1) / blockNorm1.x + 1, BlockSize) );
    norm1Step1<BlockSize, SegmentsToUse><<<gridNorm1, blockNorm1>>>(size, d_vec, d_result);
    if (gridNorm1.x > 1) {
        d_result = d_temp + gridNorm1.x;
        dim3 blockAddReduction(BlockSize);
        dim3 gridAddReduction(1);
        addReduction<BlockSize><<<gridAddReduction, blockAddReduction>>>(gridNorm1.x, d_temp, d_result);
    }
    OriginalPointerType res;
    cudaMemcpy(&res, d_result, 1 * sizeof(OriginalPointerType), cudaMemcpyDeviceToHost);
    return res;
}


template<unsigned int BlockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType cudaNorm1(
    const SplitPointer<NumberSegments, OriginalPointerType> &d_vec, IndexType size,
    const SplitPointer<NumberSegments, OriginalPointerType> &d_tempSplit,
    OriginalPointerType *d_temp)
{
    dim3 blockNorm1(BlockSize);
    dim3 gridNorm1( std::min((size - 1) / blockNorm1.x + 1, BlockSize) );
    norm1Step1<BlockSize, SegmentsToUse><<<gridNorm1, blockNorm1>>>(size, d_vec, d_tempSplit);
    
    dim3 blockAddReduction(BlockSize);
    dim3 gridAddReduction(1);
    addReduction<BlockSize, SegmentsToUse><<<gridAddReduction, blockAddReduction>>>(gridNorm1.x, d_tempSplit, d_temp);
    
    OriginalPointerType res;
    cudaMemcpy(&res, d_temp, 1 * sizeof(OriginalPointerType), cudaMemcpyDeviceToHost);
    return res;
}


namespace detail {
namespace {

template<unsigned int BlockSize, MySegmentsNumberType NumberSegments, int SegmentsToUse, typename FloatType, typename... InferredArgs>
void cudaNorm1SelectorHelper(
    macros::compile_int_list<SegmentsToUse>,
    FloatType *res,
    InferredArgs... args)
{
    constexpr int segmentsToUse_ = (SegmentsToUse < NumberSegments) ? SegmentsToUse : NumberSegments;
    *res = cudaNorm1<BlockSize, segmentsToUse_>(std::forward<InferredArgs>(args)...);
}

GKO_ENABLE_IMPLEMENTATION_SELECTION(cudaNorm1Selector, cudaNorm1SelectorHelper)

}   // anonymout namespace
}   // detail namespace


template<unsigned int BlockSize, int NumberSegments, typename FloatType, typename... InferredArgs>
void cudaNorm1Selector(FloatType *result, int readParameter, InferredArgs... args)
{
    detail::cudaNorm1Selector(
            macros::compile_int_list<1,2,3,4,5,6,7,8>{},    // TODO: Change to range
            [readParameter](int readParameter_) { return readParameter_ == readParameter; },
            macros::compile_int_list<BlockSize, NumberSegments>{},
            macros::compile_type_list<>{},
            result,
            std::forward<InferredArgs>(args)...);
}


template<unsigned int BlockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType cudaNorm1Diff(
    const SplitPointer<NumberSegments, OriginalPointerType> &d_vec1,
    const SplitPointer<NumberSegments, OriginalPointerType> &d_vec2,
    IndexType size,
    //const SplitPointer<NumberSegments, OriginalPointerType> &d_tempSplit,
    OriginalPointerType *d_temp)
{
    OriginalPointerType *d_result = d_temp;
    dim3 blockNorm1(BlockSize);
    dim3 gridNorm1( std::min((size - 1) / blockNorm1.x + 1, BlockSize) );
    norm1DiffStep1<BlockSize, SegmentsToUse><<<gridNorm1, blockNorm1>>>(size, d_vec1, d_vec2, d_result);
    if (gridNorm1.x > 1) {
        d_result = d_temp + gridNorm1.x;
        dim3 blockAddReduction(BlockSize);
        dim3 gridAddReduction(1);
        addReduction<BlockSize><<<gridAddReduction, blockAddReduction>>>(gridNorm1.x, d_temp, d_result);
    }
    OriginalPointerType res;
    cudaMemcpy(&res, d_result, 1 * sizeof(OriginalPointerType), cudaMemcpyDeviceToHost);
    return res;
}



namespace detail {
namespace {

template<unsigned int BlockSize, MySegmentsNumberType NumberSegments, int SegmentsToUse, typename FloatType, typename... InferredArgs>
void cudaNorm1DiffSelectorHelper(
    macros::compile_int_list<SegmentsToUse>,
    FloatType *res,
    InferredArgs... args)
{
    constexpr int segmentsToUse_ = (SegmentsToUse < NumberSegments) ? SegmentsToUse : NumberSegments;
    *res = cudaNorm1Diff<BlockSize, segmentsToUse_>(std::forward<InferredArgs>(args)...);
}

GKO_ENABLE_IMPLEMENTATION_SELECTION(cudaNorm1DiffSelector, cudaNorm1DiffSelectorHelper)

}   // anonymout namespace
}   // detail namespace


template<unsigned int BlockSize, int NumberSegments, typename FloatType, typename... InferredArgs>
void cudaNorm1DiffSelector(FloatType *result, int readParameter, InferredArgs... args)
{
    detail::cudaNorm1DiffSelector(
            macros::compile_int_list<1,2,3,4,5,6,7,8>{},    // TODO: Change to range
            [readParameter](int readParameter_) { return readParameter_ == readParameter; },
            macros::compile_int_list<BlockSize, NumberSegments>{},
            macros::compile_type_list<>{},
            result,
            std::forward<InferredArgs>(args)...);
}

template<unsigned int BlockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType cudaSelectiveNorm1(
    const SplitPointer<NumberSegments, OriginalPointerType> &d_vec,
    IndexType indexVectorSize,
    const IndexType *d_vecIndices,
    OriginalPointerType *d_temp)
{
    OriginalPointerType *d_result = d_temp;
    dim3 blockNorm1(BlockSize);
    dim3 gridNorm1( std::min((indexVectorSize - 1) / blockNorm1.x + 1, BlockSize) );
    selectiveNorm1Step1<BlockSize, SegmentsToUse><<<gridNorm1, blockNorm1>>>(
        indexVectorSize,
        d_vecIndices,
        d_vec,
        d_result);
    if (gridNorm1.x > 1) {
        d_result = d_temp + gridNorm1.x;
        cudaDeviceSynchronize();
        dim3 blockAddReduction(BlockSize);
        dim3 gridAddReduction(1);
        addReduction<BlockSize><<<gridAddReduction, blockAddReduction>>>(gridNorm1.x, d_temp, d_result);
    }
    OriginalPointerType res;
    cudaMemcpy(&res, d_result, 1 * sizeof(OriginalPointerType), cudaMemcpyDeviceToHost);
    return res;
}

// ------------------------------------------------------------------

namespace detail {
namespace {

template<unsigned int BlockSize, MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void normalizeVectorHelper(
    macros::compile_int_list<SegmentsToUse>,
    const SplitPointer<NumberSegments, OriginalPointerType> &d_vec, IndexType size,
    OriginalPointerType *d_temp)
{
    constexpr MySegmentsNumberType segmentsToUse_ = (SegmentsToUse < NumberSegments) ? SegmentsToUse : NumberSegments;
    OriginalPointerType norm = cudaNorm1<BlockSize, segmentsToUse_>(d_vec, size, d_temp);
    
    dim3 blockNormalize(BlockSize);
    dim3 gridNormalize((size - 1) / blockNormalize.x + 1);
    normalizeVector<segmentsToUse_><<<gridNormalize, blockNormalize>>>(size, d_vec, norm);
}

GKO_ENABLE_IMPLEMENTATION_SELECTION(normalizeVectorSelector, normalizeVectorHelper)

}   // anonymout namespace
}   // detail namespace


template<unsigned int BlockSize, typename... InferredArgs>
void normalizeVectorSelector(int readParameter, InferredArgs... args)
{
    detail::normalizeVectorSelector(
            macros::compile_int_list<1,2,3,4,5,6,7,8>{},    // TODO: Change to range
            [readParameter](int readParameter_) { return readParameter_ == readParameter; },
            macros::compile_int_list<BlockSize>{},
            macros::compile_type_list<>{},
            std::forward<InferredArgs>(args)...);
}


#endif  // SPLIT_POINTER_HELPER_CUH
