#ifndef CPU_KERNELS_HPP
#define CPU_KERNELS_HPP

#include "split_pointer.hpp"
#include "omp_helper.hpp"
#include "cpu_kernel_selector.hpp"

namespace kernels {



template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void rwTest(IndexType n, SplitPointer<NumberSegments, OriginalPointerType> rw, OriginalPointerType factor)
{
    detail::KernelSelector<NumberSegments, OriginalPointerType>::template rwTest<SegmentsToUse, IndexType>(n, rw, factor);
}


template<MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void setVector(IndexType n, SplitPointer<NumberSegments, OriginalPointerType> vec, OriginalPointerType initData)
{
    detail::KernelSelector<NumberSegments, OriginalPointerType>::template setVector<IndexType>(n, vec, initData);
}

template<typename ValueType, typename IndexType>
void setVector(IndexType n, ValueType *__restrict__ vec, ValueType initData)
{
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(dynamic, ompHelper::chunkSize)
    #endif
    for (IndexType idx = 0; idx < n; ++idx) {
        vec[idx] = initData;
    }
}


template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType norm1(IndexType n, const SplitPointer<NumberSegments, OriginalPointerType> in)
{
    return detail::KernelSelector<NumberSegments, OriginalPointerType>::template norm1<SegmentsToUse, IndexType>(n, in);
}


template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType norm1Diff(IndexType n, const SplitPointer<NumberSegments, OriginalPointerType> in1,
    const SplitPointer<NumberSegments, OriginalPointerType> in2)
{
    return detail::KernelSelector<NumberSegments, OriginalPointerType>::template norm1Diff<SegmentsToUse, IndexType>(n, in1, in2);
}

template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType selectiveNorm1(
    IndexType indexVectorSize,
    const IndexType * __restrict__ vecIndices,
    const SplitPointer<NumberSegments, OriginalPointerType> vecSplit)
{
    return detail::KernelSelector<NumberSegments, OriginalPointerType>::template selectiveNorm1<SegmentsToUse, IndexType>(indexVectorSize, vecIndices, vecSplit);
}


// PageRank as Described by matlab: https://www.mathworks.com/help/matlab/examples/use-page-rank-algorithm-to-rank-websites.html
template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void ellPageRank(IndexType n, IndexType k, const IndexType * __restrict__ indices,
    const SplitPointer<NumberSegments, OriginalPointerType> valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPointer<NumberSegments, OriginalPointerType> oldRankSplit,
   SplitPointer<NumberSegments, OriginalPointerType> newRankSplit)
{
    detail::KernelSelector<NumberSegments, OriginalPointerType>::template ellPageRank<partPointersToRead, partPointersToWrite, IndexType>(n, k, indices, valsSplit, dampingFactor, correction, oldRankSplit, newRankSplit);
}

// PageRank as Described by matlab: https://www.mathworks.com/help/matlab/examples/use-page-rank-algorithm-to-rank-websites.html
template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void csrPageRank(IndexType m, IndexType n, const IndexType * __restrict__ colIndices,
    const IndexType *__restrict__ rowStartIndices,
    const SplitPointer<NumberSegments, OriginalPointerType> valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPointer<NumberSegments, OriginalPointerType> oldRankSplit,
   SplitPointer<NumberSegments, OriginalPointerType> newRankSplit)
{
    detail::KernelSelector<NumberSegments, OriginalPointerType>::template csrPageRank<partPointersToRead, partPointersToWrite, IndexType>(m, n, colIndices, rowStartIndices, valsSplit, dampingFactor, correction, oldRankSplit, newRankSplit);
}

template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void normalizeVector(IndexType n,
       SplitPointer<NumberSegments, OriginalPointerType> vecValSplit,
       OriginalPointerType norm)
{
    detail::KernelSelector<NumberSegments, OriginalPointerType>::template normalizeVector<SegmentsToUse, IndexType>(n, vecValSplit, norm);
}


/* Thread assignment should be in a way so one thread block handles an integer amount of memory blocks
   After execution, the segment pointer of vec should only be used as a pointer of the type OriginalPointerType.
*/
template<MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void convertToOriginal(IndexType n, SplitPointer<NumberSegments, OriginalPointerType> vec)
{
    detail::KernelSelector<NumberSegments, OriginalPointerType>::template convertToOriginal<IndexType>(n, vec);
}


// Thread assignment should be in a way so one thread block handles an integer amount of memory blocks
template<MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void convertToSplitPointer(IndexType n, SplitPointer<NumberSegments, OriginalPointerType> vec)
{
    detail::KernelSelector<NumberSegments, OriginalPointerType>::template convertToSplitPointer<IndexType>(n, vec);
}
/*
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
void testForPrecision(
        IndexType n,
        OriginalPointerType relPrecisionEpsilon,
        OriginalPointerType relConvergenceEpsilon,
        SplitPointer<NumberSegments, OriginalPointerType> oldVec,
        SplitPointer<NumberSegments, OriginalPointerType> newVec,
        bool *changePrecision,
        bool *isConverged)
{
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(dynamic, ompHelper::chunkSize)
    #endif
    for (std::size_t idx = 0; idx < n; ++idx) {
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
*/

}

#endif //CPU_KERNELS_HPP
