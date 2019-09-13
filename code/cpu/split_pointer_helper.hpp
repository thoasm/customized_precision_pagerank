#ifndef SPLIT_POINTER_HELPER_HPP
#define SPLIT_POINTER_HELPER_HPP

#include "split_pointer.hpp"
#include "cpu_kernels.hpp"
#include "macros.hpp"

#include <cstdlib>
#include <cstddef>  // for std::size_t
#include <type_traits>
#include <algorithm>
#include <cstring> // std::memcpy

template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
void newSplitPtr(SplitPointer<NumberSegments, OriginalPointerType> *splitPtr, std::size_t size) {
    //auto numSegs = SplitPointer<NumberSegments, OriginalPointerType>::getMinimumNumberSegmentsFor(size);   // TODO find reason why more space is needed here!
    auto numBytes = SplitPointer<NumberSegments, OriginalPointerType>::getMinimumBytesFor(size);
    //splitPtr->segmentPtr = new typename SplitPointer<NumberSegments, OriginalPointerType>::SegmentType [numSegs];
    //splitPtr->size = numSegs;
    auto errCodeMemAlign = posix_memalign(reinterpret_cast<void **>(&splitPtr->segmentPtr), memoryBlockSize, numBytes);
    if (errCodeMemAlign) {
        throw "Aligned memory allocation failed!";
    }
}

template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
void deleteSplitPtr(SplitPointer<NumberSegments, OriginalPointerType> &splitPtr) {
    //delete[] splitPtr.segmentPtr;
    free(splitPtr.segmentPtr);
    splitPtr.segmentPtr = nullptr;
}

template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
void copy(SplitPointer<NumberSegments, OriginalPointerType> &dest, const SplitPointer<NumberSegments, OriginalPointerType> &src, std::size_t size) {
    std::memcpy(dest.segmentPtr, src.segmentPtr, src.getMinimumBytesFor(size));
}


template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType norm1(
    const SplitPointer<NumberSegments, OriginalPointerType> &vec, IndexType size)
{
    return kernels::norm1<SegmentsToUse, NumberSegments, OriginalPointerType, IndexType>(size, vec);
}


namespace detail {
namespace {

template<MySegmentsNumberType NumberSegments, int SegmentsToUse, typename FloatType, typename... InferredArgs>
void norm1SelectorHelper(
    macros::compile_int_list<SegmentsToUse>,
    FloatType *res,
    InferredArgs... args)
{
    constexpr int segmentsToUse_ = (SegmentsToUse < NumberSegments) ? SegmentsToUse : NumberSegments;
    *res = norm1<segmentsToUse_>(std::forward<InferredArgs>(args)...);
}

GKO_ENABLE_IMPLEMENTATION_SELECTION(norm1Selector, norm1SelectorHelper)

}   // anonymout namespace
}   // detail namespace


template<int NumberSegments, typename FloatType, typename... InferredArgs>
void norm1Selector(FloatType *result, int readParameter, InferredArgs... args)
{
    detail::norm1Selector(
            macros::compile_int_list<1,2,3,4,5,6,7,8>{},    // TODO: Change to range
            [readParameter](int readParameter_) { return readParameter_ == readParameter; },
            macros::compile_int_list<NumberSegments>{},
            macros::compile_type_list<>{},
            result,
            std::forward<InferredArgs>(args)...);
}


template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType norm1Diff(
    const SplitPointer<NumberSegments, OriginalPointerType> &vec1,
    const SplitPointer<NumberSegments, OriginalPointerType> &vec2,
    IndexType size)
{
    return kernels::norm1Diff<SegmentsToUse>(size, vec1, vec2);
}



namespace detail {
namespace {

template<MySegmentsNumberType NumberSegments, int SegmentsToUse, typename FloatType, typename... InferredArgs>
void norm1DiffSelectorHelper(
    macros::compile_int_list<SegmentsToUse>,
    FloatType *res,
    InferredArgs... args)
{
    constexpr int segmentsToUse_ = (SegmentsToUse < NumberSegments) ? SegmentsToUse : NumberSegments;
    *res = norm1Diff<segmentsToUse_>(std::forward<InferredArgs>(args)...);
}

GKO_ENABLE_IMPLEMENTATION_SELECTION(norm1DiffSelector, norm1DiffSelectorHelper)

}   // anonymout namespace
}   // detail namespace


template<int NumberSegments, typename FloatType, typename... InferredArgs>
void norm1DiffSelector(FloatType *result, int readParameter, InferredArgs... args)
{
    detail::norm1DiffSelector(
            macros::compile_int_list<1,2,3,4,5,6,7,8>{},    // TODO: Change to range
            [readParameter](int readParameter_) { return readParameter_ == readParameter; },
            macros::compile_int_list<NumberSegments>{},
            macros::compile_type_list<>{},
            result,
            std::forward<InferredArgs>(args)...);
}

template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
OriginalPointerType selectiveNorm1(
    const SplitPointer<NumberSegments, OriginalPointerType> &vec,
    IndexType indexVectorSize,
    const IndexType *vecIndices)
{
    return kernels::selectiveNorm1<SegmentsToUse>(indexVectorSize, vecIndices, vec);
}

// ------------------------------------------------------------------

namespace detail {
namespace {

template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void normalizeVectorHelper(
    macros::compile_int_list<SegmentsToUse>,
    const SplitPointer<NumberSegments, OriginalPointerType> &vec, IndexType size)
{
    constexpr MySegmentsNumberType segmentsToUse_ = (SegmentsToUse < NumberSegments) ? SegmentsToUse : NumberSegments;
    OriginalPointerType norm = norm1<segmentsToUse_>(vec, size);
    
    kernels::normalizeVector<segmentsToUse_>(size, vec, norm);
}

GKO_ENABLE_IMPLEMENTATION_SELECTION(normalizeVectorSelector, normalizeVectorHelper)

}   // anonymout namespace
}   // detail namespace


template<typename... InferredArgs>
void normalizeVectorSelector(int readParameter, InferredArgs... args)
{
    detail::normalizeVectorSelector(
            macros::compile_int_list<1,2,3,4,5,6,7,8>{},    // TODO: Change to range
            [readParameter](int readParameter_) { return readParameter_ == readParameter; },
            macros::compile_int_list<>{},
            macros::compile_type_list<>{},
            std::forward<InferredArgs>(args)...);
}


#endif  // SPLIT_POINTER_HELPER_HPP
