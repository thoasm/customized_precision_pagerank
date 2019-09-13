#ifndef AVX2_HELPER_HPP
#define AVX2_HELPER_HPP

#include "bit_level_functions.hpp"
#include "split_pointer.hpp"

#include <immintrin.h>
#include <emmintrin.h>

#include <limits>

namespace avx2Helper {
/*
void printIdxs(__m128i indices32) {
    int32_t *pIdx = reinterpret_cast<int32_t *>(&indices32);
    for (std::size_t j = 0; j < 4; ++j) {
        std::cout << pIdx[j] << " ";
    } std::cout << std::endl;
}
*/
template<std::size_t SegmentsToUse, std::size_t NumberSegments, typename FloatType>
inline __m256d readPack(const SplitPointer<NumberSegments, FloatType> &val, __m128i indices32) {
    union {__m128i indices32_; int indices[4];};
    indices32_ = indices32;
    
    return _mm256_setr_pd(val.template read<SegmentsToUse>(indices[0]),
                          val.template read<SegmentsToUse>(indices[1]),
                          val.template read<SegmentsToUse>(indices[2]),
                          val.template read<SegmentsToUse>(indices[3])
                          );
    //throw "This is only implemented for double and (2 or 1) segment(s)!\n";
}

template<std::size_t SegmentsToUse, std::size_t NumberSegments, typename FloatType>
inline __m256d directReadPack(const SplitPointer<NumberSegments, FloatType> &val, int directHeadStartIdx) {
    return _mm256_setr_pd(val.template internalRead<SegmentsToUse>(directHeadStartIdx),
                          val.template internalRead<SegmentsToUse>(directHeadStartIdx + 1),
                          val.template internalRead<SegmentsToUse>(directHeadStartIdx + 2),
                          val.template internalRead<SegmentsToUse>(directHeadStartIdx + 3)
                          );
    //throw "This is only implemented for double and (2 or 1) segment(s)!\n";
}


//It should never be written over a block-boarder!!!
template<std::size_t SegmentsToUse, std::size_t NumberSegments, typename FloatType>
void writePack(SplitPointer<NumberSegments, FloatType> &var, __m256d dblPack, int startIndex) {
    union {__m256d dblPack_; double dblValue[4];};
    dblPack_ = dblPack;
    
    var.template write<SegmentsToUse>(startIndex, dblValue[0]);
    var.template write<SegmentsToUse>(startIndex + 1, dblValue[1]);
    var.template write<SegmentsToUse>(startIndex + 2, dblValue[2]);
    var.template write<SegmentsToUse>(startIndex + 3, dblValue[3]);
    //throw "This is only implemented for double and (2 or 1) segment(s)!\n";
}

//It should never be written over a block-boarder!!!
template<std::size_t SegmentsToUse, std::size_t NumberSegments, typename FloatType>
void directWritePack(SplitPointer<NumberSegments, FloatType> &var, __m256d dblPack, int directHeadStartIdx) {
    union {__m256d dblPack_; double dblValue[4];};
    dblPack_ = dblPack;
    
    var.template internalWrite<SegmentsToUse>(directHeadStartIdx, dblValue[0]);
    var.template internalWrite<SegmentsToUse>(directHeadStartIdx + 1, dblValue[1]);
    var.template internalWrite<SegmentsToUse>(directHeadStartIdx + 2, dblValue[2]);
    var.template internalWrite<SegmentsToUse>(directHeadStartIdx + 3, dblValue[3]);
    //throw "This is only implemented for double and (2 or 1) segment(s)!\n";
}

template<std::size_t SegmentsToUse, std::size_t NumberSegments, typename FloatType>
inline __m256d readPack(const SplitPointer<NumberSegments, FloatType> &val, int startIndex) {
    constexpr int localIndexMask = bitLevelFunctions::createMask<int, 
                    bitLevelFunctions::getTwoExponent<int, 
                        SplitPointer<NumberSegments, FloatType>::segmentsPerMemoryBlock>()>();
    constexpr int blockMask = ~localIndexMask;
    
    const int firstHeadStartIdx = (startIndex & blockMask) 
        << bitLevelFunctions::template getTwoExponent<int, NumberSegments>();
    const int localIdx = startIndex & localIndexMask;
    
    if (localIdx < SplitPointer<NumberSegments, FloatType>::segmentsPerMemoryBlock - 3) {
        return directReadPack<SegmentsToUse, NumberSegments, FloatType>(val, localIdx + firstHeadStartIdx);
    }
    else {
        __m128i indices32 = _mm_setr_epi32(startIndex, startIndex+1, startIndex+2, startIndex+3);
        return readPack<SegmentsToUse, NumberSegments, FloatType>(val, indices32);
    }
}

//calculates the absolute value of a packed double precision vector by removing the sign bit on bit level
inline __m256d absPack(__m256d val) {
    return _mm256_andnot_pd(_mm256_set1_pd(-0.0), val);
}
        

template<>
inline __m256d directReadPack<1,1, double>(const SplitPointer<1, double> &val, int directHeadStartIdx)
{
    return _mm256_loadu_pd(val.segmentPtr + directHeadStartIdx);
}


template<>
inline __m256d readPack<1,1, double>(const SplitPointer<1, double> &val, int startIndex)
{
    return directReadPack<1,1,double>(val, startIndex);
}

template<>
inline __m256d readPack<1,1, double>(const SplitPointer<1, double> &val, __m128i indices32)
{
    return _mm256_i32gather_pd(val.segmentPtr, indices32, 8);
}


template<>
inline void directWritePack<1,1, double>(SplitPointer<1, double> &var, __m256d dblPack, int directHeadStartIdx) {
    _mm256_storeu_pd(var.segmentPtr + directHeadStartIdx, dblPack);
}

template<>
inline void writePack<1,1, double>(SplitPointer<1, double> &var, __m256d dblPack, int startIndex) {
    directWritePack<1,1,double>(var, dblPack, startIndex);
}


namespace detail {


template<MySegmentsNumberType numberSegments, std::size_t segmentsPerMemoryBlock>
inline __m128i recalculateIndices(__m128i indices32) {
    constexpr int numberRightShifts = bitLevelFunctions::getTwoExponent<int, segmentsPerMemoryBlock>();
    __m128i blockNumberPack = _mm_srli_epi32(indices32, numberRightShifts);
    __m128i numberSegmentsUntilRepeatPack 
            = _mm_set1_epi32(numberSegments * segmentsPerMemoryBlock);
    __m128i blockStartPack = _mm_mullo_epi32(blockNumberPack, numberSegmentsUntilRepeatPack);
    
    constexpr int indexMask = bitLevelFunctions::createMask<int, numberRightShifts>();
    __m128i localIndexMaskPack = _mm_set1_epi32(indexMask);
    __m128i localIndexPack = _mm_and_si128(indices32, localIndexMaskPack);
    
    __m128i actualHeadIndexPack = _mm_add_epi32(blockStartPack, localIndexPack);
    return actualHeadIndexPack;
}

template<MySegmentsNumberType numberSegments, std::size_t segmentsPerMemoryBlock>
inline int recalculateIndices(int index) {
    constexpr int numberRightShifts = bitLevelFunctions::getTwoExponent<int, segmentsPerMemoryBlock>();
    const int blockNumber = index >> numberRightShifts;
    const int numberSegmentsUntilRepeat = numberSegments * segmentsPerMemoryBlock;
    const int blockStart = blockNumber * numberSegmentsUntilRepeat;
    
    constexpr int indexMask = bitLevelFunctions::createMask<int, numberRightShifts>();
    int localIndex = index & indexMask;
    
    return blockStart + localIndex;
}

}


template<>
inline __m256d readPack<1,2, double>(const SplitPointer<2, double> &val, __m128i indices32)
{
    constexpr MySegmentsNumberType numberSegments = 2;
    constexpr std::size_t segmentsPerMemoryBlock = SplitPointer<numberSegments, double>::segmentsPerMemoryBlock;
    
    __m128i actualHeadIndexPack = detail::recalculateIndices<numberSegments, segmentsPerMemoryBlock>(indices32);
    
    __m128i headPack = _mm_i32gather_epi32(val.segmentPtr, actualHeadIndexPack, 4);
    __m256i exHeadPack = _mm256_cvtepu32_epi64(headPack);
    exHeadPack = _mm256_slli_epi64(exHeadPack, 32);
    //the cast has zero overhead, just there to make the compiler happy
    return _mm256_castsi256_pd(exHeadPack);
}



template<>
inline __m256d readPack<2,2>(const SplitPointer<2, double> &val, __m128i indices32)
{
    constexpr MySegmentsNumberType numberSegments = 2;
    constexpr std::size_t segmentsPerMemoryBlock = SplitPointer<numberSegments, double>::segmentsPerMemoryBlock;
    
    __m128i actualHeadIndexPack = detail::recalculateIndices<numberSegments, segmentsPerMemoryBlock>(indices32);
    __m128i segmentsPerBlockPack = _mm_set1_epi32(val.segmentsPerMemoryBlock);
    __m128i actualTailIndexPack = _mm_add_epi32(actualHeadIndexPack, segmentsPerBlockPack);
    
    __m128i headPack = _mm_i32gather_epi32(val.segmentPtr, actualHeadIndexPack, 4);
    __m256i exHeadPack = _mm256_cvtepu32_epi64(headPack);
    exHeadPack = _mm256_slli_epi64(exHeadPack, 32);
    
    __m128i tailPack = _mm_i32gather_epi32(val.segmentPtr, actualTailIndexPack, 4);
    __m256i exTailPack = _mm256_cvtepu32_epi64(tailPack);

    //the cast has zero overhead, just there to make compiler happy
    return _mm256_castsi256_pd(_mm256_or_si256(exHeadPack, exTailPack));
}


template<>
inline __m256d directReadPack<1,2, double>(const SplitPointer<2, double> &val, int directHeadStartIdx)
{
    constexpr MySegmentsNumberType numberSegments = 2;
    constexpr std::size_t segmentsPerMemoryBlock = SplitPointer<numberSegments, double>::segmentsPerMemoryBlock;
    
    __m128i headPack = _mm_loadu_si128(reinterpret_cast<__m128i *>(val.segmentPtr + directHeadStartIdx));
    __m256i exHeadPack = _mm256_cvtepu32_epi64(headPack);
    exHeadPack = _mm256_slli_epi64(exHeadPack, 32);
    //the cast has zero overhead, just there to make the compiler happy
    return _mm256_castsi256_pd(exHeadPack);
}


template<>
inline __m256d directReadPack<2,2, double>(const SplitPointer<2, double> &val, int directHeadStartIdx)
{
    constexpr MySegmentsNumberType numberSegments = 2;
    constexpr std::size_t segmentsPerMemoryBlock = SplitPointer<numberSegments, double>::segmentsPerMemoryBlock;
    
    __m128i headPack = _mm_loadu_si128(reinterpret_cast<__m128i *>(val.segmentPtr + directHeadStartIdx));
    __m256i exHeadPack = _mm256_cvtepu32_epi64(headPack);
    exHeadPack = _mm256_slli_epi64(exHeadPack, 32);
    //the cast has zero overhead, just there to make the compiler happy
    __m128i tailPack = _mm_loadu_si128(reinterpret_cast<__m128i *>(val.segmentPtr + directHeadStartIdx + segmentsPerMemoryBlock));
    __m256i exTailPack = _mm256_cvtepu32_epi64(tailPack);

    //the cast has zero overhead, just there to make compiler happy
    return _mm256_castsi256_pd(_mm256_or_si256(exHeadPack, exTailPack));
}


template<>
inline __m256d readPack<1,2, double>(const SplitPointer<2, double> &val, int startIdx)
{
    constexpr MySegmentsNumberType numberSegments = 2;
    constexpr int localIndexMask = bitLevelFunctions::createMask<int, 
                    bitLevelFunctions::getTwoExponent<int, 
                        SplitPointer<2, double>::segmentsPerMemoryBlock>()>();
    constexpr int blockMask = ~localIndexMask;
    
    const int firstHeadStartIdx = (startIdx & blockMask) << 1;
    const int localIdx = startIdx & localIndexMask;
    
    if (localIdx < SplitPointer<2, double>::segmentsPerMemoryBlock - 3) {
        return directReadPack<1,2,double>(val, localIdx + firstHeadStartIdx);
    }
    else {
        __m128i indices32 = _mm_setr_epi32(startIdx, startIdx+1, startIdx+2, startIdx+3);
        return readPack<1,2,double>(val, indices32);
    }
}

template<>
inline __m256d readPack<2,2, double>(const SplitPointer<2, double> &val, int startIdx)
{
    constexpr MySegmentsNumberType numberSegments = 2;
    constexpr int localIndexMask = bitLevelFunctions::createMask<int, 
                    bitLevelFunctions::getTwoExponent<int, 
                        SplitPointer<2, double>::segmentsPerMemoryBlock>()>();
    constexpr int blockMask = ~localIndexMask;
    
    const int firstHeadStartIdx = (startIdx & blockMask) << 1;
    const int localIdx = startIdx & localIndexMask;
    
    if (localIdx < SplitPointer<2, double>::segmentsPerMemoryBlock - 3) {
        return directReadPack<2,2,double>(val, localIdx + firstHeadStartIdx);
    }
    else {
        __m128i indices32 = _mm_setr_epi32(startIdx, startIdx+1, startIdx+2, startIdx+3);
        return readPack<2,2,double>(val, indices32);
    }
}


template<>
inline void writePack<1,2>(SplitPointer<2, double> &var, __m256d dblPack, int startIndex) {
    constexpr MySegmentsNumberType numberSegments = 2;
    constexpr std::size_t segmentsPerMemoryBlock = SplitPointer<numberSegments, double>::segmentsPerMemoryBlock;
    __m256i index32Pack = _mm256_setr_epi32(0,2,4,6,1,3,5,7);

    int actualIndex = detail::recalculateIndices<numberSegments, segmentsPerMemoryBlock>(startIndex);
    
    __m256i mixPack = _mm256_permutevar8x32_epi32(_mm256_castpd_si256(dblPack), index32Pack);
    //now: tail | head in mixPack
    __m128i tailPack = _mm256_extracti128_si256(mixPack, 0);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(var.segmentPtr+actualIndex), tailPack);
}

template<>
inline void writePack<2,2>(SplitPointer<2, double> &var, __m256d dblPack, int startIndex) {
    constexpr MySegmentsNumberType numberSegments = 2;
    constexpr std::size_t segmentsPerMemoryBlock = SplitPointer<numberSegments, double>::segmentsPerMemoryBlock;
    __m256i index32Pack = _mm256_setr_epi32(0,2,4,6,1,3,5,7);

    int actualIndex = detail::recalculateIndices<numberSegments, segmentsPerMemoryBlock>(startIndex);

    __m256i mixPack = _mm256_permutevar8x32_epi32(_mm256_castpd_si256(dblPack), index32Pack);
    //now: tail | head in mixPack
    __m128i tailPack = _mm256_extracti128_si256(mixPack, 0);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(var.segmentPtr+actualIndex), tailPack);
    actualIndex += segmentsPerMemoryBlock;
    __m128i headPack = _mm256_extracti128_si256(mixPack, 1);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(var.segmentPtr+actualIndex), headPack);
}


/*
static inline __m256d getMergedMask(const SplitPtr<NumberSegments> &val, int startIndex, __m128i tailMask)
{
    __m128i tailPack = _mm_maskload_epi32((val.ptr[0].data() + startIndex), tailMask);
    __m256i exTailPack = _mm256_cvtepu32_epi64(tailPack);

    __m128i headPack = _mm_loadu_si128((__m128i *)(val.ptr[1].data() + startIndex));
    __m256i exHeadPack = _mm256_cvtepu32_epi64(headPack);
    exHeadPack = _mm256_slli_epi64(exHeadPack, 32);
    //the cast has zero overhead, just there to make compiler happy
    return _mm256_castsi256_pd(_mm256_or_si256(exHeadPack, exTailPack));
}*/

}




#endif
