#ifndef SPLIT_POINTER_CUH
#define SPLIT_POINTER_CUH

#include "macros.hpp"
#include "bit_level_functions.hpp"
#include "types.hpp"

#include <cstdint>
#include <type_traits>

#include <iostream>

/*
Names for clarification:
    Segment ^= part of the original type
    memoryBlock ^= Block of memory that is read at once (currently focus on gpu)
*/

using MySegmentsNumberType = int;  //std::uint8_t;

//MUST be a power of two
constexpr std::size_t memoryBlockSize = 8 * 1024; // in Bytes

static_assert(1 << bitLevelFunctions::getTwoExponent<std::size_t, memoryBlockSize>() == memoryBlockSize, "memoryBlockSize MUST be a power of 2!");


//needs to be allocated with at least NumberSegments * memoryBlockSize Bytes
template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
struct SplitPointer {
    static_assert(NumberSegments > 0, "There needs to be at least one segment!");
    static_assert(sizeof(OriginalPointerType) >= NumberSegments, "Out of efficiency, only splitting into bytes supported!");
    static_assert(sizeof(OriginalPointerType) % NumberSegments == 0, "OriginalPointerType needs to be dividable by the number of splits!");
    static_assert(sizeof(OriginalPointerType) <= memoryBlockSize, "Original type is too big, it does not fit into one segment!");
    static_assert(memoryBlockSize % (sizeof(OriginalPointerType) / NumberSegments) == 0, "The segemnts need to fit into one memory block!");
    
    using OriginalType = OriginalPointerType;
    using SegmentType = typename type_of_size<sizeof(OriginalPointerType) / NumberSegments>::type;
    using IntType = typename type_of_size<sizeof(OriginalType)>::type;
    
    static constexpr MySegmentsNumberType numberSegments = NumberSegments;
    //MUST be a power of two
    static constexpr std::size_t segmentsPerMemoryBlock
                = memoryBlockSize * NumberSegments / (sizeof(OriginalPointerType));
    static constexpr std::size_t minimumSize = NumberSegments * memoryBlockSize;   // size in bytes
    static constexpr std::size_t getMinimumBytesFor(std::size_t numElements) {
        return (numElements * sizeof(OriginalPointerType)) + (minimumSize - (numElements * sizeof(OriginalPointerType)) % minimumSize);
    }
    
    static constexpr MY_ATTRIBUTES std::size_t getMinimumNumberSegmentsFor(std::size_t numElements) {
        return getMinimumBytesFor(numElements) / sizeof(SegmentType);
        //return numElements * NumberSegments + NumberSegments * segmentsPerMemoryBlock - (numElements * NumberSegments) % (NumberSegments * segmentsPerMemoryBlock);
    }
    
    template<typename IndexType>
    static constexpr MY_ATTRIBUTES IndexType getHeadIndex(IndexType index) {
       return ((index
                & (~bitLevelFunctions::createMask<IndexType, 
                    bitLevelFunctions::getTwoExponent<IndexType, segmentsPerMemoryBlock>()>()) )
                << bitLevelFunctions::getTwoExponent<IndexType, NumberSegments>())
            + (index 
                    & bitLevelFunctions::createMask<IndexType,
                          bitLevelFunctions::getTwoExponent<IndexType, segmentsPerMemoryBlock>()>())
            ;
    }
    
    SegmentType *segmentPtr;
    
    
    template <MySegmentsNumberType SegmentsToRead, typename IndexType>
    MY_ATTRIBUTES MY_INLINE OriginalPointerType internalRead(IndexType idx) const;
    
    template <MySegmentsNumberType SegmentsToWrite, typename IndexType>
    MY_ATTRIBUTES MY_INLINE void internalWrite(IndexType idx, OriginalPointerType value);
    
    
    
    template <MySegmentsNumberType SegmentsToRead, typename IndexType>
    MY_ATTRIBUTES MY_INLINE OriginalPointerType read(IndexType idx) const {
        return internalRead<SegmentsToRead>(getHeadIndex(idx));
    }
    
    template <MySegmentsNumberType SegmentsToWrite, typename IndexType>
    MY_ATTRIBUTES MY_INLINE void write(IndexType idx, OriginalPointerType value) {
        internalWrite<SegmentsToWrite>(getHeadIndex(idx), value);
    }
    
    template <typename IndexType>
    MY_ATTRIBUTES MY_INLINE OriginalPointerType readAll(IndexType idx) const {
        return read<NumberSegments>(idx);
    }
    
    template <typename IndexType>
    MY_ATTRIBUTES MY_INLINE void writeAll(IndexType idx, OriginalPointerType value) {
        write<NumberSegments>(idx, value);
    }
/*
protected:
    union Conv {
        OriginalType flt;
        IntType it;
    };*/
};


template<typename OriginalPointerType>
struct SplitPointer<1, OriginalPointerType> {
    
    using OriginalType = OriginalPointerType;
    using SegmentType = OriginalPointerType;
    using IntType = typename type_of_size<sizeof(OriginalType)>::type;
    
    //static constexpr IntType lowMask = bitLevelFunctions::createMask<IntType, sizeof(SegmentType) * 8>();
    
    static constexpr MySegmentsNumberType numberSegments = 1;
    static constexpr std::uint32_t segmentsPerMemoryBlock = 1;
    static constexpr std::size_t minimumSize = sizeof(OriginalPointerType);   // size in bytes
    static constexpr std::size_t getMinimumBytesFor(std::size_t numElements) {
        return numElements * sizeof(OriginalPointerType);
    }
    static constexpr MY_ATTRIBUTES std::size_t getMinimumNumberSegmentsFor(std::size_t numElements) {
        return numElements;
    }
    
    template<typename IndexType>
    static constexpr MY_ATTRIBUTES IndexType getHeadIndex(IndexType index) {
        return index;
    }
    
    SegmentType *segmentPtr;
    //std::size_t size;
    
    template <MySegmentsNumberType , typename IndexType>
    MY_ATTRIBUTES MY_INLINE OriginalPointerType internalRead(IndexType idx) const {
        return segmentPtr[idx];
    }
    
    template <MySegmentsNumberType , typename IndexType>
    MY_ATTRIBUTES MY_INLINE void internalWrite(IndexType idx, OriginalPointerType value) {
        segmentPtr[idx] = value;
    }
    
    template <MySegmentsNumberType SegmentsToRead, typename IndexType>
    MY_ATTRIBUTES MY_INLINE OriginalPointerType read(IndexType idx) const {
        return internalRead<SegmentsToRead>(getHeadIndex(idx));
    }
    
    template <MySegmentsNumberType SegmentsToWrite, typename IndexType>
    MY_ATTRIBUTES MY_INLINE void write(IndexType idx, OriginalPointerType value) {
        internalWrite<SegmentsToWrite>(getHeadIndex(idx), value);
    }
    
    template <typename IndexType>
    MY_ATTRIBUTES MY_INLINE OriginalPointerType readAll(IndexType idx) const {
        return read<1>(idx);
    }
    
    template <typename IndexType>
    MY_ATTRIBUTES MY_INLINE void writeAll(IndexType idx, OriginalPointerType value) {
        write<1>(idx, value);
    }
};



namespace detail {
namespace {


template<typename ValueIntType, MySegmentsNumberType SegmentsToCopy, MySegmentsNumberType CurrentIteration = SegmentsToCopy>
struct CopyHelper {
    static_assert(CurrentIteration <= SegmentsToCopy, "Iteration must be lower than total iteration count!");
    
    
    template<std::size_t SegmentsPerMemoryBlock, MySegmentsNumberType NumberSegmentsPerValue, typename SegmentPointerType>
    static MY_ATTRIBUTES MY_INLINE ValueIntType read(const SegmentPointerType *MY_RESTRICT src) {
        static_assert(1 <= SegmentsToCopy && SegmentsToCopy <= NumberSegmentsPerValue, "SegmentsToCopy must be between 1 and NumberSegmentsPerValue!");
        static_assert(NumberSegmentsPerValue > 0, "Number of segments has to be > 0 !");
        
        //numerated from most significant (NumberSegmentsPerValue, left) to least significant (1, right)
        constexpr int currentSegment = NumberSegmentsPerValue - (SegmentsToCopy - CurrentIteration);
        constexpr int bitsPerSegment = sizeof(SegmentPointerType) * 8;
        constexpr ValueIntType lowMask = bitLevelFunctions::createMask<ValueIntType, bitsPerSegment>();
        constexpr int shiftNumber = bitsPerSegment * (currentSegment - 1);
        
        ValueIntType curVal = (static_cast<ValueIntType>(*src) & lowMask) << shiftNumber;
        
        return curVal | 
                CopyHelper<ValueIntType, SegmentsToCopy, CurrentIteration - 1>::template read<SegmentsPerMemoryBlock, NumberSegmentsPerValue>(src + SegmentsPerMemoryBlock);
    }
    
    
    template<std::size_t SegmentsPerMemoryBlock, MySegmentsNumberType NumberSegmentsPerValue, typename SegmentPointerType>
    static MY_ATTRIBUTES MY_INLINE void write(ValueIntType src, SegmentPointerType *MY_RESTRICT dest) {
        static_assert(1 <= SegmentsToCopy && SegmentsToCopy <= NumberSegmentsPerValue, "SegmentsToCopy must be between 1 and NumberSegmentsPerValue!");
        static_assert(NumberSegmentsPerValue > 0, "Number of segments has to be > 0 !");
        
        //numerated from most significant (NumberSegmentsPerValue, left) to least significant (1, right)
        constexpr int currentSegment = NumberSegmentsPerValue - (SegmentsToCopy - CurrentIteration);
        constexpr int bitsPerSegment = sizeof(SegmentPointerType) * 8;
        constexpr ValueIntType lowMask = bitLevelFunctions::createMask<ValueIntType, bitsPerSegment>();
        constexpr int shiftNumber = bitsPerSegment * (currentSegment - 1);
        
        ValueIntType curVal = (src >> shiftNumber) & lowMask;
        
        *dest = curVal;
        CopyHelper<ValueIntType, SegmentsToCopy, CurrentIteration - 1>::template write<SegmentsPerMemoryBlock, NumberSegmentsPerValue>(src, dest + SegmentsPerMemoryBlock);
    }
};

//end of recursion
template<typename ValueIntType, MySegmentsNumberType SegmentsToCopy>
struct CopyHelper<ValueIntType, SegmentsToCopy, 0> {
    template<std::size_t SegmentsPerMemoryBlock, MySegmentsNumberType NumberSegmentsPerValue, typename SegmentPointerType>
    static MY_ATTRIBUTES MY_INLINE ValueIntType read(const SegmentPointerType *MY_RESTRICT src) {
        return 0;
    }
    template<std::size_t SegmentsPerMemoryBlock, MySegmentsNumberType NumberSegmentsPerValue, typename SegmentPointerType>
    static MY_ATTRIBUTES MY_INLINE void write(ValueIntType src, SegmentPointerType *MY_RESTRICT dest) {
        return;
    }
};

template<typename SourceType, typename DestType>
MY_INLINE MY_ATTRIBUTES DestType convert(SourceType source) {
    union { SourceType src; DestType dest;  };
    src = source;
    return dest;
}

/*
template<>
MY_INLINE double convert<std::int64_t, double>(std::int64_t source) {
    double ret;
    __asm__(""//"movq %[in_int], %[out_dbl]"
              : [out_dbl] "=r" (ret)
              : [in_int] "r" (source)
              //: ""
              );
    return ret;
}

template<>
MY_INLINE std::int64_t convert<double, std::int64_t>(double source) {
    std::int64_t ret;
    __asm__(""//"movq %[in_int], %[out_dbl]"
              : [out_dbl] "=r" (ret)
              : [in_int] "r" (source)
              //: ""
              );
    return ret;

}
*/
} // unnamed namespace
} // detail namespace

template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
template <MySegmentsNumberType SegmentsToRead, typename IndexType>
MY_ATTRIBUTES MY_INLINE OriginalPointerType SplitPointer<NumberSegments, OriginalPointerType>::internalRead(IndexType idx) const
{
    IntType dest;
    dest = 
        detail::CopyHelper<IntType, SegmentsToRead>::template read<segmentsPerMemoryBlock, NumberSegments>(segmentPtr + idx);
    
    return detail::convert<IntType, OriginalPointerType>(dest);
}


template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
template <MySegmentsNumberType SegmentsToWrite, typename IndexType>
MY_ATTRIBUTES MY_INLINE void SplitPointer<NumberSegments, OriginalPointerType>::internalWrite(IndexType idx, OriginalPointerType src)
{
    //Conv source;
    IntType source = detail::convert<OriginalPointerType, IntType>(src);
    
    detail::CopyHelper<IntType, SegmentsToWrite>::template write<segmentsPerMemoryBlock, NumberSegments>(source, segmentPtr + idx);
}


#endif
