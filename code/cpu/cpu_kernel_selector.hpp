#ifndef CPU_KERNEL_SELECTOR_HPP
#define CPU_KERNEL_SELECTOR_HPP

#include "bit_level_functions.hpp"
#include "split_pointer.hpp"
#include "avx2_helper.hpp"
#include "omp_helper.hpp"
#include "thread_buffer.hpp"

#include <cinttypes>
#include <type_traits>
#include <algorithm>
#include <cstdlib>  //for std::abs


#define USE_OWN_CSR_SPMV_VECTORIZATION
#define USE_OWN_ELL_SPMV_VECTORIZATION
#define USE_VECTORIZATION
#define USE_SINGLE_FOR_BLOCK



namespace kernels {

namespace detail {

template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
struct BasicKernels
{

using SplitPtr = SplitPointer<NumberSegments, OriginalPointerType>;

static constexpr std::size_t segmentsPerBlock = SplitPtr::segmentsPerMemoryBlock;

static constexpr std::size_t potentialChunk = ompHelper::blas1ChunkSize / segmentsPerBlock;

static constexpr int blockChunk = (1 < potentialChunk) ? potentialChunk : 1;

template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static void rwTest(IndexType n, SplitPtr rw, OriginalPointerType factor)
{
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(static, blockChunk)
    #endif
    for (IndexType idx = 0; idx < n; idx += segmentsPerBlock) {
        const IndexType headIndexStart = SplitPtr::getHeadIndex(idx);
        const IndexType headIndexStop = headIndexStart + segmentsPerBlock;
        const IndexType headLoopStop = headIndexStop;
                                        //(n - idx < segmentsPerBlock) 
                                        //? headIndexStart + (n - idx)
                                        //: headIndexStop;
        for (IndexType headIdx = headIndexStart; headIdx < headLoopStop; ++headIdx) {
            OriginalPointerType val = rw.template internalRead<SegmentsToUse>(headIdx);
            rw.template internalWrite<SegmentsToUse>(headIdx, val * factor);
        }
    }
}

template<typename IndexType>
static void setVector(IndexType n, SplitPtr vec, OriginalPointerType initData)
{
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(static, blockChunk)
    #endif
    for (IndexType idx = 0; idx < n; idx += segmentsPerBlock) {
        const IndexType headIndexStart = SplitPtr::getHeadIndex(idx);
        const IndexType headIndexStop = headIndexStart + segmentsPerBlock;
        const IndexType headLoopStop = headIndexStop;
                                        //(n - idx < segmentsPerBlock) 
                                        //? headIndexStart + (n - idx)
                                        //: headIndexStop;
        for (IndexType headIdx = headIndexStart; headIdx < headLoopStop; ++headIdx) {
            vec.template internalWrite<NumberSegments>(headIdx, initData);
        }
    }
}

template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static OriginalPointerType norm1(IndexType n, const SplitPtr in)
{
    #if !defined(USE_SINGLE_FOR_BLOCK)
    const IndexType undoneIterations = n % segmentsPerBlock;
    const IndexType endParallelLoop = n - undoneIterations;
    #else
    const IndexType endParallelLoop = n;
    #endif
    
    using std::abs;
    OriginalPointerType sum{0};
    #ifdef OMP_ENABLE_
    #pragma omp parallel for reduction(+:sum) schedule(static, blockChunk)
    #endif
    for (std::size_t i = 0; i < endParallelLoop; i += segmentsPerBlock) {
        const IndexType headIndexStart = SplitPtr::getHeadIndex(i);
        const IndexType headIndexStop = headIndexStart + segmentsPerBlock;
        #if !defined(USE_SINGLE_FOR_BLOCK)
        const IndexType headLoopStop = headIndexStop;
        #else
        const IndexType headLoopStop = (n - i < segmentsPerBlock) 
                                        ? headIndexStart + (n - i)
                                        : headIndexStop;
        #endif
        for (IndexType headIdx = headIndexStart; headIdx < headLoopStop; ++headIdx) {
            sum += abs(in.template internalRead<SegmentsToUse>(headIdx));
        }
    }
    
    #if !defined(USE_SINGLE_FOR_BLOCK)
    OriginalPointerType localSum{0};
    const IndexType headIndexStart = SplitPtr::getHeadIndex(endParallelLoop);
    const IndexType headLoopStop = headIndexStart + undoneIterations;
    for (IndexType headIdx = headIndexStart; headIdx < headLoopStop; ++headIdx) {
        localSum += abs(in.template internalRead<SegmentsToUse>(headIdx));
    }
    
    sum += localSum;
    #endif
    
    return sum;
}


template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static OriginalPointerType norm1Diff(IndexType n, const SplitPtr in1,
    const SplitPtr in2)
{
    #if !defined(USE_SINGLE_FOR_BLOCK)
    const IndexType undoneIterations = n % segmentsPerBlock;
    const IndexType endParallelLoop = n - undoneIterations;
    #else
    const IndexType endParallelLoop = n;
    #endif
    
    using std::abs;
    OriginalPointerType diff{0};
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for reduction(+:diff) schedule(static, blockChunk)
    #endif
    for (IndexType i = 0; i < endParallelLoop; i += segmentsPerBlock) {
        const IndexType headIndexStart = SplitPtr::getHeadIndex(i);
        const IndexType headIndexStop = headIndexStart + segmentsPerBlock;
        #if !defined(USE_SINGLE_FOR_BLOCK)
        const IndexType headLoopStop = headIndexStop;
        #else
        const IndexType headLoopStop = (n - i < segmentsPerBlock) 
                                        ? headIndexStart + (n - i)
                                        : headIndexStop;
        #endif
        for (IndexType headIdx = headIndexStart; headIdx < headLoopStop; ++headIdx) {
            diff += abs(in1.template internalRead<SegmentsToUse>(headIdx) 
                    - in2.template internalRead<SegmentsToUse>(headIdx));
        }
    }
    
    #if !defined(USE_SINGLE_FOR_BLOCK)
    OriginalPointerType localDiff{0};
    const IndexType headIndexStart = SplitPtr::getHeadIndex(endParallelLoop);
    const IndexType headLoopStop = headIndexStart + undoneIterations;
    for (IndexType headIdx = headIndexStart; headIdx < headLoopStop; ++headIdx) {
        localDiff += abs(in1.template internalRead<SegmentsToUse>(headIdx) 
                - in2.template internalRead<SegmentsToUse>(headIdx));
    }
    diff += localDiff;
    #endif
    
    return diff;
}

template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static OriginalPointerType selectiveNorm1(
    IndexType indexVectorSize,
    const IndexType * __restrict__ vecIndices,
    const SplitPtr vecSplit)
{
    using std::abs;
    OriginalPointerType sum = 0.0;
    #ifdef OMP_ENABLE_
    #pragma omp parallel for reduction(+:sum) schedule(static, ompHelper::blas1ChunkSize)
    #endif
    for (std::size_t i = 0; i < indexVectorSize; ++i) {
        sum += abs(vecSplit.template read<SegmentsToUse>(vecIndices[i]));
    }
    return sum;
}


template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, typename IndexType>
static void ellPageRank(IndexType n, IndexType k, const IndexType * __restrict__ indices,
    const SplitPtr valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPtr oldRankSplit,
   SplitPtr newRankSplit)
{
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(guided, ompHelper::chunkSize)
    #endif
    for (std::size_t row = 0; row < n; ++row) {
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        for (IndexType i = 0; i < k; ++i) {
            const IndexType matIdx = i*n + row;
            const IndexType rankIdx = indices[matIdx];
            const OriginalPointerType matValue = valsSplit.template read<partPointersToRead>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.template read<partPointersToRead>(rankIdx);
            localRes += matValue * rankValue;
        }

        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * (localRes + correction);
        newRankSplit.template write<partPointersToWrite>(row, newRank);
    }
}

template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, typename IndexType>
static void csrPageRank(IndexType m, IndexType n, const IndexType * __restrict__ colIndices,
    const IndexType *__restrict__ rowStartIndices,
    const SplitPtr valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPtr oldRankSplit,
   SplitPtr newRankSplit)
{
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(guided, ompHelper::chunkSize)
    #endif
    for (std::size_t row = 0; row < m; ++row) {
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        for (IndexType i = rowStartIndices[row]; i < rowStartIndices[row + 1]; i+= 1) {
            const IndexType matIdx = i;
            const IndexType rankIdx = colIndices[i];
            const OriginalPointerType matValue = valsSplit.template read<partPointersToRead>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.template read<partPointersToRead>(rankIdx);
            localRes += matValue * rankValue;
        }
        
        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * (localRes + correction);
        newRankSplit.template write<partPointersToWrite>(row, newRank);
    }
}


template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static void normalizeVector(IndexType n,
       SplitPtr vecValSplit,
       OriginalPointerType norm)
{
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(static, blockChunk)
    #endif
    for (std::size_t idx = 0; idx < n; idx += segmentsPerBlock) {
        const IndexType headIndexStart = SplitPtr::getHeadIndex(idx);
        const IndexType headIndexStop = headIndexStart + segmentsPerBlock;
        const IndexType headLoopStop = headIndexStop;
        for (IndexType headIdx = headIndexStart; headIdx < headLoopStop; ++headIdx) {
            const OriginalPointerType vecValue 
                    = vecValSplit.template internalRead<SegmentsToUse>(headIdx);
            vecValSplit.template internalWrite<SegmentsToUse>(headIdx, vecValue / norm);
        }
    }
}


/* Thread assignment should be in a way so one thread block handles an integer amount of memory blocks
   After execution, the segment pointer of vec should only be used as a pointer of the type OriginalPointerType.
*/
template<typename IndexType>
static void convertToOriginal(IndexType n, SplitPtr vec)
{
    if (NumberSegments <= 1) {
        return;
    }
    
    ThreadBuffer::getInstance().reserveBufferSizePerThread(NumberSegments * memoryBlockSize);
    
    constexpr std::size_t segmentsPerBlock = SplitPtr::segmentsPerMemoryBlock;
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(static, blockChunk)
    #endif
    for (std::size_t idx = 0; idx < n; idx += segmentsPerBlock) {
        
        //OriginalPointerType tempValues[segmentsPerBlock];
        OriginalPointerType *__restrict__ tempValues = reinterpret_cast<OriginalPointerType *>(ThreadBuffer::getInstance().getBuffer());
        OriginalPointerType *convertedPtr = reinterpret_cast<OriginalPointerType *>(vec.segmentPtr);
        
        const IndexType headIndexStart = SplitPtr::getHeadIndex(idx);
        
        // read values temporarily
        for (std::size_t i = 0; i < segmentsPerBlock; ++i) {
            tempValues[i] = vec.template internalRead<NumberSegments>(headIndexStart + i);
        }
        
        // write values
        for (std::size_t i = 0; i < segmentsPerBlock; ++i) {
            convertedPtr[idx + i] = tempValues[i];
        }
    }
}


// Thread assignment should be in a way so one thread block handles an integer amount of memory blocks
template<typename IndexType>
static void convertToSplitPointer(IndexType n, SplitPtr vec)
{
    if (NumberSegments <= 1) {
        return;
    }
    
    ThreadBuffer::getInstance().reserveBufferSizePerThread(NumberSegments * memoryBlockSize);
    
    constexpr std::size_t segmentsPerBlock = SplitPtr::segmentsPerMemoryBlock;
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(static, blockChunk)
    #endif
    for (std::size_t idx = 0; idx < n; idx += segmentsPerBlock) {
        
        //OriginalPointerType tempValues[segmentsPerBlock];
        OriginalPointerType *__restrict__ tempValues = reinterpret_cast<OriginalPointerType *>(ThreadBuffer::getInstance().getBuffer());
        OriginalPointerType *convertedPtr = reinterpret_cast<OriginalPointerType *>(vec.segmentPtr);
        
        // read values temporarily
        for (std::size_t i = 0; i < segmentsPerBlock; ++i) {
            tempValues[i] = convertedPtr[idx + i];
        }
        
        const IndexType headIndexStart = SplitPtr::getHeadIndex(idx);
        
        // write values
        for (std::size_t i = 0; i < segmentsPerBlock; ++i) {
            vec.template internalWrite<NumberSegments>(headIndexStart + i, tempValues[i]);
        }
    }
}

};


template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
struct EllPageRankStruct {

};

template<>
struct EllPageRankStruct<1, double> {

static constexpr MySegmentsNumberType NumberSegments = 1;
using OriginalPointerType = double;
using SplitPtr = SplitPointer<NumberSegments, OriginalPointerType>;

union ConvertAvx {__m256d avx; double dbl[4]; };

//vectorize over row
template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, typename IndexType>
static void ellPageRank(IndexType n, IndexType k, const IndexType * __restrict__ indices,
    const SplitPtr valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPtr oldRankSplit,
   SplitPtr newRankSplit)
{
#ifndef USE_OWN_ELL_SPMV_VECTORIZATION
    BasicKernels<NumberSegments, OriginalPointerType>::template ellPageRank<partPointersToRead, partPointersToWrite>(
                n, k, indices, valsSplit, dampingFactor, correction, oldRankSplit,
                newRankSplit);
#else
    static_assert(partPointersToRead <= NumberSegments, "It is not possible to read more segments than are available!");
    static_assert(partPointersToWrite <= NumberSegments, "It is not possible to write more segments than are available!");
    static_assert(std::is_same<IndexType, unsigned int>::value 
                    || std::is_same<IndexType, int>::value
                    || std::is_same<IndexType, std::int32_t>::value
                    || std::is_same<IndexType, std::uint32_t>::value, 
                    "IndexType must be a 32 bit integer (or unsigned integet)");
    
    const IndexType avxUpperLimit = n - (4-1);
    IndexType row = 0;
    const IndexType rowAfterAvx = n - (n % 4);
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for lastprivate(row) schedule(guided, ompHelper::chunkSize / 4)
    #endif
    for (row = 0; row < avxUpperLimit; row += 4) {
        __m256d resPack = _mm256_setzero_pd();
        
        for (IndexType i = 0; i < k; i += 1) {
            const IndexType matStartIdx = i*n + row;
            
            __m128i rankIndicesPack = _mm_loadu_si128((__m128i *)(indices + matStartIdx));
            
            __m256d matValPack = 
                avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (valsSplit, matStartIdx);
            
            __m256d rankPack = 
                avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (oldRankSplit, rankIndicesPack);
            
            __m256d mulPack = _mm256_mul_pd(matValPack, rankPack);
            resPack = _mm256_add_pd(resPack, mulPack);
        }
        
        ConvertAvx result;
        result.avx = resPack;
        
        #pragma unroll
        for (std::size_t i = 0; i < 4; ++i) {
            const OriginalPointerType newRank =
                (static_cast<OriginalPointerType>(1) - dampingFactor)
                 / static_cast<OriginalPointerType>(n)
                + dampingFactor * (result.dbl[i] + correction);
            
            newRankSplit.template write<partPointersToWrite>(row + i, newRank);
        }
    }
    
    //if (row != rowAfterAvx) {
    //    std::cerr << "Mismatch! n = " << n << "; Row = " << row << ", expected: " << rowAfterAvx << "\n";
    //}
    
    // After the loop above, row should be:
    // row = rowAfterAvx
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    //row must be the same as before!
    for (row = rowAfterAvx; row < n; ++row) {
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        for (IndexType i = 0; i < k; ++i) {
            const IndexType matIdx = i*n + row;
            const IndexType rankIdx = indices[matIdx];
            const OriginalPointerType matValue = valsSplit.template read<partPointersToRead>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.template read<partPointersToRead>(rankIdx);
            localRes += matValue * rankValue;
        }

        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * (localRes + correction);
        newRankSplit.template write<partPointersToWrite>(row, newRank);
    }
#endif
}

};


template<>
struct EllPageRankStruct<2, double> {

static constexpr MySegmentsNumberType NumberSegments = 2;
using OriginalPointerType = double;
using SplitPtr = SplitPointer<NumberSegments, OriginalPointerType>;

union ConvertAvx {__m256d avx; double dbl[4]; };

//vectorize over row
template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, typename IndexType>
static void ellPageRank(IndexType n, IndexType k, const IndexType * __restrict__ indices,
    const SplitPtr valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPtr oldRankSplit,
   SplitPtr newRankSplit)
{
#ifndef USE_OWN_ELL_SPMV_VECTORIZATION
    BasicKernels<NumberSegments, OriginalPointerType>::template ellPageRank<partPointersToRead, partPointersToWrite>(
                n, k, indices, valsSplit, dampingFactor, correction, oldRankSplit,
                newRankSplit);
#else
    static_assert(partPointersToRead <= NumberSegments, "It is not possible to read more segments than are available!");
    static_assert(partPointersToWrite <= NumberSegments, "It is not possible to write more segments than are available!");
    static_assert(std::is_same<IndexType, unsigned int>::value 
                    || std::is_same<IndexType, int>::value
                    || std::is_same<IndexType, std::int32_t>::value
                    || std::is_same<IndexType, std::uint32_t>::value, 
                    "IndexType must be a 32 bit integer (or unsigned integet)");
    
    const IndexType avxUpperLimit = n - (8-1);
    IndexType row = 0;
    const IndexType rowAfterAvx = n - (n % 8);
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for lastprivate(row) schedule(guided, ompHelper::chunkSize / 8)
    #endif
    for (row = 0; row < avxUpperLimit; row += 8) {
        __m256d resPack[2] = { _mm256_setzero_pd(), _mm256_setzero_pd() };
        
        for (IndexType i = 0; i < k; i += 1) {
            const IndexType matStartIdx[2] = { i*n + row, i*n + row + 4 };
            
            __m128i rankIndicesPack[2] = 
                { _mm_loadu_si128((__m128i *)(indices + matStartIdx[0])), 
                  _mm_loadu_si128((__m128i *)(indices + matStartIdx[1]))  };
            
            __m256d matValPack[2] = 
                { avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (valsSplit, matStartIdx[0]),
                  avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (valsSplit, matStartIdx[1])  };
            
            __m256d rankPack[2] = 
                { avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (oldRankSplit, rankIndicesPack[0]), 
                  avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (oldRankSplit, rankIndicesPack[1])  };
            
            __m256d mulPack[2] = 
                { _mm256_mul_pd(matValPack[0], rankPack[0]),
                  _mm256_mul_pd(matValPack[1], rankPack[1]) };
            resPack[0] = _mm256_add_pd(resPack[0], mulPack[0]);
            resPack[1] = _mm256_add_pd(resPack[1], mulPack[1]);
        }
        
        ConvertAvx result[2];
        result[0].avx = resPack[0];
        result[1].avx = resPack[1];
        
        #pragma unroll
        for (std::size_t i = 0; i < 4; ++i) {
            const OriginalPointerType newRank[2] = 
                {
                (static_cast<OriginalPointerType>(1) - dampingFactor)
                 / static_cast<OriginalPointerType>(n)
                + dampingFactor * (result[0].dbl[i] + correction), 
                
                (static_cast<OriginalPointerType>(1) - dampingFactor)
                 / static_cast<OriginalPointerType>(n)
                + dampingFactor * (result[1].dbl[i] + correction)
                };
            
            newRankSplit.template write<partPointersToWrite>(row + i, newRank[0]);
            newRankSplit.template write<partPointersToWrite>(row + i + 4, newRank[1]);
        }
    }
    
    //if (row != rowAfterAvx) {
    //    std::cerr << "Mismatch! n = " << n << "; Row = " << row << ", expected: " << rowAfterAvx << "\n";
    //}
    
    // After the loop above, row should be:
    // row = rowAfterAvx
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    //row must be the same as before!
    for (row = rowAfterAvx; row < n; ++row) {
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        for (IndexType i = 0; i < k; ++i) {
            const IndexType matIdx = i*n + row;
            const IndexType rankIdx = indices[matIdx];
            const OriginalPointerType matValue = valsSplit.template read<partPointersToRead>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.template read<partPointersToRead>(rankIdx);
            localRes += matValue * rankValue;
        }

        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * (localRes + correction);
        newRankSplit.template write<partPointersToWrite>(row, newRank);
    }
#endif
}

};




template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
struct Avx2Kernels
{

union ConvertAvx {__m256d avx; double dbl[4]; };

using SplitPtr = SplitPointer<NumberSegments, OriginalPointerType>;

static constexpr std::size_t segmentsPerBlock = SplitPtr::segmentsPerMemoryBlock;

static constexpr std::size_t potentialChunk = ompHelper::blas1ChunkSize / segmentsPerBlock;

static constexpr int blockChunk = (1 < potentialChunk) ? potentialChunk : 1;

#ifdef OMP_ENABLE_
#pragma omp declare reduction \
  (vecAdd:__m256d:omp_out=_mm256_add_pd(omp_out, omp_in)) \
  initializer(omp_priv=_mm256_setzero_pd())
#endif

template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static void rwTest(IndexType n, SplitPtr rw, OriginalPointerType factor)
{
    static_assert(SegmentsToUse <= NumberSegments, "It is not possible to use more segments than are available!");
    const std::size_t avxUpperLimit = n - (4-1);
    IndexType idx = 0;
    #ifdef OMP_ENABLE_
    #pragma omp parallel for lastprivate(idx) schedule(static, ompHelper::blas1ChunkSize / 4)
    #endif
    for (idx = 0; idx < avxUpperLimit; idx += 4) {
        const __m256d factorPack = _mm256_set1_pd(factor);
        
        __m256d valPack = avx2Helper::readPack<SegmentsToUse, NumberSegments, OriginalPointerType>(rw, idx);
        __m256d multResPack = _mm256_mul_pd(valPack, factorPack);
        avx2Helper::writePack<SegmentsToUse, NumberSegments, OriginalPointerType>(rw, multResPack, idx);
    }
    //idx = n - (4 - (n % 4));
    for (; idx < n; idx += 1) {
        OriginalPointerType val = rw.template read<SegmentsToUse>(idx);
        rw.template write<SegmentsToUse>(idx, val * factor);
    }
}

//No need to vectorize because not runtime critical!
template<typename IndexType>
static void setVector(IndexType n, SplitPtr vec, OriginalPointerType initData)
{
    BasicKernels<NumberSegments, OriginalPointerType>::setVector(n, vec, initData);
}

template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static OriginalPointerType norm1(IndexType n, const SplitPtr in)
{
    static_assert(SegmentsToUse <= NumberSegments, "It is not possible to use more segments than are available!");
    using std::abs;
    const std::size_t avxUpperLimit = n - (4-1);
    
    std::size_t i = 0;
    __m256d sumPack = _mm256_setzero_pd();
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for reduction(vecAdd:sumPack) lastprivate(i) schedule(static, ompHelper::blas1ChunkSize / 4)
    #endif
    for (i = 0; i < avxUpperLimit; i += 4) {
        //__m128i indices32 = _mm_setr_epi32(i, i+1, i+2, i+3);
        //__m256d valPack = avx2Helper::readPack<SegmentsToUse, NumberSegments, OriginalPointerType>(in, indices32);
        __m256d valPack = avx2Helper::readPack<SegmentsToUse, NumberSegments, OriginalPointerType>(in, i);
        __m256d absValPack = avx2Helper::absPack(valPack);
        sumPack = _mm256_add_pd(sumPack, absValPack);
    }
    ConvertAvx result;
    result.avx = sumPack;
    
    OriginalPointerType sum = result.dbl[0] + result.dbl[1] + result.dbl[2] + result.dbl[3];
    
    for (; i < n; i += 1) {
        sum += abs(in.template read<SegmentsToUse>(i));
    }
    return sum;
}


template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static OriginalPointerType norm1Diff(IndexType n, const SplitPtr in1,
    const SplitPtr in2)
{
    static_assert(SegmentsToUse <= NumberSegments, "It is not possible to use more segments than are available!");
    using std::abs;
    const std::size_t avxUpperLimit = n - (4-1);
    
    std::size_t i = 0;
    __m256d diffPack = _mm256_setzero_pd();
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for reduction(vecAdd:diffPack) lastprivate(i) schedule(static, ompHelper::blas1ChunkSize / 4)
    #endif
    for (i = 0; i < avxUpperLimit; i += 4) {
        __m256d valPack1 = avx2Helper::readPack<SegmentsToUse, NumberSegments, OriginalPointerType>(in1, i);
        __m256d valPack2 = avx2Helper::readPack<SegmentsToUse, NumberSegments, OriginalPointerType>(in2, i);
        __m256d subResPack = _mm256_sub_pd(valPack1, valPack2);
        __m256d absValPack = avx2Helper::absPack(subResPack);
        diffPack = _mm256_add_pd(diffPack, absValPack);
    }
    ConvertAvx result;
    result.avx = diffPack;
    
    OriginalPointerType diff = result.dbl[0] + result.dbl[1] + result.dbl[2] + result.dbl[3];
    for (; i < n; ++i) {
        diff += abs(in1.template read<SegmentsToUse>(i) - in2.template read<SegmentsToUse>(i));
    }

    return diff;
}

template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static OriginalPointerType selectiveNorm1(
    IndexType indexVectorSize,
    const IndexType * __restrict__ vecIndices,
    const SplitPtr vecSplit)
{
    return BasicKernels<NumberSegments, OriginalPointerType>::template selectiveNorm1<SegmentsToUse>(
                indexVectorSize, vecIndices, vecSplit);
}



//vectorize over row
template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, typename IndexType>
static void ellPageRank(IndexType n, IndexType k, const IndexType * __restrict__ indices,
    const SplitPtr valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPtr oldRankSplit,
   SplitPtr newRankSplit)
{
#ifndef USE_OWN_ELL_SPMV_VECTORIZATION
    BasicKernels<NumberSegments, OriginalPointerType>::template ellPageRank<partPointersToRead, partPointersToWrite>(
                n, k, indices, valsSplit, dampingFactor, correction, oldRankSplit,
                newRankSplit);
#else
/*
    EllPageRankStruct<NumberSegments, OriginalPointerType>::template ellPageRank<partPointersToRead, partPointersToWrite>(
                n, k, indices, valsSplit, dampingFactor, correction, oldRankSplit, newRankSplit);
/*/
    static_assert(partPointersToRead <= NumberSegments, "It is not possible to read more segments than are available!");
    static_assert(partPointersToWrite <= NumberSegments, "It is not possible to write more segments than are available!");
    static_assert(std::is_same<IndexType, unsigned int>::value 
                    || std::is_same<IndexType, int>::value
                    || std::is_same<IndexType, std::int32_t>::value
                    || std::is_same<IndexType, std::uint32_t>::value, 
                    "IndexType must be a 32 bit integer (or unsigned integet)");
    
    const IndexType avxUpperLimit = n - (4-1);
    IndexType row = 0;
    const IndexType rowAfterAvx = n - (n % 4);
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for lastprivate(row) schedule(guided, ompHelper::chunkSize / 4)
    #endif
    for (row = 0; row < avxUpperLimit; row += 4) {
        __m256d resPack = _mm256_setzero_pd();
        
        for (IndexType i = 0; i < k; i += 1) {
            const IndexType matStartIdx = i*n + row;
            
            __m128i rankIndicesPack = _mm_loadu_si128((__m128i *)(indices + matStartIdx));
            
            __m256d matValPack = 
                avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (valsSplit, matStartIdx);
            
            __m256d rankPack = 
                avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (oldRankSplit, rankIndicesPack);
            
            __m256d mulPack = _mm256_mul_pd(matValPack, rankPack);
            resPack = _mm256_add_pd(resPack, mulPack);
        }
        
        ConvertAvx result;
        result.avx = resPack;
        
        #pragma unroll
        for (std::size_t i = 0; i < 4; ++i) {
            const OriginalPointerType newRank =
                (static_cast<OriginalPointerType>(1) - dampingFactor)
                 / static_cast<OriginalPointerType>(n)
                + dampingFactor * (result.dbl[i] + correction);
            
            newRankSplit.template write<partPointersToWrite>(row + i, newRank);
        }
    }
    
    //if (row != rowAfterAvx) {
    //    std::cerr << "Mismatch! n = " << n << "; Row = " << row << ", expected: " << rowAfterAvx << "\n";
    //}
    
    // After the loop above, row should be:
    // row = rowAfterAvx
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    //row must be the same as before!
    for (row = rowAfterAvx; row < n; ++row) {
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        for (IndexType i = 0; i < k; ++i) {
            const IndexType matIdx = i*n + row;
            const IndexType rankIdx = indices[matIdx];
            const OriginalPointerType matValue = valsSplit.template read<partPointersToRead>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.template read<partPointersToRead>(rankIdx);
            localRes += matValue * rankValue;
        }

        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * (localRes + correction);
        newRankSplit.template write<partPointersToWrite>(row, newRank);
    }
//*/
#endif
}
/*/
// vectorize over k
template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, typename IndexType>
static void ellPageRank(IndexType n, IndexType k, const IndexType * __restrict__ indices,
    const SplitPtr valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPtr oldRankSplit,
   SplitPtr newRankSplit)
{
#ifndef USE_OWN_ELL_SPMV_VECTORIZATION
    BasicKernels<NumberSegments, OriginalPointerType>::template ellPageRank<partPointersToRead, partPointersToWrite>(
                n, k, indices, valsSplit, dampingFactor, correction, oldRankSplit,
                newRankSplit);
#else
    
    static_assert(partPointersToRead <= NumberSegments, "It is not possible to read more segments than are available!");
    static_assert(partPointersToWrite <= NumberSegments, "It is not possible to write more segments than are available!");
    static_assert(std::is_same<IndexType, unsigned int>::value 
                    || std::is_same<IndexType, int>::value
                    || std::is_same<IndexType, std::int32_t>::value
                    || std::is_same<IndexType, std::uint32_t>::value, 
                    "IndexType must be a 32 bit integer (or unsigned integet)");
    
    const IndexType avxUpperLimit = k - (4-1);
    IndexType row = 0;
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for lastprivate(row) schedule(guided, ompHelper::chunkSize)
    #endif
    for (IndexType row = 0; row < n; ++row) {
        __m256d resPack = _mm256_setzero_pd();
        IndexType i = 0;
        
        for (i = 0; i < avxUpperLimit; i += 4) {
            __m128i matIndicesPack = _mm_setr_epi32(i, i+1, i+2, i+3);
            matIndicesPack = _mm_mullo_epi32(matIndicesPack, _mm_set1_epi32(n));
            matIndicesPack = _mm_add_epi32(matIndicesPack, _mm_set1_epi32(row));
            
            __m128i rankIndicesPack = _mm_i32gather_epi32(reinterpret_cast<const int *>(indices), matIndicesPack, 4);
            
            __m256d matValPack = 
                avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (valsSplit, matIndicesPack);
            
            __m256d rankPack = 
                avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                    (oldRankSplit, rankIndicesPack);
            
            __m256d mulPack = _mm256_mul_pd(matValPack, rankPack);
            resPack = _mm256_add_pd(resPack, mulPack);
        }
        
        OriginalPointerType localRes = static_cast<OriginalPointerType>(0.0);
        
        for (; i < k; ++i) {
            const IndexType matStartIdx = i*n + row;
            
            const IndexType rankIdx = indices[matStartIdx];
            
            const OriginalPointerType matValue = valsSplit.template read<partPointersToRead>(matStartIdx);
            const OriginalPointerType rankValue = oldRankSplit.template read<partPointersToRead>(rankIdx);
            localRes += matValue * rankValue;
        }
        
        //REDUCTION
        ConvertAvx result;
        result.avx = resPack;
        localRes = localRes + result.dbl[0] + result.dbl[1]
                    + result.dbl[2] + result.dbl[3];
                
        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * (localRes + correction);
        newRankSplit.template write<partPointersToWrite>(row, newRank);
        
    }
#endif
}
//*/

template<MySegmentsNumberType partPointersToRead, MySegmentsNumberType partPointersToWrite, typename IndexType>
static void csrPageRank(IndexType m, IndexType n, const IndexType * __restrict__ colIndices,
    const IndexType *__restrict__ rowStartIndices,
    const SplitPtr valsSplit,
    OriginalPointerType dampingFactor,
    OriginalPointerType correction,
    const SplitPtr oldRankSplit,
   SplitPtr newRankSplit)
{
#ifndef USE_OWN_CSR_SPMV_VECTORIZATION
    BasicKernels<NumberSegments, OriginalPointerType>::template csrPageRank<partPointersToRead, partPointersToWrite>(
                m, n, colIndices, rowStartIndices, valsSplit, dampingFactor, correction, 
                oldRankSplit, newRankSplit);
#else
    static_assert(partPointersToRead <= NumberSegments, "It is not possible to read more segments than are available!");
    static_assert(partPointersToWrite <= NumberSegments, "It is not possible to write more segments than are available!");
    static_assert(std::is_same<IndexType, unsigned int>::value 
                    || std::is_same<IndexType, int>::value
                    || std::is_same<IndexType, std::int32_t>::value
                    || std::is_same<IndexType, std::uint32_t>::value, 
                    "IndexType must be a 32 bit integer (or unsigned integet)");
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(guided, ompHelper::chunkSize)
    #endif
    for (std::size_t row = 0; row < m; ++row) {
        IndexType i = rowStartIndices[row];
        __m256d currentResPack = _mm256_setr_pd(0.0, 0.0, 0.0, 0.0);
        if (rowStartIndices[row + 1] > (4-1)) {
            for (i = rowStartIndices[row]; i < rowStartIndices[row + 1] - (4-1); i += 4) {
                const __m128i rankIndices32 = _mm_loadu_si128((__m128i *) (colIndices + i));
                __m256d matValPack =
                    avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                        (valsSplit, i);
                __m256d rankPack =
                    avx2Helper::readPack<partPointersToRead, NumberSegments, OriginalPointerType>
                        (oldRankSplit, rankIndices32);
                __m256d multResPack = _mm256_mul_pd(matValPack, rankPack);
                currentResPack = _mm256_add_pd(currentResPack, multResPack);
            }
        }
        
        double currentResult = 0.0;
        for (; i < rowStartIndices[row + 1]; ++i) {
            const IndexType matIdx = i;
            const IndexType rankIdx = colIndices[i];
            const OriginalPointerType matValue = valsSplit.template read<partPointersToRead>(matIdx);
            const OriginalPointerType rankValue = oldRankSplit.template read<partPointersToRead>(rankIdx);
            currentResult += matValue * rankValue;
        }

        //REDUCTION
        ConvertAvx result;
        result.avx = currentResPack;
        currentResult = currentResult + result.dbl[0]
                + result.dbl[1] + result.dbl[2] + result.dbl[3];

        const OriginalPointerType newRank =
            (static_cast<OriginalPointerType>(1) - dampingFactor)
             / static_cast<OriginalPointerType>(n)
            + dampingFactor * (currentResult + correction);
        newRankSplit.template write<partPointersToWrite>(row, newRank);
    }
#endif
}


template<MySegmentsNumberType SegmentsToUse, typename IndexType>
static void normalizeVector(IndexType n,
       SplitPtr vecValSplit,
       OriginalPointerType norm)
{
    BasicKernels<NumberSegments, OriginalPointerType>::template 
            normalizeVector<SegmentsToUse>(n, vecValSplit, norm);
}


//template<typename IndexType>
//static void convertToOriginal(IndexType n, SplitPtr vec)
//{
//    BasicKernels<NumberSegments, OriginalPointerType>::convertToOriginal(n, vec);
//}

/* Thread assignment should be in a way so one thread block handles an integer amount of memory blocks
   After execution, the segment pointer of vec should only be used as a pointer of the type OriginalPointerType.
*/

template<typename IndexType>
static void convertToOriginal(IndexType n, SplitPtr vec)
{
    if (NumberSegments <= 1) {
        return;
    }
    
    ThreadBuffer::getInstance().reserveBufferSizePerThread(NumberSegments * memoryBlockSize);
    
    constexpr std::size_t segmentsPerBlock = SplitPtr::segmentsPerMemoryBlock;
    
    #ifdef OMP_ENABLE_
    #pragma omp parallel for schedule(static, blockChunk)
    #endif
    for (IndexType idx = 0; idx < n; idx += segmentsPerBlock) {
        
        //OriginalPointerType tempValues[segmentsPerBlock];
        __m256d *__restrict__ tempValues = reinterpret_cast<__m256d *>(ThreadBuffer::getInstance().getBuffer());
        OriginalPointerType *convertedPtr = reinterpret_cast<OriginalPointerType *>(vec.segmentPtr);
        
        const IndexType headIndexStart = SplitPtr::getHeadIndex(idx);
        
        // read values temporarily
        for (IndexType i = 0, vecIdx = 0; i < segmentsPerBlock; i += 4, ++vecIdx) {
            tempValues[vecIdx] = avx2Helper::directReadPack<NumberSegments, NumberSegments, OriginalPointerType>(vec, headIndexStart + i);
        }
        
        // write values
        OriginalPointerType *origTempValues = reinterpret_cast<OriginalPointerType *>(tempValues);
        for (IndexType i = 0; i < segmentsPerBlock; ++i) {
            convertedPtr[idx + i] = origTempValues[i];
        }
    }
}


// Thread assignment should be in a way so one thread block handles an integer amount of memory blocks
template<typename IndexType>
static void convertToSplitPointer(IndexType n, SplitPtr vec)
{
    BasicKernels<NumberSegments, OriginalPointerType>::convertToSplitPointer(n, vec);
}



};


template<MySegmentsNumberType NumberSegments, typename OriginalPointerType>
struct KernelSelector : BasicKernels<NumberSegments, OriginalPointerType> {};

//Make the actual specializations here, so AVX is used for the ones below.
#if defined(USE_VECTORIZATION) && defined(__AVX2__)
template<> struct KernelSelector<1, double> : public Avx2Kernels<1, double> {};
template<> struct KernelSelector<2, double> : public Avx2Kernels<2, double> {};
template<> struct KernelSelector<4, double> : public Avx2Kernels<4, double> {};

#endif


}   // namespace detail

}   // namespace kernels


#if defined(USE_VECTORIZATION)
#define USE_VECTORIZATION_STRING "vectorized "
#else
#define USE_VECTORIZATION_STRING ""
#endif

#if defined(USE_OWN_CSR_SPMV_VECTORIZATION) && defined(USE_VECTORIZATION)
#define USE_OWN_CSR_SPMV_VECTORIZATION_STRING "own_csr "
#else
#define USE_OWN_CSR_SPMV_VECTORIZATION_STRING ""
#endif

#if defined(USE_OWN_ELL_SPMV_VECTORIZATION) && defined(USE_VECTORIZATION)
#define USE_OWN_ELL_SPMV_VECTORIZATION_STRING "own_ell "
#else
#define USE_OWN_ELL_SPMV_VECTORIZATION_STRING ""
#endif

#if defined(USE_SINGLE_FOR_BLOCK)
#define USE_SINGLE_FOR_BLOCK_STRING "single_block "
#else
#define USE_SINGLE_FOR_BLOCK_STRING ""
#endif

#define KERNEL_ENABLED_MACROS USE_VECTORIZATION_STRING USE_OWN_CSR_SPMV_VECTORIZATION_STRING USE_OWN_ELL_SPMV_VECTORIZATION_STRING USE_SINGLE_FOR_BLOCK_STRING


#undef USE_OWN_CSR_SPMV_VECTORIZATION
#undef USE_OWN_ELL_SPMV_VECTORIZATION
#undef USE_VECTORIZATION
#undef USE_SINGLE_FOR_BLOCK



#endif //CPU_KERNEL_SELECTOR_HPP
