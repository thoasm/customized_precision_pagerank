#ifndef CSR_BENCHMARK_CUH
#define CSR_BENCHMARK_CUH
#include "split_pointer.hpp"
#include "split_pointer_helper.cuh"
#include "coo_matrix.hpp"
#include "csr_matrix.hpp"
#include "gpu_kernels.cuh"
#include "constants.hpp"

#include <cstdint>
#include <vector>
#include <chrono>
#include <iostream>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <utility>
#include <iomanip>


#define SPLIT_POINTER_TEMPLATE_HEAD MySegmentsNumberType NumberSegments, typename OriginalPointerType
#define SPLIT_POINTER_TEMPLATES NumberSegments, OriginalPointerType

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
class CsrBenchmark {
public:
    CsrBenchmark(std::size_t loopCount = 100) : loopCount_(loopCount), vecInitialized_(false), matInitialized_(false),
                  normResult_{0.0}, m_{0}, n_{0}, nnz_{0}, vecN_{0}, d_buffer_{nullptr}, d_colIndices_(nullptr), d_rowStartIndices_(nullptr), d_emptyColIndices_(nullptr), emptyColNumber_(0)  {
        srand(time(NULL));    
    }
    
    void setUpVectors(std::size_t n);
    void setVectors(OriginalPointerType initValues);
    void setUpVectors(const std::vector<OriginalPointerType> &n);
    void setUpMatrix(const CooMatrix<IndexType, OriginalPointerType> &mtx);
    void setUpMatrix(const CsrMatrix<IndexType, OriginalPointerType> &mtx);
    
    // returns ns for one execution, averaged over <loopCount_> runs
    template<int ReadPartsNumber, unsigned int BlockSize = cnst::blockSize>
    std::uint64_t runOwnNorm();
    
    template<int ReadPartsNumber, unsigned int BlockSize = cnst::blockSize>
    std::uint64_t runAltOwnNorm();
    
    //works only if calculated with normal doubles
    std::uint64_t runCuBlasNorm();
    
    template<int ReadPartsNumber, unsigned int BlockSize = cnst::blockSize>
    std::uint64_t runRWTest();
    
    template<int ReadPartsNumber, unsigned int BlockSize = cnst::blockSize>
    std::uint64_t runSpmv();
    
    template<unsigned int BlockSize = cnst::blockSize>
    std::uint64_t runConversion();

    template<unsigned int BlockSize = cnst::blockSize>
    std::uint64_t runNormalization(int readPrec);
    
    template<int ReadPartsNumber, unsigned int BlockSize = cnst::blockSize>
    std::uint64_t runDiffNorm();
    
    template<int ReadPartsNumber, unsigned int BlockSize = cnst::blockSize>
    std::uint64_t runSelectiveNorm();
    
    
    void testConversion();
    
    void clearVec() {
        if (vecInitialized_) {
            clearHostSplitVec();
            vecN_ = 0;
            freeSplitPtr(d_x_);
            freeSplitPtr(d_y_);
            freeSplitPtr(d_split_buffer_);
            cudaFree(d_buffer_);
            d_buffer_ = nullptr;
            vecInitialized_ = false;
        }
    }
    
    void clearMat() {
        if (matInitialized_) {
            m_ = 0;
            n_ = 0;
            nnz_ = 0;
            emptyColNumber_ = 0;
            
            freeSplitPtr(d_vals_);
            cudaFree(d_colIndices_);    d_colIndices_ = nullptr;
            cudaFree(d_rowStartIndices_);   d_rowStartIndices_ = nullptr;
            cudaFree(d_emptyColIndices_);   d_emptyColIndices_ = nullptr;
            
            matInitialized_ = false;
        }
    }
    void clear() {
        clearVec();
        clearMat();
    }
    
    ~CsrBenchmark() {
        clear();
    }
    
    const std::vector<OriginalPointerType> &getLatestVectorResult() const {
        return vectorResult_;
    }
    OriginalPointerType getLatestNormResult() const {   return normResult_; }

private:
    void mallocVectors(std::size_t vecSize);
    
    std::size_t loopCount_;
    static constexpr std::size_t maxRndNumber = 13;
    static constexpr int numberThreadsPerRow = 1;
    using SplitType = SplitPointer<SPLIT_POINTER_TEMPLATES>;
    
    SplitType hostSplitVec_;
    std::vector<typename SplitType::SegmentType> hostVecSegments_;
    
    void setUpHostSplitVec(std::size_t n) {
        hostVecSegments_.resize(SplitType::getMinimumBytesFor(n) / sizeof(typename SplitType::SegmentType));
        hostSplitVec_.segmentPtr = hostVecSegments_.data();
    }
    void clearHostSplitVec() {
        hostVecSegments_.clear();
        hostSplitVec_.segmentPtr = nullptr;
    }
    
    std::vector<OriginalPointerType> vectorResult_;
    OriginalPointerType normResult_;
    
    bool vecInitialized_;
    bool matInitialized_;
    OriginalPointerType *d_buffer_;
    SplitType d_x_;
    SplitType d_y_, d_split_buffer_;
    
    SplitType d_vals_;
    IndexType *d_colIndices_;
    IndexType *d_rowStartIndices_;
    IndexType *d_emptyColIndices_;
    
    IndexType emptyColNumber_;
    IndexType m_, n_, nnz_, vecN_;
};

#include "csr_benchmark.inc.cu"


#undef SPLIT_POINTER_TEMPLATE_HEAD

#undef SPLIT_POINTER_TEMPLATES

#endif  // CSR_BENCHMARK_CUH
