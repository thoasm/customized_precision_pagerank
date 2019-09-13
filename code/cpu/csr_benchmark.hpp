#ifndef CSR_BENCHMARK_HPP
#define CSR_BENCHMARK_HPP
#include "split_pointer.hpp"
#include "split_pointer_helper.hpp"
#include "coo_matrix.hpp"
#include "csr_matrix.hpp"
#include "cpu_kernels.hpp"
#include "constants.hpp"

#include <cstdint>
#include <vector>
#include <chrono>
#include <iostream>
//#include <cstdlib>
//#include <ctime>
#include <algorithm>
#include <utility>
#include <iomanip>


#define SPLIT_POINTER_TEMPLATE_HEAD MySegmentsNumberType NumberSegments, typename OriginalPointerType
#define SPLIT_POINTER_TEMPLATES NumberSegments, OriginalPointerType

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
class CsrBenchmark {
public:
    CsrBenchmark(std::size_t loopCount = 100) : loopCount_(loopCount), vecInitialized_(false), matInitialized_(false),
                  normResult_{0.0}, m_{0}, n_{0}, nnz_{0}, vecN_{0}, emptyColNumber_(0)  {
        srand(time(NULL));    
    }
    
    void setUpVectors(std::size_t n);
    void setVectors(OriginalPointerType initValues);
    void setUpVectors(const std::vector<OriginalPointerType> &n);
    void setUpMatrix(const CooMatrix<IndexType, OriginalPointerType> &mtx);
    void setUpMatrix(const CsrMatrix<IndexType, OriginalPointerType> &mtx);
    
    // returns ns for one execution, averaged over <loopCount_> runs
    template<int ReadPartsNumber>
    std::uint64_t runOwnNorm();
    
    template<int ReadPartsNumber>
    std::uint64_t runRWTest();
    
    template<int ReadPartsNumber>
    std::uint64_t runSpmv();
    
    std::uint64_t runConversion();

    std::uint64_t runNormalization(int readPrec);
    
    template<int ReadPartsNumber>
    std::uint64_t runDiffNorm();
    
    template<int ReadPartsNumber>
    std::uint64_t runSelectiveNorm();
    
    
    void testConversion();
    
    void clearVec() {
        if (vecInitialized_) {
            vecN_ = 0;
            deleteSplitPtr(x_);
            deleteSplitPtr(y_);
            
            vecInitialized_ = false;
        }
    }
    
    void clearMat() {
        if (matInitialized_) {
            m_ = 0;
            n_ = 0;
            nnz_ = 0;
            emptyColNumber_ = 0;
            
            deleteSplitPtr(vals_);
            colIndices_.clear();
            rowStartIndices_.clear();
            emptyColIndices_.clear();
            
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
    void newVectors(std::size_t vecSize);
    
    std::size_t loopCount_;
    static constexpr std::size_t maxRndNumber = 13;
    static constexpr int numberThreadsPerRow = 1;
    using SplitType = SplitPointer<SPLIT_POINTER_TEMPLATES>;
    
    
    std::vector<OriginalPointerType> vectorResult_;
    OriginalPointerType normResult_;
    
    bool vecInitialized_;
    bool matInitialized_;
    
    SplitType x_;
    SplitType y_;
    
    SplitType vals_;
    std::vector<IndexType> colIndices_;
    std::vector<IndexType> rowStartIndices_;
    std::vector<IndexType> emptyColIndices_;
    
    IndexType emptyColNumber_;
    IndexType m_, n_, nnz_, vecN_;
};

#include "csr_benchmark.inc.cpp"


#undef SPLIT_POINTER_TEMPLATE_HEAD

#undef SPLIT_POINTER_TEMPLATES

#endif  // CSR_BENCHMARK_CUH
