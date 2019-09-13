#ifndef CSR_PAGE_RANK_CUH
#define CSR_PAGE_RANK_CUH

#include "csr_matrix.hpp"
#include "split_pointer.hpp"
#include "split_pointer_helper.cuh"
#include "gpu_kernels.cuh"
#include "macros.hpp"
#include "precision.hpp"

#include <vector>
#include <cstdint>
#include <cmath>
#include <chrono>

#include <iostream>
#include <iomanip>

#define SPLIT_POINTER_TEMPLATE_HEAD MySegmentsNumberType NumberSegments, typename OriginalPointerType
#define SPLIT_POINTER_TEMPLATES NumberSegments, OriginalPointerType

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
class CsrPageRank {
public:
    CsrPageRank() : vectorInitialized_{false}, matrixInitialized_{false}, curXIdx_{0}, usedIterations_{0}, resultUpdated_{false}, precisionChanges_(NumberSegments+1, 0), precisionChanesTimePoint_(NumberSegments+1) {}
    
    void solve(OriginalPointerType dampingFactor, OriginalPointerType epsilon = 1e-10, std::size_t maxIterations = 1000);
    
    bool isInit() const {   return vectorInitialized_ && matrixInitialized_;    }
    
    bool initMatrix(const CsrMatrix<IndexType, OriginalPointerType> &csrMtx);
    void resetVector();
    
    const std::vector<OriginalPointerType> &getResult() const;
    std::size_t getUsedIterations() const   {   return usedIterations_; }
    
    const std::vector<std::size_t> &getPrecisionChanges() const { return precisionChanges_;   }
    const std::vector<std::chrono::high_resolution_clock::time_point> &getPrecisionTimePoints() const { return precisionChanesTimePoint_;   }
    
    void clearMtx();
    void clearVector();
    
    void clear() {
        clearMtx();
        clearVector();
    }
    
    ~CsrPageRank() {
        clear();
    }

    
private:
    using SplitType = SplitPointer<SPLIT_POINTER_TEMPLATES>;
    using RangeRW = macros::compile_int_list<1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16>;
    using RangeWrite = RangeRW;
    
    static constexpr std::uint32_t blockSize{cnst::blockSize};
    static constexpr std::int32_t threadsPerRow{cnst::numberThreadsPerRow};
    
    static constexpr std::size_t normCheckFrequency{cnst::normCheckFrequency};
    
    bool initVector();
    std::vector<std::size_t> precisionChanges_;
    std::vector<std::chrono::high_resolution_clock::time_point> precisionChanesTimePoint_;
    
    bool vectorInitialized_;
    bool matrixInitialized_;
    
    IndexType usedIterations_;
    
    IndexType m_, n_, nnz_;
    SplitType d_matVals_;
    IndexType *d_matColIndices_;
    IndexType *d_matRowStartIndices_;
    IndexType *d_matEmptyColIndices_;
    IndexType matEmptyColNumber_;

    OriginalPointerType *d_buffer_;
    
    static constexpr std::size_t numberXVecs = 2;
    std::size_t curXIdx_;
    SplitType d_x_[numberXVecs];
    bool isFormatSwitched_;
    
    MySegmentsNumberType readPrec_, writePrec_;
    
    SplitType &getNewX(std::size_t newVersionNumber) {
        std::size_t newIdx = (curXIdx_ + newVersionNumber) % numberXVecs;
        return d_x_[newIdx];
    }
    SplitType &getOldX(std::size_t oldVersionNumber) {
        std::size_t oldIdx = (curXIdx_ + numberXVecs - oldVersionNumber % numberXVecs) % numberXVecs;
        return d_x_[oldIdx];
    }
    SplitType &getCurrentX() {
        return d_x_[curXIdx_];
    }
    void IncreaseXIdx() {
        curXIdx_ = (curXIdx_ + 1) % numberXVecs;
    }
    
    std::vector<OriginalPointerType> x_;
    bool resultUpdated_;
};


#include "csr_page_rank.inc.cu"

#undef SPLIT_POINTER_TEMPLATE_HEAD
#undef SPLIT_POINTER_TEMPLATES


#endif // CSR_PAGE_RANK_CUH
