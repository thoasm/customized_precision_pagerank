

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::clearMtx()
{
    if (matrixInitialized_) {
        matColIndices_.clear();
        matRowStartIndices_.clear();
        matEmptyColIndices_.clear();
        
        deleteSplitPtr(matVals_);
        matEmptyColNumber_ = 0;
    }
    matrixInitialized_ = false;
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::clearVector()
{
    if (vectorInitialized_) {
        for (std::size_t i = 0; i < numberXVecs; ++i) {
            deleteSplitPtr(x_[i]);
        }
    }
    result_.clear();
    vectorInitialized_ = false;
}


template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
bool CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::initMatrix(const CsrMatrix<IndexType, OriginalPointerType> &csrMtx)
{
    if (!csrMtx.isValid() || csrMtx.getVals().size() != csrMtx.getColIndices().size()
        || csrMtx.getM() != csrMtx.getN()
        || csrMtx.getRowStartIndices().size() != (csrMtx.getM() + 1)) {
        return false;
    }
    
    int status = 0;
    
    clearMtx();
    
    m_ = csrMtx.getM();
    n_ = csrMtx.getN();
    nnz_ = csrMtx.getNnz();
    
    const IndexType tempSize = std::max(nnz_, n_);
    
    // Copy matrix to device
    newSplitPtr(&matVals_, nnz_);
    
    for (size_t i = 0; i < nnz_; ++i) {
        matVals_.writeAll(i, csrMtx.getVals()[i]);
    }
    
    matColIndices_ = csrMtx.getColIndices();
    matRowStartIndices_ = csrMtx.getRowStartIndices();
    matEmptyColIndices_ = csrMtx.getEmptyColIndices();
    matEmptyColNumber_ = matEmptyColIndices_.size();
    
    
    matrixInitialized_ = true;
    isFormatSwitched_ = false;
    
    if (initVector()) {
        return matrixInitialized_;
    } else {
        clearMtx();
        return matrixInitialized_;
    }
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
bool CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::initVector()
{
    if (!matrixInitialized_) {
        return false;
    }
    clearVector();
    
    for (std::size_t i = 0; i < numberXVecs; ++i) {
        newSplitPtr(&(x_[i]), n_);
    }
    
    vectorInitialized_ = true;
    resetVector();
    
    return vectorInitialized_;
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::resetVector()
{
    if (!vectorInitialized_) {
        return;
    }
    kernels::setVector(n_, getCurrentX(), static_cast<OriginalPointerType>(1) / static_cast<OriginalPointerType>(n_));
    
    if (isFormatSwitched_ && matrixInitialized_) {
        kernels::convertToSplitPointer(nnz_, matVals_);
    }
    isFormatSwitched_ = false;
}


namespace csr {
namespace detail {

template<int RwParameter, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void pageRankIteration(
    macros::compile_int_list<RwParameter>,
    IndexType m,
    IndexType n,
    const IndexType *matColIndices,
    const IndexType *matRowStartIndices,
    const SplitPointer<NumberSegments, OriginalPointerType> &matVals,
    IndexType matNumberEmptyCols,
    const IndexType *matEmptyColIndices,
    OriginalPointerType dampingFactor,
    const SplitPointer<NumberSegments, OriginalPointerType> &xOld,
    const SplitPointer<NumberSegments, OriginalPointerType> &xNew)
{
    constexpr MySegmentsNumberType readSegments_ = (RwParameter / 2 < NumberSegments) ? RwParameter / 2 : NumberSegments;
    constexpr MySegmentsNumberType writeSegments_ = (readSegments_ == NumberSegments || (RwParameter & int{1} == 0)) ? readSegments_ : readSegments_ + 1;
    
    OriginalPointerType correction = static_cast<OriginalPointerType>(0);
    
    if (matNumberEmptyCols > 0) {
        const OriginalPointerType selectiveVecNorm = selectiveNorm1<readSegments_>(xOld, matNumberEmptyCols, matEmptyColIndices);
        correction = selectiveVecNorm / static_cast<OriginalPointerType>(n);
    }
    kernels::csrPageRank<readSegments_, writeSegments_>(
        m,
        n,
        matColIndices,
        matRowStartIndices,
        matVals,
        dampingFactor,
        correction,
        xOld,
        xNew);
}

GKO_ENABLE_IMPLEMENTATION_SELECTION(selectPageRankIteration, pageRankIteration);


}   // namespace detail
}   // namespace csr




template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::solve(OriginalPointerType dampingFactor, OriginalPointerType epsilon, std::size_t maxIterations)
{
    if (!isInit()) {
        return;
    }
    using intParameterList = macros::compile_int_list<>;
    using ConvergedRWParameterList = macros::compile_int_list<1+1>;
    using ConvergedRParameterList = macros::compile_int_list<1>;
    using TypeArgs = macros::compile_type_list<>;
    
    using ConvertedSplitType = SplitPointer<1, OriginalPointerType>;
    
    bool doVectorNormalization{false};
    
    isFormatSwitched_ = false;
    
    OriginalPointerType lastDiff{1};
    OriginalPointerType normDiff{1};

    //OriginalPointerType startingDiffNorm{};
    //OriginalPointerType relativeGrowth{};
    
    // start with full precision to determine the most accurate norm differences
    writePrec_ = 1; //NumberSegments;
    readPrec_ = writePrec_;
    
    precisionChanges_[0] = 0;
    precisionChanesTimePoint_[0] = std::chrono::high_resolution_clock::now();
    
    for (usedIterations_ = 0; usedIterations_ < maxIterations; ++usedIterations_, IncreaseXIdx()) {
        
        auto &curX = getCurrentX();
        auto &newX = getNewX(1);
        
        ConvertedSplitType convertedCurX;
        convertedCurX.segmentPtr = reinterpret_cast<OriginalPointerType *>(curX.segmentPtr);
        ConvertedSplitType convertedNewX;
        convertedNewX.segmentPtr = reinterpret_cast<OriginalPointerType *>(newX.segmentPtr);
        
        ConvertedSplitType convertedMatVals;
        convertedMatVals.segmentPtr = reinterpret_cast<OriginalPointerType *>(matVals_.segmentPtr);
        
        if (!isFormatSwitched_) {
            csr::detail::selectPageRankIteration(RangeRW{},
                [this](int rwParameter) { return rwParameter == readPrec_ + writePrec_; },
                intParameterList{},
                TypeArgs{},
                m_,
                n_,
                matColIndices_.data(),
                matRowStartIndices_.data(),
                matVals_,
                matEmptyColNumber_,
                matEmptyColIndices_.data(),
                dampingFactor,
                curX,
                newX);
        } else {
            csr::detail::pageRankIteration(ConvergedRWParameterList{},
                m_,
                n_,
                matColIndices_.data(),
                matRowStartIndices_.data(),
                convertedMatVals,
                matEmptyColNumber_,
                matEmptyColIndices_.data(),
                dampingFactor,
                convertedCurX,
                convertedNewX);
        }

        readPrec_ = writePrec_;
        
        // Normalization is only done when everything is separated, so before the conversion
        if (doVectorNormalization) {
            normalizeVectorSelector(readPrec_, newX, n_);
            doVectorNormalization = false;
        }
        
        if (usedIterations_ > 0 && usedIterations_ % normCheckFrequency == 0) {
            
            if (!isFormatSwitched_) {
                norm1DiffSelector<NumberSegments>(&normDiff, readPrec_, curX, newX, n_);
            } else {
                norm1DiffSelector<1>(&normDiff, readPrec_, convertedCurX, convertedNewX, n_);
            }
            
            if (writePrec_ < NumberSegments) {
                if (lastDiff <= normDiff) {
                    ++writePrec_;
                }
                else if (normDiff <= currentPrecision<OriginalPointerType>(writePrec_, NumberSegments) * 10) {
                    ++writePrec_;
                }
            }
            // if it was increased, store the iteration it changed
            if (writePrec_ > readPrec_) {
                doVectorNormalization = true;
                
                precisionChanges_[writePrec_ - 1] = usedIterations_;
                precisionChanesTimePoint_[writePrec_ - 1] = std::chrono::high_resolution_clock::now();
                //std::cout << "Changed write to " << writePrec_ << " at: " << usedIterations_ << " with diffNorm: " << normDiff << std::endl;
            }
            
            lastDiff = normDiff;
        }
        
        if (NumberSegments > 1 && !isFormatSwitched_ && writePrec_ == NumberSegments && readPrec_ == NumberSegments) {
            //std::cout << "Switching at iteration " << usedIterations_ << std::endl;
            //for (std::size_t i = 0; i < numberXVecs; ++i) {
            //    convertToOriginal<<<gridSwitchVector, blockSwitchVector>>>(n_, d_x_[i]);
            //}
            kernels::convertToOriginal(n_, newX);
            
            kernels::convertToOriginal(nnz_, matVals_);
            isFormatSwitched_ = true;
        }
        
        if (normDiff < epsilon) {
            IncreaseXIdx();
            break;
        }
    }
    // the last part of precisionChanges states the End
    precisionChanges_[NumberSegments] = usedIterations_;
    precisionChanesTimePoint_[NumberSegments] = std::chrono::high_resolution_clock::now();
    
    resultUpdated_ = false;
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
const std::vector<OriginalPointerType> &CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::getResult()
{
    if (!resultUpdated_ && isInit()) {
        
        result_.resize(n_);
        
        if (!isFormatSwitched_) {
            // TODO: include readPrec_ in here!!!
            for (std::size_t i = 0; i < n_; ++i) {
                result_[i] = getCurrentX().readAll(i);
            }
        } else {
            std::memcpy(result_.data(), getCurrentX().segmentPtr, n_*sizeof(OriginalPointerType));
        }
        
        resultUpdated_ = true;
    }
    return result_;
}



