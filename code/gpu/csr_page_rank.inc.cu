

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::clearMtx()
{
    if (matrixInitialized_) {
        cudaFree(d_matColIndices_); d_matColIndices_ = nullptr;
        cudaFree(d_matRowStartIndices_); d_matRowStartIndices_ = nullptr;
        cudaFree(d_matEmptyColIndices_); d_matEmptyColIndices_ = nullptr;
        freeSplitPtr(d_matVals_);
        matEmptyColNumber_ = 0;
    }
    matrixInitialized_ = false;
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::clearVector()
{
    if (vectorInitialized_) {
        for (std::size_t i = 0; i < numberXVecs; ++i) {
            freeSplitPtr((d_x_[i]));
        }
        cudaFree(d_buffer_);
    }
    vectorInitialized_ = false;
    x_.clear();
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
    
    SplitType tmpSplit;
    
    // Set up temp variable (on host system)
    std::vector<std::int8_t> tmpHostSegments(SplitType::getMinimumBytesFor(tempSize), 0);
    tmpSplit.segmentPtr = reinterpret_cast<typename SplitType::SegmentType *>(tmpHostSegments.data());
    
    // Copy matrix to device
    status |= mallocSplitPtr(&d_matVals_, nnz_);
    if (status) {
        return false;
    }
    for (size_t i = 0; i < nnz_; ++i) {
        tmpSplit.writeAll(i, csrMtx.getVals()[i]);
    }
    copy(d_matVals_, tmpSplit, nnz_, cudaMemcpyHostToDevice);
    
    status |= cudaMalloc(&d_matColIndices_, nnz_ * sizeof(IndexType));
    if (status) {
        freeSplitPtr(d_matVals_);
        return false;
    }
    cudaMemcpy(d_matColIndices_, csrMtx.getColIndices().data(), nnz_ * sizeof(IndexType), cudaMemcpyHostToDevice);
    
    status |= cudaMalloc(&d_matRowStartIndices_, (m_ + 1) * sizeof(IndexType));
    if (status) {
        freeSplitPtr(d_matVals_);
        cudaFree(d_matColIndices_);
        return false;
    }
    cudaMemcpy(d_matRowStartIndices_, csrMtx.getRowStartIndices().data(), (m_ + 1) * sizeof(IndexType), cudaMemcpyHostToDevice);
    
    const auto &emptyColIdxs = csrMtx.getEmptyColIndices();
    matEmptyColNumber_ = emptyColIdxs.size();
    
    status |= cudaMalloc(&d_matEmptyColIndices_, matEmptyColNumber_ * sizeof(IndexType));
    if (status) {
        freeSplitPtr(d_matVals_);
        cudaFree(d_matColIndices_);
        cudaFree(d_matRowStartIndices_);
        return false;
    }
    cudaMemcpy(d_matEmptyColIndices_, emptyColIdxs.data(), matEmptyColNumber_ * sizeof(IndexType), cudaMemcpyHostToDevice);
    
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
    int status = 0;
    
    for (std::size_t i = 0; i < numberXVecs; ++i) {
        status |= mallocSplitPtr(&(d_x_[i]), n_);
        if (status) {
            for (std::size_t j = 0; j < i; ++j) {
                freeSplitPtr(d_x_[j]);
            }
            return false;
        }
    }
    status |= cudaMalloc(&d_buffer_, n_ * sizeof(OriginalPointerType));
    if (status) {
        for (std::size_t i = 0; i < numberXVecs; ++i) {
            freeSplitPtr(d_x_[i]);
        }
        return false;
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
    dim3 blockSetVector(blockSize);
    dim3 gridSetVector((n_ - 1) / blockSetVector.x + 1);
    setVector<<<gridSetVector, blockSetVector>>>(n_, getCurrentX(), static_cast<OriginalPointerType>(1) / static_cast<OriginalPointerType>(n_));
    
    if (isFormatSwitched_ && matrixInitialized_) {
        dim3 blockSwitchMatVals(blockSize);
        dim3 gridSwitchMatVals((nnz_ - 1) / blockSwitchMatVals.x + 1);
        convertToSplitPointer<<<gridSwitchMatVals, blockSwitchMatVals>>>(nnz_, d_matVals_);
    }
    isFormatSwitched_ = false;
    cudaDeviceSynchronize();
}


namespace csr {
namespace detail {

template<int BlockSize, int ThreadsPerRow, int RwParameter, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void pageRankIteration(
    macros::compile_int_list<RwParameter>,
    IndexType m,
    IndexType n,
    const IndexType *d_matColIndices,
    const IndexType *d_matRowStartIndices,
    const SplitPointer<NumberSegments, OriginalPointerType> &d_matVals,
    IndexType matNumberEmptyCols,
    const IndexType *d_matEmptyColIndices,
    OriginalPointerType dampingFactor,
    const SplitPointer<NumberSegments, OriginalPointerType> &d_xOld,
    const SplitPointer<NumberSegments, OriginalPointerType> &d_xNew,
    OriginalPointerType *d_buffer)
{
    const dim3 blockPageRank(BlockSize);
    const dim3 gridPageRank((m * ThreadsPerRow - 1) / blockPageRank.x + 1);
    constexpr MySegmentsNumberType readSegments_ = (RwParameter / 2 < NumberSegments) ? RwParameter / 2 : NumberSegments;
    constexpr MySegmentsNumberType writeSegments_ = (readSegments_ == NumberSegments || (RwParameter & int{1} == 0)) ? readSegments_ : readSegments_ + 1;
    
    OriginalPointerType correction = static_cast<OriginalPointerType>(0);
    
    if (matNumberEmptyCols > 0) {
        const OriginalPointerType selectiveVecNorm = cudaSelectiveNorm1<BlockSize, readSegments_>(d_xOld, matNumberEmptyCols, d_matEmptyColIndices, d_buffer);
        correction = selectiveVecNorm / static_cast<OriginalPointerType>(n);
    }
    csrPageRank<ThreadsPerRow, readSegments_, writeSegments_><<<gridPageRank, blockPageRank>>>(
        m,
        n,
        d_matColIndices,
        d_matRowStartIndices,
        d_matVals,
        dampingFactor,
        correction,
        d_xOld,
        d_xNew);
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
    using intParameterList = macros::compile_int_list<blockSize, threadsPerRow>;
    // Since CUDA 11.6, the blockSize must be
    using blockSizeParameterList = macros::compile_int_list<cnst::blockSize>;
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
        convertedMatVals.segmentPtr = reinterpret_cast<OriginalPointerType *>(d_matVals_.segmentPtr);
        
        if (!isFormatSwitched_) {
            csr::detail::selectPageRankIteration(RangeRW{},
                [this](int rwParameter) { return rwParameter == readPrec_ + writePrec_; },
                intParameterList{},
                TypeArgs{},
                m_,
                n_,
                d_matColIndices_,
                d_matRowStartIndices_,
                d_matVals_,
                matEmptyColNumber_,
                d_matEmptyColIndices_,
                dampingFactor,
                curX,
                newX,
                d_buffer_);
        } else {
            csr::detail::pageRankIteration<blockSize, threadsPerRow>(ConvergedRWParameterList{},
                m_,
                n_,
                d_matColIndices_,
                d_matRowStartIndices_,
                convertedMatVals,
                matEmptyColNumber_,
                d_matEmptyColIndices_,
                dampingFactor,
                convertedCurX,
                convertedNewX,
                d_buffer_);
        }

        readPrec_ = writePrec_;
        
        // Normalization is only done when everything is separated, so before the conversion
        if (doVectorNormalization) {
            normalizeVectorSelector<blockSize>(readPrec_, newX, n_, d_buffer_);
            doVectorNormalization = false;
        }
        
        if (usedIterations_ > 0 && usedIterations_ % normCheckFrequency == 0) {
            
            if (!isFormatSwitched_) {
                cudaNorm1DiffSelector<blockSize, NumberSegments>(&normDiff, readPrec_, curX, newX, n_, d_buffer_);
            } else {
                cudaNorm1DiffSelector<blockSize, NumberSegments>(&normDiff, readPrec_, convertedCurX, convertedNewX, n_, d_buffer_);
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
            dim3 blockSwitchVector(blockSize);
            dim3 gridSwitchVector((n_ - 1) / blockSwitchVector.x + 1);
            
            //std::cout << "Switching at iteration " << usedIterations_ << std::endl;
            //for (std::size_t i = 0; i < numberXVecs; ++i) {
            //    convertToOriginal<<<gridSwitchVector, blockSwitchVector>>>(n_, d_x_[i]);
            //}
            convertToOriginal<<<gridSwitchVector, blockSwitchVector>>>(n_, newX);
            
            dim3 blockSwitchMatVals(blockSize);
            dim3 gridSwitchMatVals((nnz_ - 1) / blockSwitchMatVals.x + 1);
            convertToOriginal<<<gridSwitchMatVals, blockSwitchMatVals>>>(nnz_, d_matVals_);
            isFormatSwitched_ = true;
        }
        
        if (normDiff < epsilon) {
            IncreaseXIdx();
            break;
        }
    }
    // the last part of precisionChanges states the End
    cudaDeviceSynchronize();
    precisionChanges_[NumberSegments] = usedIterations_;
    precisionChanesTimePoint_[NumberSegments] = std::chrono::high_resolution_clock::now();
    
    resultUpdated_ = false;
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
const std::vector<OriginalPointerType> &CsrPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::getResult() const
{
    if (!resultUpdated_ && isInit()) {
        x_.resize(n_);
        
        if (!isFormatSwitched_) {
            std::vector<std::uint8_t> tempData(SplitType::getMinimumBytesFor(n_));
            SplitType tempSplit;
            tempSplit.segmentPtr = reinterpret_cast<typename SplitType::SegmentType>(tempData.data());
            
            copy(tempSplit, getCurrentX(), n_);
            // TODO: include readPrec_ in here!!!
            for (std::size_t i = 0; i < n_; ++i) {
                x_[i] = tempSplit.readAll(i);
            }
        } else {
            cudaMemcpy(x_.data(), getCurrentX().segmentPtr, n_ * sizeof(OriginalPointerType));
        }
        
        resultUpdated_ = true;
    }
    return x_;
}



