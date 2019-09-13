

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::clearMtx()
{
    if (matrixInitialized_) {
        matIndices_.clear();
        matEmptyColIndices_.clear();
        deleteSplitPtr(matVals_);
        matEmptyColNumber_ = 0;
    }
    matrixInitialized_ = false;
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::clearVector()
{
    if (vectorInitialized_) {
        for (std::size_t i = 0; i < numberXVecs; ++i) {
            deleteSplitPtr((x_[i]));
        }
    }
    vectorInitialized_ = false;
    result_.clear();
}


template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
bool EllPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::initMatrix(const EllMatrix<IndexType, OriginalPointerType> &ellMtx)
{
    if (!ellMtx.isValid() || ellMtx.getVals().size() != ellMtx.getIndices().size()
        || ellMtx.getM() != ellMtx.getN() || ellMtx.getVals().size() != ellMtx.getM() * ellMtx.getK()) {
        return false;
    }
    
    clearMtx();

    m_ = ellMtx.getM();
    n_ = ellMtx.getN();
    k_ = ellMtx.getK();

    const IndexType ellMtxSize = ellMtx.getVals().size();   // == m_ * k_
    //const IndexType tempSize = std::max(ellMtxSize, n_);

    // Set up temp variable (on host system)
    newSplitPtr(&matVals_, ellMtxSize);

    for (size_t i = 0; i < ellMtxSize; ++i) {
        matVals_.writeAll(i, ellMtx.getVals()[i]);
    }
    matIndices_ = ellMtx.getIndices();

    matEmptyColIndices_ = ellMtx.getEmptyColIndices();
    matEmptyColNumber_ = matEmptyColIndices_.size();

    matrixInitialized_ = true;
    isFormatSwitched_ = false;

    return initVector();
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
bool EllPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::initVector()
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
void EllPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::resetVector()
{
    if (!vectorInitialized_) {
        return;
    }
    kernels::setVector(n_, getCurrentX(), static_cast<OriginalPointerType>(1) / static_cast<OriginalPointerType>(n_));

    if (isFormatSwitched_ && matrixInitialized_) {
        kernels::convertToSplitPointer(m_ * k_, matVals_);
    }
    isFormatSwitched_ = false;
}


namespace ell {
namespace detail {

template<int RwParameter, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void pageRankIteration(
    macros::compile_int_list<RwParameter>,
    IndexType vectorSize,
    IndexType matK,
    const IndexType *matIndices,
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
        correction = selectiveVecNorm / static_cast<OriginalPointerType>(vectorSize);
    }
    kernels::ellPageRank<readSegments_, writeSegments_>(
        vectorSize,
        matK,
        matIndices,
        matVals,
        dampingFactor,
        correction,
        xOld,
        xNew);
}

GKO_ENABLE_IMPLEMENTATION_SELECTION(selectPageRankIteration, pageRankIteration);


}   // detail namespace
}   // ell namespace


template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::solve(OriginalPointerType dampingFactor, OriginalPointerType epsilon, std::size_t maxIterations)
{
    if (!isInit()) {
        return;
    }
    using BlockSizeList = macros::compile_int_list<>;
    using ConvergedRWParameterList = macros::compile_int_list<1+1>;
    //using ConvergedRParameterList = macros::compile_int_list<1>;
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
        //convertedCurX.size = curX.size;
        convertedCurX.segmentPtr = reinterpret_cast<OriginalPointerType *>(curX.segmentPtr);
        ConvertedSplitType convertedNewX;
        //convertedNewX.size = newX.size;
        convertedNewX.segmentPtr = reinterpret_cast<OriginalPointerType *>(newX.segmentPtr);

        ConvertedSplitType convertedMatVals;
        //convertedMatVals.size = matVals_.size;
        convertedMatVals.segmentPtr = reinterpret_cast<OriginalPointerType *>(matVals_.segmentPtr);

        if (!isFormatSwitched_) {
            ell::detail::selectPageRankIteration(RangeRW{},
                [this](int rwParameter) { return rwParameter == readPrec_ + writePrec_; },
                BlockSizeList{},
                TypeArgs{},
                m_,
                k_,
                matIndices_.data(),
                matVals_,
                matEmptyColNumber_,
                matEmptyColIndices_.data(),
                dampingFactor,
                curX,
                newX);
        } else {
            ell::detail::pageRankIteration(ConvergedRWParameterList{},
                m_,
                k_,
                matIndices_.data(),
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
            }

            //std::cout << usedIterations_ << ": normDiff: " << normDiff << "; estimated: " << startingDiffNorm * std::pow(relativeGrowth, usedIterations_ - growthCalculationIteration + 1) <<  "; write: " << writePrec_ << std::endl;

            lastDiff = normDiff;
        }

        if (NumberSegments > 1 && !isFormatSwitched_ && writePrec_ == NumberSegments && readPrec_ == NumberSegments) {

            //std::cout << "Switching at iteration " << usedIterations_ << std::endl;

            //for (std::size_t i = 0; i < numberXVecs; ++i) {
            //    convertToOriginal<<<gridSwitchVector, blockSwitchVector>>>(n_, d_x_[i]);
            //}
            kernels::convertToOriginal(n_, newX);

            kernels::convertToOriginal(m_ * k_, matVals_);
            isFormatSwitched_ = true;
        }

        if (normDiff < epsilon) {
            IncreaseXIdx();
            break;
        }
    }
    precisionChanges_[writePrec_] = usedIterations_;
    precisionChanesTimePoint_[writePrec_] = std::chrono::high_resolution_clock::now();

    //std::cout << "Used Iterations: " << usedIterations_ << "; normDiff: " << normDiff << std::endl;
    resultUpdated_ = false;
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
const std::vector<OriginalPointerType> &EllPageRank<SPLIT_POINTER_TEMPLATES, IndexType>::getResult()
{
    if (!resultUpdated_ && isInit()) {
        result_.resize(n_);

        if (!isFormatSwitched_) {
            // TODO: include readPrec_ in here!!!
            for (std::size_t i = 0; i < n_; ++i) {
                result_[i] = getCurrentX().readAll(i);
            }
        } else {
            std::memcpy(result_.data(), getCurrentX().segmentPtr, n_ * sizeof(OriginalPointerType));
        }

        resultUpdated_ = true;
    }
    return result_;
}



