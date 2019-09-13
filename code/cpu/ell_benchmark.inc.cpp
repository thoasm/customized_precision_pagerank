template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::newVectors(std::size_t vecSize)
{
    if (vecInitialized_ && vecN_ != vecSize) {
        clearVec();
    }
    if (!vecInitialized_) {
        newSplitPtr(&x_, vecSize);
        newSplitPtr(&y_, vecSize);
    }
    vecN_ = vecSize;
    vecInitialized_ = true;
}


template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setUpVectors(std::size_t n)
{
    if (vecInitialized_ && vecN_ < n) {
        clearVec();
    }
    if (!vecInitialized_) {
        newVectors(n);
        setVectors(static_cast<OriginalPointerType>(1) / static_cast<OriginalPointerType>(n));
    }
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setUpVectors(const std::vector<OriginalPointerType> &n)
{
    if (vecInitialized_ && vecN_ < n.size()) {
        clearVec();
    }
    if (!vecInitialized_) {

        newVectors(n.size());
        if (!vecInitialized_) {
            return;
        }
        OriginalPointerType sum = 0;
        for (std::size_t i = 0; i < n.size(); ++i) {
            x_.writeAll(i, n[i]);
            y_.writeAll(i, n[i]);
        }
    }
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setVectors(OriginalPointerType initValues)
{
    if (!vecInitialized_) {
        return;
    }

    kernels::setVector(vecN_, x_, initValues);
    kernels::setVector(vecN_, y_, initValues);
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setUpMatrix(const CooMatrix<IndexType, OriginalPointerType> &mtx)
{
    if (matInitialized_) {
        clearMat();
    }
    if (!matInitialized_) {
        EllMatrix<IndexType, OriginalPointerType> ellMtx(mtx);

        setUpMatrix(ellMtx);
    }
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setUpMatrix(const EllMatrix<IndexType, OriginalPointerType> &mtx)
{
    if (matInitialized_) {
        clearMat();
    }
    if (!matInitialized_ && mtx.isValid()) {
        const std::vector<OriginalPointerType> &entries = mtx.getVals();
        newSplitPtr(&vals_, entries.size());

        for (std::size_t i = 0; i < entries.size(); ++i) {
            vals_.writeAll(i, entries[i]);
        }
        indices_ = mtx.getIndices();
        emptyColIndices_ = mtx.getEmptyColIndices();

        emptyColNumber_ = emptyColIndices_.size();
        m_ = mtx.getM();
        n_ = mtx.getN();
        k_ = mtx.getK();
        if (!vecInitialized_ || vecN_ != n_) {
            setUpVectors(n_);
        }
        matInitialized_ = true;
    }
}



// returns ns for one execution, averaged over <loopCount_> runs
template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber>
std::uint64_t EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runOwnNorm()
{
    if (!vecInitialized_) {
        return 0;
    }
    // run the Iterations at least once before measuring
    OriginalPointerType result{}, lastResult{};
    lastResult = result = norm1<ReadPartsNumber>(x_, vecN_);
    //std::cout << "calculating Norm of size: "<< vecN_ << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        result = norm1<ReadPartsNumber>(x_, vecN_);
        //if (result != lastResult /*|| abs(result - lastResult) > 1e-19*/ ) {
        //    std::cerr << "Error while calculating own norm! " << lastResult << " vs. " << result << "\n";
        //}
        lastResult = result;
    }
    auto end = std::chrono::high_resolution_clock::now();
    normResult_ = result;
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    return totalRuntime / loopCount_;
}



template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber>
std::uint64_t EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runSpmv()
{
    if (!matInitialized_ || !vecInitialized_ || vecN_ < m_ || m_ != n_) {
        std::cout << "Spmv error: " << matInitialized_ << ' ' << vecInitialized_ << ' ' << vecN_ << ' ' << m_ << 'x' << n_ << std::endl;
        return 0;
    }

    // run the Iterations at least once before measuring
    const OriginalPointerType dampingFactor = static_cast<OriginalPointerType>(0.85);
    const OriginalPointerType correction = static_cast<OriginalPointerType>(0);

    kernels::ellPageRank<ReadPartsNumber, ReadPartsNumber>(
        m_,
        k_,
        indices_.data(),
        vals_,
        dampingFactor,
        correction,
        x_,
        y_);
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        auto oldX = x_;
        x_ = y_;
        y_ = oldX;
        kernels::ellPageRank<ReadPartsNumber, ReadPartsNumber>(
            m_,
            k_,
            indices_.data(),
            vals_,
            dampingFactor,
            correction,
            x_,
            y_);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    vectorResult_.resize(m_);
    for (std::size_t i = 0; i < m_; ++i) {
        vectorResult_[i] = y_.readAll(i);
    }
    return totalRuntime / loopCount_;
}



template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber>
std::uint64_t EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runRWTest()
{
    if (!vecInitialized_) {
        return 0;
    }

    //std::cout << "Running RW test with: " << gridRW.x << " x " << blockRW.x << ", Reading " << ReadPartsNumber << " Parts..." << std::endl;
    kernels::rwTest<ReadPartsNumber>(vecN_, x_, 2.0);

    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        kernels::rwTest<ReadPartsNumber>(vecN_, x_, 2.0);
    }
    auto end = std::chrono::high_resolution_clock::now();

    vectorResult_.resize(m_);
    for (std::size_t i = 0; i < m_; ++i) {
        vectorResult_[i] = x_.readAll(i);
    }
    normResult_ = norm1<ReadPartsNumber>(x_, vecN_);

    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    return totalRuntime / loopCount_;
}



template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::testConversion()
{
    if (!vecInitialized_) {
        return;
    }
    constexpr std::size_t doubleValuesIn2Blocks = 8 * memoryBlockSize / sizeof(typename SplitType::SegmentType);
    OriginalPointerType hostVals[doubleValuesIn2Blocks];
    std::size_t testSize = std::min(doubleValuesIn2Blocks, static_cast<std::size_t>(vecN_));

    std::vector<typename SplitType::SegmentType> tempVec(SplitType::getMinimumNumberSegmentsFor(testSize));
    SplitType tempSplit;
    tempSplit.segmentPtr = tempVec.data();

    auto vectorValue = [](std::size_t i) {    return static_cast<OriginalPointerType>(i) + OriginalPointerType{1} / static_cast<OriginalPointerType>(i); };
    for (std::size_t i = 0; i < testSize; ++i) {
        tempSplit.writeAll(i, vectorValue(i));
    }

    kernels::convertToOriginal(testSize, tempSplit);
    OriginalPointerType *dblPtr = reinterpret_cast<OriginalPointerType *>(tempSplit.segmentPtr);

    std::cout << std::scientific << std::setprecision(15);
    for (std::size_t i = 0; i < testSize; ++i) {
        if (dblPtr[i] != vectorValue(i)) {
            std::cout << std::setw(3) << i << ": Expected: " << vectorValue(i)
                << "\n     " << "Actual:   " << dblPtr[i] << '\n';
            std::cerr << "BIG ERROR!!! UNEQUAL!!!\n\n";
        }
    }

    kernels::convertToSplitPointer(testSize, tempSplit);

    for (std::size_t i = 0; i < testSize; ++i) {
        OriginalPointerType curVal = tempSplit.readAll(i);
        if (curVal != vectorValue(i)) {
            std::cout << std::setw(3) << i << ": Expected: " << vectorValue(i)
                << "\n     " << "Actual:   " << curVal << '\n';
            std::cerr << "BIG ERROR!!! UNEQUAL!!!\n\n";
        }
    }
}


template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
std::uint64_t EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runConversion()
{
    if (!vecInitialized_ || !matInitialized_) {
        return 0;
    }

    // run the Iterations at least once before measuring

    for (int readPrec = 2; readPrec < NumberSegments; ++readPrec) {
        normalizeVectorSelector(readPrec, x_, vecN_);
    }
    kernels::convertToOriginal(vecN_, x_);
    kernels::convertToOriginal(k_*m_, vals_);

    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        for (int readPrec = 2; readPrec < NumberSegments; ++readPrec) {
            normalizeVectorSelector(readPrec, x_, vecN_);
        }
        kernels::convertToOriginal(vecN_, x_);
        kernels::convertToOriginal(k_*m_, vals_);
    }
    auto end = std::chrono::high_resolution_clock::now();

    if ((loopCount_ + 1) % 2 == 1) {
        kernels::convertToOriginal(vecN_, x_);
        kernels::convertToOriginal(k_*m_, vals_);
    }
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    return totalRuntime / loopCount_;
}



template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
std::uint64_t EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runNormalization(int readPrec)
{
    if (!vecInitialized_ || !matInitialized_) {
        return 0;
    }

    // run the Iterations at least once before measuring

    normalizeVectorSelector(readPrec, x_, vecN_);

    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        normalizeVectorSelector(readPrec, x_, vecN_);
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    return totalRuntime / loopCount_;
}





template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber>
std::uint64_t EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runDiffNorm()
{
    if (!vecInitialized_) {
        return 0;
    }
    // run the Iterations at least once before measuring
    OriginalPointerType result{}, lastResult{};
    lastResult = result = norm1Diff<ReadPartsNumber>(x_, y_, vecN_);
    //std::cout << "calculating Norm of size: "<< vecN_ << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        result = norm1Diff<ReadPartsNumber>(x_, y_, vecN_);
        //if (result != lastResult) {
        //    std::cerr << "Error while calculating diff norm! " << lastResult << " vs. " << result << "\n";
        //}
        lastResult = result;
    }
    auto end = std::chrono::high_resolution_clock::now();
    normResult_ = result;
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    return totalRuntime / loopCount_;
}


template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber>
std::uint64_t EllBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runSelectiveNorm()
{
    if (!vecInitialized_) {
        return 0;
    }
    // run the Iterations at least once before measuring
    OriginalPointerType result{}, lastResult{};
    lastResult = result = selectiveNorm1<ReadPartsNumber>(x_, emptyColNumber_, emptyColIndices_.data());
    //std::cout << "calculating Norm of size: "<< vecN_ << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        result = selectiveNorm1<ReadPartsNumber>(x_, emptyColNumber_, emptyColIndices_.data());
        //if (result != lastResult) {
        //    std::cerr << "Error while calculating selective norm! " << lastResult << " vs. " << result << "\n";
        //}
        lastResult = result;
    }
    auto end = std::chrono::high_resolution_clock::now();
    normResult_ = result;
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    return totalRuntime / loopCount_;
}
