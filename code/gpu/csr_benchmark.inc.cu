template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::mallocVectors(std::size_t vecSize)
{
    if (vecInitialized_ && vecN_ != vecSize) {
        clearVec();
    }
    if (!vecInitialized_) {
        int status = 0;
        status |= mallocSplitPtr(&d_x_, vecSize);
        if (status) {
            return;
        }
        status |= mallocSplitPtr(&d_y_, vecSize);
        if (status) {
            freeSplitPtr(d_x_);
            return;
        }
        status |= mallocSplitPtr(&d_split_buffer_, vecSize);
        if (status) {
            freeSplitPtr(d_x_);
            freeSplitPtr(d_y_);
            freeSplitPtr(d_split_buffer_);
            return;
        }
        status |= cudaMalloc(&d_buffer_, sizeof(OriginalPointerType) * vecSize);
        if (status) {
            freeSplitPtr(d_x_);
            freeSplitPtr(d_y_);
            freeSplitPtr(d_split_buffer_);
            return;
        }
    }
    vecN_ = vecSize;
    vecInitialized_ = true;
}



template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setUpVectors(std::size_t n)
{
    if (vecInitialized_ && vecN_ < n) {
        clearVec();
    }
    if (!vecInitialized_) {
        mallocVectors(n);
        setVectors(static_cast<OriginalPointerType>(1) / static_cast<OriginalPointerType>(n));
        /*
        std::vector<OriginalPointerType> rndNumbers(n);
        for (size_t i = 0; i < n; ++i) {
            rndNumbers[i] = rand() % maxRndNumber;
        }
        setUpVectors(rndNumbers);
        */
        /*mallocSplitPtr(&d_x_, n);
        mallocSplitPtr(&d_y_, n);
        cudaMalloc(&d_buffer_, sizeof(OriginalPointerType) * n);
        vecN_ = n;
        vecInitialized_ = true;*/
    }
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setUpVectors(const std::vector<OriginalPointerType> &n)
{
    if (vecInitialized_ && vecN_ < n.size()) {
        clearVec();
    }
    if (!vecInitialized_) {
        mallocVectors(n.size());
        if (!vecInitialized_) {
            return;
        }
        //setVector<<<gridSet, blockSet>>>(d_buffer_, static_cast<OriginalPointerType>(rand() % maxRndNumber));
        cudaMemcpy(d_buffer_, n.data(), sizeof(OriginalPointerType) * n.size(), cudaMemcpyHostToDevice);
        
        setUpHostSplitVec(n.size());
        OriginalPointerType sum = 0;
        for (std::size_t i = 0; i < n.size(); ++i) {
            hostSplitVec_.writeAll(i, n[i]);
            sum += n[i];
            /*
            if (i < 10) {
                std::cout << "[" << i << "] vec number = " << n[i] << std::endl;
            }*/
        }
        //std::cout << "Expected Norm: " << sum << '\n';
        copy(d_x_, hostSplitVec_, n.size(), cudaMemcpyHostToDevice);
        copy(d_y_, hostSplitVec_, n.size(), cudaMemcpyHostToDevice);
        
        clearHostSplitVec();
    }
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setVectors(OriginalPointerType initValues)
{
    if (!vecInitialized_) {
        return;
    }
    const dim3 blockSet(512);
    const dim3 gridSet((vecN_ - 1) / blockSet.x + 1);
    
    setVector<<<gridSet, blockSet>>>(vecN_, d_x_, initValues);
    setVector<<<gridSet, blockSet>>>(vecN_, d_y_, initValues);
    cudaDeviceSynchronize();
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setUpMatrix(const CooMatrix<IndexType, OriginalPointerType> &mtx)
{
    if (matInitialized_) {
        clearMat();
    }
    if (!matInitialized_) {
        CsrMatrix<IndexType, OriginalPointerType> csrMtx(mtx);
        
        setUpMatrix(csrMtx);
    }
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::setUpMatrix(const CsrMatrix<IndexType, OriginalPointerType> &mtx)
{
    if (matInitialized_) {
        clearMat();
    }
    if (!matInitialized_ && mtx.isValid()) {
        const std::vector<OriginalPointerType> &entries = mtx.getVals();
        setUpHostSplitVec(entries.size());
        
        for (std::size_t i = 0; i < entries.size(); ++i) {
            hostSplitVec_.writeAll(i, entries[i]);
        }
        int status = 0;
        status |= mallocSplitPtr(&d_vals_, entries.size());
        if (status) {
            std::cerr << "Error while allocating matrix values (size: " << entries.size() << ") with error: " << status << "!\n";
            return;
        }
        copy(d_vals_, hostSplitVec_, entries.size(), cudaMemcpyHostToDevice);
        status |= cudaMalloc(&d_colIndices_, mtx.getColIndices().size() * sizeof(IndexType));
        if (status) {
            freeSplitPtr(d_vals_);
            std::cerr << "Error while allocating column indices (size: " << mtx.getColIndices().size() << ") with error: " << status << "!\n";
            return;
        }
        cudaMemcpy(d_colIndices_, mtx.getColIndices().data(), mtx.getColIndices().size() * sizeof(IndexType), cudaMemcpyHostToDevice);
        
        status |= cudaMalloc(&d_rowStartIndices_, mtx.getRowStartIndices().size() * sizeof(IndexType));
        if (status) {
            freeSplitPtr(d_vals_);
            cudaFree(d_colIndices_);
            std::cerr << "Error while allocating column indices (size: " << mtx.getColIndices().size() << ") with error: " << status << "!\n";
            return;
        }
        cudaMemcpy(d_rowStartIndices_, mtx.getRowStartIndices().data(), mtx.getRowStartIndices().size() * sizeof(IndexType), cudaMemcpyHostToDevice);
        
        status |= cudaMalloc(&d_emptyColIndices_, mtx.getEmptyColIndices().size() * sizeof(IndexType));
        if (status) {
            freeSplitPtr(d_vals_);
            cudaFree(d_colIndices_);
            cudaFree(d_rowStartIndices_);
            std::cerr << "Error while allocating empty column indices (size: " << mtx.getEmptyColIndices().size() << ") with error: " << status << "!\n";
            return;
        }
        cudaMemcpy(d_emptyColIndices_, mtx.getEmptyColIndices().data(), mtx.getEmptyColIndices().size() * sizeof(IndexType), cudaMemcpyHostToDevice);
        
        emptyColNumber_ = mtx.getEmptyColIndices().size();
        m_ = mtx.getM();
        n_ = mtx.getN();
        nnz_ = mtx.getNnz();
        if (!vecInitialized_ || vecN_ != n_) {
            setUpVectors(n_);
        }
        matInitialized_ = true;
        clearHostSplitVec();
    }
}



// returns ns for one execution, averaged over <loopCount_> runs
template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber, unsigned int BlockSize>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runOwnNorm()
{
    if (!vecInitialized_) {
        return 0;
    }
    // run the Iterations at least once before measuring
    OriginalPointerType result{}, lastResult{};
    lastResult = result = cudaNorm1<BlockSize, ReadPartsNumber>(d_x_, vecN_, d_buffer_);
    //std::cout << "calculating Norm of size: "<< vecN_ << std::endl;
    cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        result = cudaNorm1<BlockSize, ReadPartsNumber>(d_x_, vecN_, d_buffer_);
        if (result != lastResult /*|| abs(result - lastResult) > 1e-19*/ ) {
            std::cerr << "Error while calculating own norm! " << lastResult << " vs. " << result << "\n";
        }
        lastResult = result;
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    normResult_ = result;
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    return totalRuntime / loopCount_;
}


// returns ns for one execution, averaged over <loopCount_> runs
template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber, unsigned int BlockSize>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runAltOwnNorm()
{
    if (!vecInitialized_) {
        return 0;
    }
    // run the Iterations at least once before measuring
    OriginalPointerType result{}, lastResult{};
    lastResult = result = cudaNorm1<BlockSize, ReadPartsNumber>(d_x_, vecN_, d_split_buffer_, d_buffer_);
    //std::cout << "calculating Norm of size: "<< vecN_ << std::endl;
    cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        result = cudaNorm1<BlockSize, ReadPartsNumber>(d_x_, vecN_, d_split_buffer_, d_buffer_);
        if (result != lastResult /*|| abs(result - lastResult) > 1e-19*/ ) {
            std::cerr << "Error while calculating own norm! " << lastResult << " vs. " << result << "\n";
        }
        lastResult = result;
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    normResult_ = result;
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    return totalRuntime / loopCount_;
}


//works only if calculated with normal doubles
template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runCuBlasNorm()
{
    if (!vecInitialized_ || NumberSegments != 1 || std::is_same<OriginalPointerType, double>::value == false) {
        return 0;
    }
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "INIT error with code " << status << std::endl;
        /*
        std::cout << "Error List:\n" << CUBLAS_STATUS_NOT_INITIALIZED << '\n'
                << CUBLAS_STATUS_ALLOC_FAILED << '\n'
                << CUBLAS_STATUS_INVALID_VALUE << '\n'
                << CUBLAS_STATUS_ARCH_MISMATCH << '\n'
                << CUBLAS_STATUS_MAPPING_ERROR << '\n'
                << CUBLAS_STATUS_EXECUTION_FAILED << '\n'
                << CUBLAS_STATUS_INTERNAL_ERROR << '\n'
                << CUBLAS_STATUS_NOT_SUPPORTED << '\n'
                << CUBLAS_STATUS_LICENSE_ERROR << '\n';
        */
        return 0;
    }
    cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
    
    OriginalPointerType result{}, lastResult;
    status = cublasDnrm2(handle, vecN_, d_buffer_, 1, &result);
    lastResult = result;
    if (status != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "Norm error with code " << status << std::endl;
        return 0;
    }
    cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        status = cublasDnrm2(handle, vecN_, d_buffer_, 1, &result);
        if (lastResult != result) {
            std::cerr << "Big Error!\n";
        }
        lastResult = result;
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    normResult_ = result;
    cublasDestroy(handle);
    return totalRuntime / loopCount_;
}

template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber, unsigned int BlockSize>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runSpmv()
{
    if (!matInitialized_ || !vecInitialized_ || vecN_ < m_ || m_ != n_) {
        std::cerr << "Matrix initialized: " << matInitialized_ 
        << " Vector Initialized: "<< vecInitialized_ 
        << " Vector size: " << vecN_ 
        << " Matrix size: " << m_ << 'x' << n_ << std::endl;
        return 0;
    }
    
    const dim3 blockPageRank(BlockSize);
    const dim3 gridPageRank((m_ * numberThreadsPerRow - 1) / blockPageRank.x + 1);
    // run the Iterations at least once before measuring
    const OriginalPointerType dampingFactor = static_cast<OriginalPointerType>(0.85);
    const OriginalPointerType correction = static_cast<OriginalPointerType>(0);
    std::uint64_t status = 0;
    csrPageRank<numberThreadsPerRow, ReadPartsNumber, ReadPartsNumber><<<gridPageRank, blockPageRank>>>(
        m_,
        n_,
        d_colIndices_,
        d_rowStartIndices_,
        d_vals_,
        dampingFactor,
        correction,
        d_x_,
        d_y_);
    cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        auto oldX = d_x_;
        d_x_ = d_y_;
        d_y_ = oldX;
        csrPageRank<numberThreadsPerRow, ReadPartsNumber, ReadPartsNumber><<<gridPageRank, blockPageRank>>>(
            m_,
            n_,
            d_colIndices_,
            d_rowStartIndices_,
            d_vals_,
            dampingFactor,
            correction,
            d_x_,
            d_y_);
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    status |= cudaDeviceSynchronize();
    if (status != 0) {
        std::cerr << "Error with spmv: " << status << std::endl;
    }
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    setUpHostSplitVec(m_);
    copy(hostSplitVec_, d_y_, m_, cudaMemcpyDeviceToHost);
    
    vectorResult_.resize(m_);
    for (std::size_t i = 0; i < m_; ++i) {
        vectorResult_[i] = hostSplitVec_.readAll(i);
    }
    return totalRuntime / loopCount_;
}



template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber, unsigned int BlockSize>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runRWTest()
{
    if (!vecInitialized_) {
        return 0;
    }
    const dim3 blockRW(BlockSize);
    const dim3 gridRW((m_ - 1) / blockRW.x + 1);
    
    //std::cout << "Running RW test with: " << gridRW.x << " x " << blockRW.x << ", Reading " << ReadPartsNumber << " Parts..." << std::endl;
    rwTest<ReadPartsNumber><<<gridRW, blockRW>>>(vecN_, d_x_, 2.0);
    
    cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        rwTest<ReadPartsNumber><<<gridRW, blockRW>>>(vecN_, d_x_, 2.0);
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    
    setUpHostSplitVec(m_);
    copy(hostSplitVec_, d_y_, m_, cudaMemcpyDeviceToHost);
    
    vectorResult_.resize(m_);
    for (std::size_t i = 0; i < m_; ++i) {
        vectorResult_[i] = hostSplitVec_.readAll(i);
    }
    normResult_ = cudaNorm1<BlockSize, ReadPartsNumber>(d_x_, vecN_, d_buffer_);
    
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    return totalRuntime / loopCount_;
}



template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
void CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::testConversion()
{
    if (!vecInitialized_) {
        return;
    }
    constexpr std::size_t doubleValuesIn2Blocks = 8 * 128 / sizeof(typename SplitType::SegmentType);
    OriginalPointerType hostVals[doubleValuesIn2Blocks];
    std::size_t testSize = std::min(doubleValuesIn2Blocks, static_cast<std::size_t>(vecN_));
    setUpHostSplitVec(testSize);
    
    auto vectorValue = [](std::size_t i) {    return static_cast<OriginalPointerType>(i) + OriginalPointerType{1} / static_cast<OriginalPointerType>(i); };
    for (std::size_t i = 0; i < testSize; ++i) {
        hostSplitVec_.writeAll(i, vectorValue(i));
    }
    
    copy(d_split_buffer_, hostSplitVec_, m_, cudaMemcpyHostToDevice);
    
    convertToOriginal<<<1, testSize>>>(testSize, d_split_buffer_);
    
    cudaMemcpy(hostVals, d_split_buffer_.segmentPtr, sizeof(OriginalPointerType) * testSize, cudaMemcpyDeviceToHost);
    std::cout << std::scientific << std::setprecision(15);
    for (std::size_t i = 0; i < testSize; ++i) {
        if (hostVals[i] != vectorValue(i)) {
            std::cout << std::setw(3) << i << ": Expected: " << vectorValue(i)
                << "\n     " << "Actual:   " << hostVals[i] << '\n';
            std::cerr << "BIG ERROR!!! UNEQUAL!!!\n\n";
        }
    }
    
    convertToSplitPointer<<<1,testSize>>>(testSize, d_split_buffer_);
    copy(hostSplitVec_, d_split_buffer_, testSize, cudaMemcpyDeviceToHost);
    
    for (std::size_t i = 0; i < testSize; ++i) {
        OriginalPointerType curVal = hostSplitVec_.readAll(i);
        if (curVal != vectorValue(i)) {
            std::cout << std::setw(3) << i << ": Expected: " << vectorValue(i)
                << "\n     " << "Actual:   " << curVal << '\n';
            std::cerr << "BIG ERROR!!! UNEQUAL!!!\n\n";
        }
    }
}


template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<unsigned int BlockSize>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runConversion()
{
    if (!vecInitialized_ || !matInitialized_) {
        return 0;
    }
    
    const dim3 blockVectorConversion(BlockSize);
    const dim3 gridVectorConversion((vecN_ - 1) / blockVectorConversion.x + 1);
    const dim3 blockMatrixConversion(BlockSize);
    const dim3 gridMatrixConversion((nnz_ - 1) / blockMatrixConversion.x + 1);
    
    // run the Iterations at least once before measuring
    std::uint64_t status = 0;
    
    for (int readPrec = 2; readPrec < NumberSegments; ++readPrec) {
        normalizeVectorSelector<BlockSize>(readPrec, d_x_, vecN_, d_buffer_);
    }
    convertToOriginal<<<gridVectorConversion, blockVectorConversion>>>(vecN_, d_x_);
    convertToOriginal<<<gridMatrixConversion, blockMatrixConversion>>>(nnz_, d_vals_);
    
    status |= cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        for (int readPrec = 2; readPrec < NumberSegments; ++readPrec) {
            normalizeVectorSelector<BlockSize>(readPrec, d_x_, vecN_, d_buffer_);
        }
        convertToOriginal<<<gridVectorConversion, blockVectorConversion>>>(vecN_, d_x_);
        convertToOriginal<<<gridMatrixConversion, blockMatrixConversion>>>(nnz_, d_vals_);
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    
    if ((loopCount_ + 1) % 2 == 1) {
        convertToOriginal<<<gridVectorConversion, blockVectorConversion>>>(vecN_, d_x_);
        convertToOriginal<<<gridMatrixConversion, blockMatrixConversion>>>(nnz_, d_vals_);
    }
    status |= cudaDeviceSynchronize();
    if (status != 0) {
        std::cerr << "Error with conversion: " << status << std::endl;
    }
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    return totalRuntime / loopCount_;
}



template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<unsigned int BlockSize>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runNormalization(int readPrec)
{
    if (!vecInitialized_ || !matInitialized_) {
        return 0;
    }
    
    // run the Iterations at least once before measuring
    std::uint64_t status = 0;
    
    normalizeVectorSelector<BlockSize>(readPrec, d_x_, vecN_, d_buffer_);
    
    status |= cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        normalizeVectorSelector<BlockSize>(readPrec, d_x_, vecN_, d_buffer_);
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    
    status |= cudaDeviceSynchronize();
    if (status != 0) {
        std::cerr << "Error with vector normalization: " << status << std::endl;
    }
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    return totalRuntime / loopCount_;
}




template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber, unsigned int BlockSize>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runDiffNorm()
{
    if (!vecInitialized_) {
        return 0;
    }
    // run the Iterations at least once before measuring
    OriginalPointerType result{}, lastResult{};
    lastResult = result = cudaNorm1Diff<BlockSize, ReadPartsNumber>(d_x_, d_y_, vecN_, d_buffer_);
    //std::cout << "calculating Norm of size: "<< vecN_ << std::endl;
    cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        result = cudaNorm1Diff<BlockSize, ReadPartsNumber>(d_x_, d_y_, vecN_, d_buffer_);
        if (result != lastResult) {
            std::cerr << "Error while calculating diff norm! " << lastResult << " vs. " << result << "\n";
        }
        lastResult = result;
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    normResult_ = result;
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    return totalRuntime / loopCount_;
}


template<SPLIT_POINTER_TEMPLATE_HEAD, typename IndexType>
template<int ReadPartsNumber, unsigned int BlockSize>
std::uint64_t CsrBenchmark<SPLIT_POINTER_TEMPLATES, IndexType>::runSelectiveNorm()
{
    if (!vecInitialized_) {
        return 0;
    }
    // run the Iterations at least once before measuring
    OriginalPointerType result{}, lastResult{};
    lastResult = result = cudaSelectiveNorm1<BlockSize, ReadPartsNumber>(d_x_, emptyColNumber_, d_emptyColIndices_, d_buffer_);
    //std::cout << "calculating Norm of size: "<< vecN_ << std::endl;
    cudaDeviceSynchronize();
    auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < loopCount_; ++i) {
        result = cudaSelectiveNorm1<BlockSize, ReadPartsNumber>(d_x_, emptyColNumber_, d_emptyColIndices_, d_buffer_);
        if (result != lastResult) {
            std::cerr << "Error while calculating selective norm! " << lastResult << " vs. " << result << "\n";
        }
        lastResult = result;
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    normResult_ = result;
    std::uint64_t totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    return totalRuntime / loopCount_;
}
