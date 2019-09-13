#include "split_pointer.hpp"
#include "split_pointer_helper.cuh"

#include "gpu_kernels.cuh"
#include "coo_matrix.hpp"
#include "ell_matrix.hpp"
#include "csr_matrix.hpp"
#include "ell_page_rank.cuh"
#include "csr_page_rank.cuh"
#include "ell_benchmark.cuh"
#include "csr_benchmark.cuh"
#include "constants.hpp"

#include <iostream>
#include <ios>
#include <iomanip>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <chrono>
#include <sstream>
#include <type_traits> // for std::is_same & std::conditional


//#define USE_ELL true
//#define USE_ELL false

#if defined(USE_ELL) && USE_ELL
#define MatrixType EllMatrix
#define BenchmarkType EllBenchmark
#define PageRankType EllPageRank
#else
#define MatrixType CsrMatrix
#define BenchmarkType CsrBenchmark
#define PageRankType CsrPageRank
#endif


template<typename IndexType, typename FloatType>
void printEllMtx(const EllMatrix<IndexType, FloatType> &ellMat)
{
    const size_t k = ellMat.getK();
    const size_t m = ellMat.getM();
    const size_t n = ellMat.getN();
    const std::string spacer(" ");
    constexpr size_t noLastColumn = std::numeric_limits<size_t>::max();
    std::cout << "Matrix: \n";
    for (size_t row = 0; row < m; ++row) {
        size_t lastColumn = noLastColumn;
        for (size_t kIdx = 0; kIdx < k; ++kIdx) {
            const IndexType idx = kIdx * m + row;
            const FloatType curVal = ellMat.getVals()[idx];
            const IndexType curIdx = ellMat.getIndices()[idx];
            if (curVal != static_cast<FloatType>(0)) {
                if (curIdx == 0) {
                    std::cout << curVal << spacer;
                }
                for (size_t j = (lastColumn == noLastColumn) ? 0 : (lastColumn + 1); j < curIdx; ++j) {
                    std::cout << static_cast<FloatType>(0) << spacer;
                }
                if (curIdx != 0) {
                    std::cout << curVal << spacer;
                }
                lastColumn = curIdx;
            }
        }
        for (size_t j = (lastColumn == noLastColumn) ? 0 : (lastColumn + 1); j < n; ++j) {
            std::cout << static_cast<FloatType>(0) << spacer;
        }
        std::cout << '\n';
    }
}




template<MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
std::vector<std::uint64_t> runPageRank(PageRankType<NumberSegments, OriginalPointerType, IndexType> &pageRank, OriginalPointerType epsilon, OriginalPointerType dampingFactor, std::size_t maxIters, std::size_t numberMeasurements)
{
    std::vector<std::uint64_t> calculationTime(NumberSegments + 1, 0);
    std::size_t totalCalcIdx = 0;
    
    for (std::size_t measurementLoops = 0; measurementLoops <= numberMeasurements; ++measurementLoops) {
        // Reset x
        pageRank.resetVector();
        
        auto start = std::chrono::high_resolution_clock::now();
        
        pageRank.solve(dampingFactor, epsilon, maxIters);
        //std::cout << "CUDA error: " << cudaDeviceSynchronize() << std::endl;
        auto end = std::chrono::high_resolution_clock::now();
        // ignore the first iteration
        if (measurementLoops > 0) {
            calculationTime[totalCalcIdx] += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            const auto &timePoints = pageRank.getPrecisionTimePoints();
            for (std::size_t i = totalCalcIdx + 1; i < calculationTime.size(); ++i) {
                const std::size_t startIdx = i - 1;
                const std::size_t endIdx = startIdx + 1;
                
                calculationTime[i] += std::chrono::duration_cast<std::chrono::nanoseconds>(timePoints[endIdx] - timePoints[startIdx]).count();
            }
        }
    }
    for (std::size_t i = 0; i < calculationTime.size(); ++i) {
        if (numberMeasurements > 0) {
            calculationTime[i] /= numberMeasurements;
        } else {
            calculationTime[i] = 0;
        }
    }
    return calculationTime;
}


template<MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
std::vector<std::uint64_t> runVectorNormalizations(BenchmarkType<NumberSegments, OriginalPointerType, IndexType> &bench)
{
    std::vector<std::uint64_t> runtimeNormalizations(NumberSegments+1, 0);
    for (int readPrec = 2; readPrec <= NumberSegments; ++readPrec) {
        runtimeNormalizations[readPrec] = bench.runNormalization(readPrec);
    }
    return runtimeNormalizations;
}

int main(int argc, char *argv[])
{
    using IndexType = unsigned int;
    using OriginalPointerType = double;
    
    std::string pageRankVersionStr(USE_ELL ? "ELL" : "CSR");
    
    constexpr unsigned int blockSize{cnst::blockSize};
    constexpr std::size_t numberKernelRunsToAverage = 100;
    constexpr std::size_t numberPageRankRunsToAverage = 10;
    
    constexpr OriginalPointerType epsilonDef = 1e-10;
    const OriginalPointerType dampingFactorDef = 0.85;
    constexpr std::size_t maxItersDef = 1000;
    
    OriginalPointerType epsilon{epsilonDef};
    OriginalPointerType dampingFactor{dampingFactorDef};
    std::size_t maxIters{maxItersDef};

    if (argc < 2) {
        std::cerr << "You need to specify a sparse matrix!\n";
        return 1;
    }
    std::string mtxPath(argv[1]);
    CooMatrix<IndexType, OriginalPointerType> cooMtx;
    cooMtx.readMatrix(mtxPath.c_str());
    std::size_t slashPos = mtxPath.rfind('/');
    std::size_t backSlashPos = mtxPath.rfind('\\');
    
    std::size_t mtxNameStart{0};
    if (slashPos != std::string::npos && backSlashPos != std::string::npos) {
        mtxNameStart = std::max(slashPos, backSlashPos) + 1;
    } else if (slashPos == std::string::npos && backSlashPos == std::string::npos) {
        mtxNameStart = 0;
    } else {
        mtxNameStart = std::min(slashPos, backSlashPos) + 1;
    }
    std::string fileName = mtxPath.substr(mtxNameStart);
    std::string matrixMarketExtension(".mtx");
    std::size_t fileDotPos = fileName.rfind(matrixMarketExtension);
    std::string matrixName;
    if (fileDotPos != std::string::npos) {
        matrixName = fileName.substr(0, fileName.length() - matrixMarketExtension.length());
    } else {
        matrixName = fileName;
    }
    
    if (argc > 2) {
        epsilon = std::stod(std::string(argv[2]));
    }
    if (argc > 3) {
        dampingFactor = std::stod(std::string(argv[3]));
    }
    if (argc > 4) {
        maxIters = std::stoi(std::string(argv[4]));
    }
    
    if (!cooMtx.isValid()) {
        std::cerr << "Matrix " << matrixName << " is not valid!\n";
        return 1;
    }
    
    cooMtx.transpose();
    cooMtx.normalizeColumns();
    cooMtx.sortRowMajor();
    
    /*
    std::cout << "Matrix: " << matrixName
            << "\nDim:    " << cooMtx.getM() << " x " << cooMtx.getN() << "; Nnz: " << cooMtx.getNnz()
            << "\nNumber empty Cols: " << cooMtx.getEmptyColumnIndices().size() << '\n' << std::endl;
    return 0;
    */
    
    MatrixType<IndexType, OriginalPointerType> sparseMatrix(cooMtx);
    if (!sparseMatrix.isValid()) {
        std::cerr << "Matrix " << matrixName << " in " << pageRankVersionStr << " format is not valid!\n";
        return 1;
    }
    
    
        
    std::stringstream ssSettingsKernels;
    ssSettingsKernels << std::scientific << "% all kernel runtimes averaged over " << numberKernelRunsToAverage
            << " executions (and measured in ns)";
    std::string basicBlocksHeader("data_basicblocks = [");
    std::string headlineSeparator("  &  ");
    std::string dataSeparator("  ");
    std::string blocksEnd("];");
    
    std::string normHeader("data_norms = [");
    
    std::stringstream ssBasicBlocksHeader;
    std::stringstream ssBasicBlocksValues;
    
    std::stringstream ssNormHeader;
    std::stringstream ssNormValues;
    
    ssBasicBlocksHeader << std::scientific;
    //ssBasicBlocksHeader << "% " << "matrix\n% ";
    ssBasicBlocksValues << "% " << matrixName << '\n';
    
    ssNormValues << "% " << matrixName << '\n';
    
    BenchmarkType<1, OriginalPointerType, IndexType> benchDouble(numberKernelRunsToAverage);
    benchDouble.setUpMatrix(sparseMatrix);
    
    ssBasicBlocksHeader << "% " << "read/write double";
    ssBasicBlocksValues << benchDouble.runRWTest<1, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "SpMV " << pageRankVersionStr << " double";
    ssBasicBlocksValues << dataSeparator << benchDouble.runSpmv<1, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "norm double";
    ssBasicBlocksValues << dataSeparator << benchDouble.runOwnNorm<1, blockSize>();
    
    ssNormHeader << "% " << "double selective-norm";
    ssNormValues << benchDouble.runSelectiveNorm<1, blockSize>();
    ssNormHeader << headlineSeparator << "double diff-norm";
    ssNormValues << dataSeparator << benchDouble.runDiffNorm<1, blockSize>();
    
    benchDouble.clear();
    
    
    BenchmarkType<2, OriginalPointerType, IndexType> bench2Segment(numberKernelRunsToAverage);
    bench2Segment.setUpMatrix(sparseMatrix);
    
    ssBasicBlocksHeader << headlineSeparator << "inplace conversion 2segment";
    auto bench2SegmentConversion = bench2Segment.runConversion<blockSize>();
    ssBasicBlocksValues << dataSeparator << bench2SegmentConversion;
    
    bench2Segment.clear();
    bench2Segment.setUpMatrix(sparseMatrix);
    ssBasicBlocksHeader << headlineSeparator << "read/write 2segment 32 bit";
    ssBasicBlocksValues << dataSeparator << bench2Segment.runRWTest<1, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "read/write 2segment 64 bit";
    ssBasicBlocksValues << dataSeparator << bench2Segment.runRWTest<2, blockSize>();
    
    ssBasicBlocksHeader << headlineSeparator << "SpMV " << pageRankVersionStr << " 2segment 32 bit";
    ssBasicBlocksValues << dataSeparator << bench2Segment.runSpmv<1, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "SpMV " << pageRankVersionStr << " 2segment 64 bit";
    ssBasicBlocksValues << dataSeparator << bench2Segment.runSpmv<2, blockSize>();
    
    ssBasicBlocksHeader << headlineSeparator << "norm 2segment 32 bit";
    ssBasicBlocksValues << dataSeparator << bench2Segment.runOwnNorm<1, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "norm 2segment 64 bit";
    ssBasicBlocksValues << dataSeparator << bench2Segment.runOwnNorm<2, blockSize>();
    
    ssNormHeader << headlineSeparator << "2segment 32 bit selective-norm";
    ssNormValues << dataSeparator << bench2Segment.runSelectiveNorm<1, blockSize>();
    ssNormHeader << headlineSeparator << "2segment 64 bit selective-norm";
    ssNormValues << dataSeparator << bench2Segment.runSelectiveNorm<2, blockSize>();
    
    ssNormHeader << headlineSeparator << "2segment 32 bit diff-norm";
    ssNormValues << dataSeparator << bench2Segment.runDiffNorm<1, blockSize>();
    ssNormHeader << headlineSeparator << "2segment 64 bit diff-norm";
    ssNormValues << dataSeparator << bench2Segment.runDiffNorm<2, blockSize>();
    
    auto bench2SegmentNormalizationRuntime = runVectorNormalizations(bench2Segment);

    bench2Segment.clear();
    
    BenchmarkType<4, OriginalPointerType, IndexType> bench4Segment(numberKernelRunsToAverage);
    bench4Segment.setUpMatrix(sparseMatrix);
    
    ssBasicBlocksHeader << headlineSeparator << "inplace conversion 4segment";
    auto bench4SegmentConversion = bench4Segment.runConversion<blockSize>();
    ssBasicBlocksValues << dataSeparator << bench4SegmentConversion;
    
    bench4Segment.clear();
    bench4Segment.setUpMatrix(sparseMatrix);
    ssBasicBlocksHeader << headlineSeparator << "read/write 4segment 16 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runRWTest<1, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "read/write 4segment 32 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runRWTest<2, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "read/write 4segment 48 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runRWTest<3, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "read/write 4segment 64 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runRWTest<4, blockSize>();
    
    ssBasicBlocksHeader << headlineSeparator << "SpMV " << pageRankVersionStr << " 4segment 16 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runSpmv<1, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "SpMV " << pageRankVersionStr << " 4segment 32 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runSpmv<2, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "SpMV " << pageRankVersionStr << " 4segment 48 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runSpmv<3, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "SpMV " << pageRankVersionStr << " 4segment 64 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runSpmv<4, blockSize>();
    
    ssBasicBlocksHeader << headlineSeparator << "norm 4segment 16 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runOwnNorm<1, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "norm 4segment 32 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runOwnNorm<2, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "norm 4segment 48 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runOwnNorm<3, blockSize>();
    ssBasicBlocksHeader << headlineSeparator << "norm 4segment 64 bit";
    ssBasicBlocksValues << dataSeparator << bench4Segment.runOwnNorm<4, blockSize>();
    
    ssNormHeader << headlineSeparator << "4segment 16 bit selective-norm";
    ssNormValues << dataSeparator << bench4Segment.runSelectiveNorm<1, blockSize>();
    ssNormHeader << headlineSeparator << "4segment 32 bit selective-norm";
    ssNormValues << dataSeparator << bench4Segment.runSelectiveNorm<2, blockSize>();
    ssNormHeader << headlineSeparator << "4segment 48 bit selective-norm";
    ssNormValues << dataSeparator << bench4Segment.runSelectiveNorm<3, blockSize>();
    ssNormHeader << headlineSeparator << "4segment 64 bit selective-norm";
    ssNormValues << dataSeparator << bench4Segment.runSelectiveNorm<4, blockSize>();
    
    ssNormHeader << headlineSeparator << "4segment 16 bit diff-norm";
    ssNormValues << dataSeparator << bench4Segment.runDiffNorm<1, blockSize>();
    ssNormHeader << headlineSeparator << "4segment 32 bit diff-norm";
    ssNormValues << dataSeparator << bench4Segment.runDiffNorm<2, blockSize>();
    ssNormHeader << headlineSeparator << "4segment 48 bit diff-norm";
    ssNormValues << dataSeparator << bench4Segment.runDiffNorm<3, blockSize>();
    ssNormHeader << headlineSeparator << "4segment 64 bit diff-norm";
    ssNormValues << dataSeparator << bench4Segment.runDiffNorm<4, blockSize>();
    
    auto bench4SegmentNormalizationRuntime = runVectorNormalizations(bench4Segment);
    
    bench4Segment.clear();
    
    std::cout << ssSettingsKernels.str() << "\n\n";
    std::cout << basicBlocksHeader << '\n'
            << ssBasicBlocksHeader.str()  << '\n'
            << ssBasicBlocksValues.str() << '\n'
            << blocksEnd << '\n' << std::endl;
    
    std::cout << normHeader << '\n'
            << ssNormHeader.str() << '\n'
            << ssNormValues.str() << '\n'
            << blocksEnd << std::endl;
    
    std::cout << "data_normalization_conversion = [\n"
            << "% 2segment 32->64 bit" << headlineSeparator << "2segment in-place conversion" << headlineSeparator << "4segment 16->32 bit" << headlineSeparator << "4segment 32->48 bit" << headlineSeparator << "4segment 48->64 bit" << headlineSeparator << "4segment in-place conversion" << '\n'
            << "% " << matrixName << '\n'
            << bench2SegmentNormalizationRuntime[2] << dataSeparator << bench2SegmentConversion << dataSeparator << bench4SegmentNormalizationRuntime[2] << dataSeparator << bench4SegmentNormalizationRuntime[3] << dataSeparator << bench4SegmentNormalizationRuntime[4] << dataSeparator << bench4SegmentConversion << '\n'
            << blocksEnd << std::endl;
    //return 0;
    
    std::stringstream ssSettingsPageRank;
    ssSettingsPageRank << std::scientific << "% all PageRank runtimes averaged over " << numberPageRankRunsToAverage
            << " executions (and measured in ns), PageRank executed with damping factor of " << dampingFactor
            << " and an epsilon of " << epsilon << "; Sparse format: " << pageRankVersionStr;
    std::string oneIterationHeader("data_oneiteration = [");
    std::stringstream ssOneIterationHeader;
    std::stringstream ssOneIterationValues;
    std::string pageRankHeader("data_pagerank = [");
    std::stringstream ssPageRankHeader;
    std::stringstream ssPageRankValues;
    
    
    PageRankType<1, OriginalPointerType, IndexType> pageDouble;
    pageDouble.initMatrix(sparseMatrix);
    auto pageDoubleRuntime = runPageRank(pageDouble, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
    //ssOneIterationHeader << "% " << "matrix\n% ";
    ssOneIterationValues << "% " << matrixName << '\n';
    
    //ssPageRankHeader << "% " << "matrix\n% ";
    ssPageRankValues << "% " << matrixName << '\n';
    
    ssOneIterationHeader << "% " << "runtime/iter double";
    ssOneIterationValues << pageDoubleRuntime[0] / pageDouble.getUsedIterations();
    
    ssPageRankHeader << "% " << "total iterations double";
    ssPageRankValues << pageDouble.getUsedIterations();
    
    ssPageRankHeader << headlineSeparator << "total runtime PageRank double";
    ssPageRankValues << dataSeparator << pageDoubleRuntime[0];
    
    pageDouble.clear();
    
    
    PageRankType<2, OriginalPointerType, IndexType> page2Segment;
    page2Segment.initMatrix(sparseMatrix);
    auto page2SegmentRuntime = runPageRank(page2Segment, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
    auto page2IterationCount = page2Segment.getPrecisionChanges();
    
    ssOneIterationHeader << headlineSeparator << "runtime/iter 2-segment 32 bit";
    ssOneIterationValues << dataSeparator << page2SegmentRuntime[1] / (page2IterationCount[1] + 1);
    
    ssOneIterationHeader << headlineSeparator << "runtime/iter 2-segment 64 bit";
    if (page2IterationCount[2] - page2IterationCount[1] > 0) {
        ssOneIterationValues << dataSeparator << (page2SegmentRuntime[2] - bench2SegmentNormalizationRuntime[2] - bench2SegmentConversion) / (page2IterationCount[2] - page2IterationCount[1]);
    } else {
        ssOneIterationValues << dataSeparator << 0;
    }
    
    ssPageRankHeader << headlineSeparator << "total iterations 2segment";
    ssPageRankValues << dataSeparator << page2Segment.getUsedIterations();
    
    ssPageRankHeader << headlineSeparator << "switch point 2segment";
    ssPageRankValues << dataSeparator << page2IterationCount[1];
    
    ssPageRankHeader << headlineSeparator << "total runtime 2segment";
    ssPageRankValues << dataSeparator << page2SegmentRuntime[0];

    page2Segment.clear();
    
    
    PageRankType<4, OriginalPointerType, IndexType> page4Segment;
    page4Segment.initMatrix(sparseMatrix);
    auto page4SegmentRuntime = runPageRank(page4Segment, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
    auto page4IterationCount = page4Segment.getPrecisionChanges();
    
    ssOneIterationHeader << headlineSeparator << "runtime/iter 4-segment 16 bit";
    ssOneIterationValues << dataSeparator << page4SegmentRuntime[1] / (page4IterationCount[1] + 1);
    
    ssOneIterationHeader << headlineSeparator << "runtime/iter 4-segment 32 bit";
    if (page4IterationCount[2] - page4IterationCount[1] > 0) {
    ssOneIterationValues << dataSeparator << (page4SegmentRuntime[2] - bench4SegmentNormalizationRuntime[2]) / (page4IterationCount[2] - page4IterationCount[1]);
    } else {
        ssOneIterationValues << dataSeparator << 0;
    }
    
    ssOneIterationHeader << headlineSeparator << "runtime/iter 4-segment 48 bit";
    if (page4IterationCount[3] - page4IterationCount[2] > 0) {
        ssOneIterationValues << dataSeparator << (page4SegmentRuntime[3] - bench4SegmentNormalizationRuntime[3]) / (page4IterationCount[3] - page4IterationCount[2]);
    } else {
        ssOneIterationValues << dataSeparator << 0;
    }
    
    ssOneIterationHeader << headlineSeparator << "runtime/iter 4-segment 64 bit";
    if (page4IterationCount[4] - page4IterationCount[3] > 0) {
        ssOneIterationValues << dataSeparator << (page4SegmentRuntime[4] - bench4SegmentNormalizationRuntime[4] - bench4SegmentConversion) / (page4IterationCount[4] - page4IterationCount[3]);
    } else {
        ssOneIterationValues << dataSeparator << 0;
    }
    
    ssPageRankHeader << headlineSeparator << "total iterations 4segment";
    ssPageRankValues << dataSeparator << page4Segment.getUsedIterations();
    
    ssPageRankHeader << headlineSeparator << "switch point 4segment 16->32";
    ssPageRankValues << dataSeparator << page4IterationCount[1];

    ssPageRankHeader << headlineSeparator << "switch point 4segment 32->48";
    ssPageRankValues << dataSeparator << page4IterationCount[2];
    
    ssPageRankHeader << headlineSeparator << "switch point 4segment 48->64";
    ssPageRankValues << dataSeparator << page4IterationCount[3];
    
    ssPageRankHeader << headlineSeparator << "total runtime 4segment";
    ssPageRankValues << dataSeparator << page4SegmentRuntime[0];

    page4Segment.clear();
    
    std::cout << ssSettingsPageRank.str() << "\n\n"
            << oneIterationHeader << '\n'
            << ssOneIterationHeader.str() << '\n'
            << ssOneIterationValues.str() << '\n'
            << blocksEnd << "\n\n";
    
    std::cout << pageRankHeader << '\n'
            << ssPageRankHeader.str() << '\n'
            << ssPageRankValues.str() << '\n'
            << blocksEnd << "\n\n";
    
    return 0;
}

