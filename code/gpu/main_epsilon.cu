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
#include <cmath>
#include <type_traits> // for std::is_same & std::conditional


//#define USE_ELL true
/*
#define USE_ELL false
#if USE_ELL
#define MatrixType EllMatrix
#define BenchmarkType EllBenchmark
#define PageRankType EllPageRank
#else
#define MatrixType CsrMatrix
#define BenchmarkType CsrBenchmark
#define PageRankType CsrPageRank
#endif
*/


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
std::vector<std::uint64_t> runPageRank(EllPageRank<NumberSegments, OriginalPointerType, IndexType> &pageRank, OriginalPointerType epsilon, OriginalPointerType dampingFactor, std::size_t maxIters, std::size_t numberMeasurements)
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
std::vector<std::uint64_t> runPageRank(CsrPageRank<NumberSegments, OriginalPointerType, IndexType> &pageRank, OriginalPointerType epsilon, OriginalPointerType dampingFactor, std::size_t maxIters, std::size_t numberMeasurements)
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



int main(int argc, char *argv[])
{
    using IndexType = unsigned int;
    using OriginalPointerType = double;
    
    //constexpr unsigned int blockSize{cnst::blockSize};
    //constexpr std::size_t numberKernelRunsToAverage = 100;
    constexpr std::size_t numberPageRankRunsToAverage = 10;
    
    //constexpr OriginalPointerType epsilonDef = 1e-10;
    const OriginalPointerType dampingFactorDef = 0.85;
    constexpr std::size_t maxItersDef = 1000;
    
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
    
    //if (argc > 2) {
    //    epsilon = std::stod(std::string(argv[2]));
    //}
    if (argc > 2) {
        dampingFactor = std::stod(std::string(argv[2]));
    }
    if (argc > 3) {
        maxIters = std::stoi(std::string(argv[3]));
    }
    
    if (!cooMtx.isValid()) {
        std::cerr << "Matrix " << matrixName << " is not valid!\n";
        return 1;
    }
    
    cooMtx.transpose();
    cooMtx.normalizeColumns();
    cooMtx.sortRowMajor();
    
    std::vector<OriginalPointerType> epsList;
    
    OriginalPointerType base10{10};
    for (int pow10 = -2; pow10 >= -12; pow10 -= 2) {
        epsList.push_back(pow(base10, static_cast<OriginalPointerType>(pow10)));
    }
    
    /*
    std::cout << "Matrix: " << matrixName
            << "\nDim:    " << cooMtx.getM() << " x " << cooMtx.getN() << "; Nnz: " << cooMtx.getNnz()
            << "\nNumber empty Cols: " << cooMtx.getEmptyColumnIndices().size() << '\n' << std::endl;
    return 0;
    */
    
    std::string headlineSeparator("  ");
    std::string dataSeparator("  ");
    std::string blocksEnd("];");
    
    bool skip = false;
    CsrMatrix<IndexType, OriginalPointerType> csrMatrix(cooMtx);
    if (!csrMatrix.isValid()) {
        std::cerr << "Matrix " << matrixName << " in CSR format is not valid!\n";
        skip = true;
    }
    
    if (!skip) {
        std::stringstream ssSettingsPageRankCsr;
        ssSettingsPageRankCsr << std::scientific << "% all PageRank runtimes averaged over " << numberPageRankRunsToAverage
                << " executions (and measured in ns), PageRank executed with damping factor of " << dampingFactor
                << "; Sparse format: CSR";
        
        std::cout << ssSettingsPageRankCsr.str() << "\n\n";
        
        std::string pageRankCsr("runtime_pagerank_1_csr = [");
        std::stringstream ssPageRank1CsrHeader;
        ssPageRank1CsrHeader << std::scientific;
        std::stringstream ssPageRank1CsrValues;
        
        CsrPageRank<1, OriginalPointerType, IndexType> pageCsrDouble;
        pageCsrDouble.initMatrix(csrMatrix);
        
        ssPageRank1CsrHeader << "%";
        ssPageRank1CsrValues << "% " << matrixName << '\n';
        for (OriginalPointerType& epsilon : epsList) {
            ssPageRank1CsrHeader << headlineSeparator << epsilon;
            auto pageRuntime = runPageRank(pageCsrDouble, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
            ssPageRank1CsrValues << dataSeparator << pageRuntime[0];
        }
        
        std::cout << pageRankCsr << '\n'
                << ssPageRank1CsrHeader.str() << '\n'
                << ssPageRank1CsrValues.str() << '\n'
                << blocksEnd << '\n';
        pageCsrDouble.clear();
        
        
        
        std::string pageRankCsr2("runtime_pagerank_2_csr = [");
        std::stringstream ssPageRank2CsrHeader;
        ssPageRank2CsrHeader << std::scientific;
        std::stringstream ssPageRank2CsrValues;
        
        CsrPageRank<2, OriginalPointerType, IndexType> pageCsr2;
        pageCsr2.initMatrix(csrMatrix);
        
        ssPageRank2CsrHeader << "%";
        ssPageRank2CsrValues << "% " << matrixName << '\n';
        for (OriginalPointerType& epsilon : epsList) {
            ssPageRank2CsrHeader << headlineSeparator << epsilon;
            auto pageRuntime = runPageRank(pageCsr2, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
            ssPageRank2CsrValues << dataSeparator << pageRuntime[0];
        }
        
        std::cout << pageRankCsr2 << '\n'
                << ssPageRank2CsrHeader.str() << '\n'
                << ssPageRank2CsrValues.str() << '\n'
                << blocksEnd << '\n';
        pageCsr2.clear();
        
        
        std::string pageRankCsr4("runtime_pagerank_4_csr = [");
        std::stringstream ssPageRank4CsrHeader;
        ssPageRank4CsrHeader << std::scientific;
        std::stringstream ssPageRank4CsrValues;
        
        CsrPageRank<4, OriginalPointerType, IndexType> pageCsr4;
        pageCsr4.initMatrix(csrMatrix);
        
        ssPageRank4CsrHeader << "%";
        ssPageRank4CsrValues << " % " << matrixName << '\n';
        for (OriginalPointerType& epsilon : epsList) {
            ssPageRank4CsrHeader << headlineSeparator << epsilon;
            auto pageRuntime = runPageRank(pageCsr4, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
            ssPageRank4CsrValues << dataSeparator << pageRuntime[0];
        }
        
        std::cout << pageRankCsr4 << '\n'
                << ssPageRank4CsrHeader.str() << '\n'
                << ssPageRank4CsrValues.str() << '\n'
                << blocksEnd << '\n';
        pageCsr4.clear();
        
    }
    
    csrMatrix.clear();
    
    
    skip = false;
    EllMatrix<IndexType, OriginalPointerType> ellMatrix(cooMtx);
    if (!ellMatrix.isValid()) {
        std::cerr << "Matrix " << matrixName << " in ELL format is not valid!\n";
        skip = true;
    }
    
    if (!skip) {
        std::stringstream ssSettingsPageRankEll;
        ssSettingsPageRankEll << std::scientific << "% all PageRank runtimes averaged over " << numberPageRankRunsToAverage
                << " executions (and measured in ns), PageRank executed with damping factor of " << dampingFactor
                << "; Sparse format: ELL";
        
        std::cout << ssSettingsPageRankEll.str() << "\n\n";
        
        std::string pageRankEll("runtime_pagerank_1_ell = [");
        std::stringstream ssPageRank1EllHeader;
        ssPageRank1EllHeader << std::scientific;
        std::stringstream ssPageRank1EllValues;
        
        EllPageRank<1, OriginalPointerType, IndexType> pageEllDouble;
        pageEllDouble.initMatrix(ellMatrix);
        
        ssPageRank1EllHeader << "%";
        ssPageRank1EllValues << "% " << matrixName << '\n';
        for (OriginalPointerType& epsilon : epsList) {
            ssPageRank1EllHeader << headlineSeparator << epsilon;
            auto pageRuntime = runPageRank(pageEllDouble, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
            ssPageRank1EllValues << dataSeparator << pageRuntime[0];
        }
        
        std::cout << pageRankEll << '\n'
                << ssPageRank1EllHeader.str() << '\n'
                << ssPageRank1EllValues.str() << '\n'
                << blocksEnd << '\n';
        pageEllDouble.clear();
        
        
        
        std::string pageRankEll2("runtime_pagerank_2_ell = [");
        std::stringstream ssPageRank2EllHeader;
        ssPageRank2EllHeader << std::scientific;
        std::stringstream ssPageRank2EllValues;
        
        EllPageRank<2, OriginalPointerType, IndexType> pageEll2;
        pageEll2.initMatrix(ellMatrix);
        
        ssPageRank2EllHeader << "%";
        ssPageRank2EllValues << "% " << matrixName << '\n';
        for (OriginalPointerType& epsilon : epsList) {
            ssPageRank2EllHeader << headlineSeparator << epsilon;
            auto pageRuntime = runPageRank(pageEll2, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
            ssPageRank2EllValues << dataSeparator << pageRuntime[0];
        }
        
        std::cout << pageRankEll2 << '\n'
                << ssPageRank2EllHeader.str() << '\n'
                << ssPageRank2EllValues.str() << '\n'
                << blocksEnd << '\n';
        pageEll2.clear();
        
        
        std::string pageRankEll4("runtime_pagerank_4_ell = [");
        std::stringstream ssPageRank4EllHeader;
        ssPageRank4EllHeader << std::scientific;
        std::stringstream ssPageRank4EllValues;
        
        EllPageRank<4, OriginalPointerType, IndexType> pageEll4;
        pageEll4.initMatrix(ellMatrix);
        
        ssPageRank4EllHeader << "%";
        ssPageRank4EllValues << "% " << matrixName << '\n';
        for (OriginalPointerType& epsilon : epsList) {
            ssPageRank4EllHeader << headlineSeparator << epsilon;
            auto pageRuntime = runPageRank(pageEll4, epsilon, dampingFactor, maxIters, numberPageRankRunsToAverage);
            ssPageRank4EllValues << dataSeparator << pageRuntime[0];
        }
        
        std::cout << pageRankEll4 << '\n'
                << ssPageRank4EllHeader.str() << '\n'
                << ssPageRank4EllValues.str() << '\n'
                << blocksEnd << '\n';
        pageEll4.clear();
        
    }
    
    ellMatrix.clear();
    
    return 0;
}

