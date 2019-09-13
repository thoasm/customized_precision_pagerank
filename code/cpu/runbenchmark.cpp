#include "split_pointer.hpp"
#include "ell_benchmark.hpp"
#include "coo_matrix.hpp"
#include "ell_matrix.hpp"
#include "thread_buffer.hpp"


#include <iostream>
#include <cstdint>
#include <string>
#include <type_traits>
#include <iomanip>
#include <algorithm>






bool headerPrinted = false;

template<MySegmentsNumberType SegmentsToUse, MySegmentsNumberType NumberSegments, typename OriginalPointerType, typename IndexType>
void evaluateBenchmark(EllBenchmark<NumberSegments, OriginalPointerType, IndexType> &bench)
{
    const std::string separator = ";  ";
    const std::string newLine = "\n";
    
    const std::string nameTitle = "precision";
    const std::string rwTestTitle = "read-write-kernel [ns]";
    const std::string spmvTitle = "SpMV-kernel [ns]";
    const std::string ownNormTitle = "norm-kernel [ns]";
    
    const std::size_t minWidth = 10;
    const std::size_t nameWidth = std::max(nameTitle.length(), std::size_t{15});
    
    const std::size_t ownNormWidth = std::max(ownNormTitle.length(), minWidth);
    const std::size_t rwTestWidth = std::max(rwTestTitle.length(), minWidth);
    const std::size_t spmvWidth = std::max(spmvTitle.length(), minWidth);
    
    std::cout << std::setfill(' ') << std::right;
    
    if (!headerPrinted) {
        /*
        std::cout << nameTitle << ": " << nameTitle.length() << ", " << nameWidth << newLine
                 << ownNormTitle << ": " << ownNormTitle.length() << ", " << ownNormWidth << newLine
                 << rwTestTitle << ": " << rwTestTitle.length() << ", " <<  rwTestWidth << newLine
                 << spmvTitle << ": " << spmvTitle.length() << ", " <<  spmvWidth << newLine
                 << cublasNormTitle << ": " << cublasNormTitle.length() << ", " << cublasNormWidth << newLine;
        std::cout << std::setw(nameWidth) << nameWidth
                << separator << std::setw(rwTestWidth) << rwTestWidth
                << separator << std::setw(spmvWidth) << spmvWidth
                << separator << std::setw(ownNormWidth) << ownNormWidth
                << separator << std::setw(cublasNormWidth) << cublasNormWidth
                << separator << newLine;
        */
        std::cout << std::setw(nameWidth) << nameTitle
                << separator << std::setw(rwTestWidth) << rwTestTitle
                << separator << std::setw(spmvWidth) << spmvTitle
                << separator << std::setw(ownNormWidth) << ownNormTitle
                << newLine;
        headerPrinted = true;
    }
    
    
    std::string name;
    if (std::is_same<OriginalPointerType, double>::value == true) {
        name = "double";
    }
    else if (std::is_same<OriginalPointerType, float>::value == true) {
        name = "float";
    }
    if (NumberSegments > 1) {
        name += std::string(" (") + std::to_string(SegmentsToUse) + " / " + std::to_string(NumberSegments) + ")";
    }
    
    std::cout << std::setw(nameWidth) << name
            << separator << std::setw(rwTestWidth) << bench.template runRWTest<SegmentsToUse>()
            << separator << std::setw(spmvWidth) << bench.template runSpmv<SegmentsToUse>()
            << separator << std::setw(ownNormWidth) << bench.template runOwnNorm<SegmentsToUse>()
            ;
    std::cout << newLine;
}

int main(int argc, char *argv[])
{
    using IndexType = unsigned int;
    using OriginalPointerType = double;
    
    
    if (argc < 2) {
        std::cout << "You need to specify a sparse matrix!\n";
        return 1;
    }
    CooMatrix<IndexType, OriginalPointerType> cooMtx;
    if (std::string(argv[1]) == std::string("--diagonal")) {
        if (argc < 3) {
            std::cerr << "You need to specify the size of the matrix!\n";
            return 1;
        }
        cooMtx.generateDiagonalMatrix(std::stoi(std::string(argv[2])));
    } else {
        cooMtx.readMatrix(argv[1]);
    }
    
    if (!cooMtx.isValid())
        return 1;
    
    std::cout << "% Optimizations: " << KERNEL_ENABLED_MACROS << std::endl;
    
    ThreadBuffer::getInstance().reserveBufferSizePerThread(4 * memoryBlockSize);
    
    cooMtx.transpose();
    cooMtx.normalizeColumns();
    cooMtx.sortRowMajor();
    EllMatrix<IndexType, OriginalPointerType> ellMtx(cooMtx);
    
    std::cout << memoryBlockSize << " Byte segment splitting; Using thread-buffer; \n";
    std::cout << "ellMtxrix: " << ellMtx.getM() << " x " << ellMtx.getN() << " with k = " << ellMtx.getK() << '\n';
    std::cout << "Number empty columns: " << ellMtx.getEmptyColIndices().size()
            << "\nNumber of non-zeros:  " << cooMtx.getEntries().size() << '\n';
    std::cout << "\n\n";
    
    EllBenchmark<1, double, IndexType> bench1;
    EllBenchmark<2, double, IndexType> bench2;
    EllBenchmark<4, double, IndexType> bench4;
    
    std::cout << "Decl Completed!\n";
    bench1.setUpMatrix(ellMtx);
    std::cout << "SetUp Completed!\n\n";
    
    evaluateBenchmark<1>(bench1);
    
    bench1.clear();
    bench2.setUpMatrix(ellMtx);
    bench2.testConversion();
    evaluateBenchmark<1>(bench2);
    evaluateBenchmark<2>(bench2);
    
    bench2.clear();
    bench4.setUpMatrix(ellMtx);
    bench4.testConversion();
    evaluateBenchmark<1>(bench4);
    evaluateBenchmark<4>(bench4);
    
    return 0;
}
