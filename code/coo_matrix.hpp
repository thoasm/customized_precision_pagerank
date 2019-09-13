#ifndef COO_MATRIX_HPP
#define COO_MATRIX_HPP

extern "C" {
    #include "mmio.h"
}

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>


template<typename IndexType, typename FloatType>
class CooMatrix
{
public:
    struct CooEntry {
        IndexType i;
        IndexType j;
        FloatType val;
    };
    
    CooMatrix(const char *fileName, bool transposeMtx = false) {
        readMatrix(fileName);
        if (transposeMtx) {
            transpose();
        }
    }
    
    CooMatrix() : valid_(false), m_(0), n_(0), nnz_(0) {}
    
    CooMatrix(const CooMatrix<IndexType, FloatType> &other) : valid_(other.valid_), m_(other.m_), n_(other.n_), nnz_(other.nnz_), entries_(other.entries_), isSymmetric_(other.isSymmetric_) {}
    
    CooMatrix & operator=(const CooMatrix<IndexType, FloatType> &other) {
        valid_ = other.valid_;
        m_ = other.m_;
        n_ = other.n_;
        nnz_ = other.nnz_;
        entries_ = other.entries_;
        isSymmetric_ = other.isSymmetric_;
    }
    
    ~CooMatrix() {
        clear();
    }
    
    bool isValid() const      {   return valid_;    }
    IndexType getM() const    {   return m_;  }
    IndexType getN() const    {   return n_;  }
    IndexType getNnz() const  {   return entries_.size();}
    
    const std::vector<CooEntry> &getEntries() const {
        return entries_;
    }
    
    bool readMatrix(const char *fileName);
    
    bool generateDiagonalMatrix(IndexType n);
    
    void transpose();
    
    void sortColumnMajor();
    void sortRowMajor();

    void normalizeRows();
    void normalizeColumns();
    
    // Problem: After one iteration, empty columns / rows might result after removing the current ones.
    void removeEmptyColumns();
    void removeEmptyRows() {
        transpose();
        removeEmptyColumns();
        transpose();
    }
    
    
    void clear() {
        entries_.clear();
        m_ = 0;
        n_ = 0;
        valid_ = false;
    }
    
    std::vector<IndexType> getEmptyColumnIndices() const;
    std::vector<IndexType> getEmptyRowIndices() const;
    
private:
    bool valid_;
    bool isSymmetric_;
    IndexType m_;
    IndexType n_;
    IndexType nnz_;
    std::vector<CooEntry> entries_;

};


#include "coo_matrix.inc.cpp"


#endif //COO_MATRIX_HPP
