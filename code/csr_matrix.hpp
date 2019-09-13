#ifndef CSR_MATRIX_HPP
#define CSR_MATRIX_HPP

#include "coo_matrix.hpp"

#include <vector>
#include <type_traits>
#include <algorithm>
#include <limits>
#include <stdexcept>

template<typename IndexType, typename FloatType>
class CsrMatrix
{
    static_assert(std::is_floating_point<FloatType>::value, "The second template parameter needs to be a floating point type!");
    
public:
    CsrMatrix(const CooMatrix<IndexType, FloatType> &cooMtx);
    
    IndexType getM() const  { return m_; }
    IndexType getN() const  { return n_; }
    IndexType getNnz() const  { return nnz_; }
    
    const std::vector<IndexType> &getColIndices() const  {   return colIndices_; }
    const std::vector<IndexType> &getRowStartIndices() const  {   return rowStartIndices_; }
    const std::vector<FloatType> &getVals() const  {   return vals_; }

    const std::vector<IndexType> &getEmptyColIndices() const  {   return emptyColIndices_; }
    
    
    void clear() {
        valid_ = false;
        m_ = n_ = nnz_ = 0;
        vals_.clear();
        rowStartIndices_.clear();
        emptyColIndices_.clear();
    }
    
    bool isValid() const    {   return valid_;   }
    ~CsrMatrix() {
        clear();
    }
    
    //static constexpr FloatType zeroEntry = static_cast<FloatType>(0);
    static constexpr IndexType zeroIndex = static_cast<IndexType>(0);

private:
    bool valid_;
    IndexType m_;
    IndexType n_;
    IndexType nnz_;
    std::vector<FloatType> vals_;
    std::vector<IndexType> colIndices_;
    std::vector<IndexType> rowStartIndices_;
    std::vector<IndexType> emptyColIndices_;
};

#include "csr_matrix.inc.cpp"

#endif
