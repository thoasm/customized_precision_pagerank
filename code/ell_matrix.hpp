#ifndef ELL_MATRIX_HPP
#define ELL_MATRIX_HPP

#include "coo_matrix.hpp"

#include <vector>
#include <type_traits>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <limits>

template<typename IndexType, typename FloatType>
class EllMatrix
{
    static_assert(std::is_floating_point<FloatType>::value, "The second template parameter needs to be a floating point type!");
    
public:
    EllMatrix(const CooMatrix<IndexType, FloatType> &cooMtx);
    
    IndexType getM() const  { return m_; }
    IndexType getN() const  { return n_; }
    IndexType getK() const  { return k_; }
    
    const std::vector<IndexType> &getIndices() const  {   return indices_; }
    const std::vector<FloatType> &getVals() const  {   return vals_; }

    const std::vector<IndexType> &getEmptyColIndices() const  {   return emptyColIndices_; }
    
    
    void clear();
    bool isValid() const    {   return valid_;   }
    ~EllMatrix();
    
private:
    bool valid_;
    IndexType m_;
    IndexType n_;
    IndexType k_;
    std::vector<FloatType> vals_;
    std::vector<IndexType> indices_;
    std::vector<IndexType> emptyColIndices_;
};

#include "ell_matrix.inc.cpp"

#endif
