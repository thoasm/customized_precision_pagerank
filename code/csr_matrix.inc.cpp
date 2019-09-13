
template<typename IndexType, typename FloatType>
CsrMatrix<IndexType, FloatType>::CsrMatrix(const CooMatrix<IndexType, FloatType> &cooMtx) : m_(cooMtx.getM()), n_(cooMtx.getN()), nnz_(cooMtx.getNnz())
{
    valid_ = false;
    if (!cooMtx.isValid()) {
        return;
    }
    //cooMtx.normalizeColumns();    // Input matrix must already be normalized and row majored
    //cooMtx.sortRowMajor();
    /* convert the COO matrix to CSR */
    // search for widest row, and therefore k
    vals_.resize(nnz_);
    colIndices_.resize(nnz_);
    rowStartIndices_.resize(m_ + 1);
    
    
    // Write the CSR values
    IndexType lastRow{0};
    rowStartIndices_[0] = 0;
    for (IndexType i = 0; i < cooMtx.getNnz(); ++i) {
        auto currentEntry = cooMtx.getEntries()[i];
        const IndexType currentRow = currentEntry.i;
        const IndexType currentColumn = currentEntry.j;
        const FloatType currentVal = currentEntry.val;
        
        if (lastRow != currentRow) {
            for (IndexType row = lastRow + 1; row <= currentRow; ++row) {
                rowStartIndices_[row] = i;
            }
        }
        vals_[i] = currentVal;
        colIndices_[i] = currentColumn;
        
        lastRow = currentRow;
    }
    
    for (IndexType row = lastRow + 1; row <= m_; ++row) {
        rowStartIndices_[row] = nnz_;
    }
    
    emptyColIndices_ = cooMtx.getEmptyColumnIndices();
    
    valid_ = true;
}


