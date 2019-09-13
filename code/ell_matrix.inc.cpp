
template<typename IndexType, typename FloatType>
EllMatrix<IndexType, FloatType>::~EllMatrix()
{
    this->clear();
}


template<typename IndexType, typename FloatType>
void EllMatrix<IndexType, FloatType>::clear()
{
    vals_.clear();
    indices_.clear();
    m_ = n_ = k_ = 0;
    valid_ = false;
}

template<typename IndexType, typename FloatType>
EllMatrix<IndexType, FloatType>::EllMatrix(const CooMatrix<IndexType, FloatType> &cooMtx) : m_(cooMtx.getM()), n_(cooMtx.getN())
{
    constexpr FloatType zeroEntry = static_cast<FloatType>(0);
    constexpr IndexType zeroIndex = static_cast<IndexType>(0);
    
    valid_ = false;
    if (!cooMtx.isValid()) {
        return;
    }
    //cooMtx.normalizeColumns();    // Input matrix must already be normalized and Row majored
    //cooMtx.sortRowMajor();
    /* convert the COO matrix to ELL */
    // search for widest row, and therefore k
    IndexType widestRow = 0;
    IndexType currentWidth = 0;
    IndexType countingRow = 0;
    for (IndexType i = 0; i < cooMtx.getNnz(); ++i) {
        const IndexType currentRow = cooMtx.getEntries()[i].i;
        //const IndexType currentColumn = cooEntry[i].j;
        if (countingRow == currentRow) {
            ++currentWidth;
        }
        else {
            if (widestRow < currentWidth) {
                widestRow = currentWidth;
            }
            currentWidth = 0 + 1;
            countingRow = currentRow;
        }
    }
    if (widestRow < currentWidth) {
        widestRow = currentWidth;
    }

    k_ = widestRow;
    if (static_cast<std::size_t>(k_) * m_ > static_cast<std::size_t>(std::numeric_limits<IndexType>::max())) {
        return;
        //throw std::runtime_error("Indexing would exceed the current index type!");
    }
    
    // Write the ELL values
    vals_ = std::vector<FloatType>(m_ * k_, zeroEntry);
    indices_ = std::vector<IndexType>(m_ * k_, zeroIndex);
    IndexType rowValueCount = 0;
    IndexType lastRow = 0;
    for (IndexType i = 0; i < cooMtx.getNnz(); ++i) {
        auto currentEntry = cooMtx.getEntries()[i];
        const IndexType currentRow = currentEntry.i;
        const IndexType currentColumn = currentEntry.j;
        const FloatType currentVal = currentEntry.val;
        
        if (lastRow != currentRow) {
            rowValueCount = 0;
        }
        const IndexType matIdx = rowValueCount * m_ + currentRow;
        vals_[matIdx] = currentVal;
        indices_[matIdx] = currentColumn;
        
        ++rowValueCount;
        lastRow = currentRow;
    }
    
    emptyColIndices_ = cooMtx.getEmptyColumnIndices();
    
    valid_ = true;
}


