/*
template<typename IndexType, typename FloatType>
CooMatrix<IndexType, FloatType>::
*/
template<typename IndexType, typename FloatType>
bool CooMatrix<IndexType, FloatType>::readMatrix(const char *fileName)
{
    valid_ = false;
    clear();
    /**************************************************************************/
    /* We first read in the matrix from a .mtx file.                          */
    /**************************************************************************/
    
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    IndexType i;
    int M, N, nz;
    
    if ((f = fopen(fileName, "r")) == NULL) {
        printf("Unable to open file '%s'\n", fileName);
        return false;
    }
    
    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        return false;
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (! ((mm_is_real(matcode) || mm_is_coordinate(matcode)) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode)) ) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        return false;
    }

    /* find out size of sparse matrix .... */

    ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz);
    this->isSymmetric_ = mm_is_symmetric(matcode);
    
    this->nnz_ = (this->isSymmetric_) ? 2 * nz : nz;
    
    this->m_ = M;
    this->n_ = N;
    
    if (ret_code !=0) {
        printf("Problem with the ret_code (%d)!\n", ret_code);
        return false;
    }
    //printf("%d x %d, %d nzs\n", M, N, nz);
    
    const CooEntry emptyCoo{static_cast<IndexType>(0), static_cast<IndexType>(0), static_cast<FloatType>(0)};
    
    entries_ = std::vector<CooEntry>(this->nnz_, emptyCoo);
    
    if (mm_is_coordinate(matcode)) {
        for (i=0; i<this->nnz_; ++i) {
            int fscRet = fscanf(f, "%d %d\n", &(entries_[i].i), &(entries_[i].j));
            if (fscRet < 2 || fscRet == EOF) {
                return false;
            }
            --entries_[i].i;  /* adjust from 1-based to 0-based */
            --entries_[i].j;
            entries_[i].val = static_cast<FloatType>(1);
            if (this->isSymmetric_) {
                entries_[i+1] = CooEntry{entries_[i].j, entries_[i].i, entries_[i].val};
                ++i;
            }
        }
    }
    else /* (mm_is_real(matcode)) */ {
        for (i=0; i<this->nnz_; ++i) {
            int fscRet = fscanf(f, "%d %d %lg\n", &(entries_[i].i), &(entries_[i].j), &(entries_[i].val));
            if (fscRet < 2 || fscRet == EOF) {
                return false;
            }
            --entries_[i].i;  /* adjust from 1-based to 0-based */
            --entries_[i].j;
            if (this->isSymmetric_) {
                entries_[i+1] = CooEntry{entries_[i].j, entries_[i].i, entries_[i].val};
                ++i;
            }
        }
    }
    
    if (f != stdin) fclose(f);
    
    valid_ = true;
    return true;
}


template<typename IndexType, typename FloatType>
bool CooMatrix<IndexType, FloatType>::generateDiagonalMatrix(IndexType n)
{
    clear();
    m_ = n_ = n;
    entries_ = std::vector<CooEntry>(n);
    for (IndexType i = 0; i < n; ++i) {
        entries_[i].i = entries_[i].j = i;
        entries_[i].val = static_cast<FloatType>(1);
    }
    valid_ = true;
    return valid_;
}


template<typename IndexType, typename FloatType>
void CooMatrix<IndexType, FloatType>::transpose()
{
    if (!valid_ || isSymmetric_) {
        return;
    }
    for (std::size_t i = 0; i < entries_.size(); ++i) {
        const IndexType oldI = entries_[i].i;
        entries_[i].i = entries_[i].j;
        entries_[i].j = oldI;
    }
}

template<typename IndexType, typename FloatType>
void CooMatrix<IndexType, FloatType>::sortColumnMajor()
{
    auto compareFuncColMajor = [this](const CooEntry &left, const CooEntry &right) -> bool {
        return left.j < right.j || (left.j == right.j && left.i < right.i);
    };
    std::sort(entries_.begin(), entries_.end(), compareFuncColMajor);
}

template<typename IndexType, typename FloatType>
void CooMatrix<IndexType, FloatType>::sortRowMajor()
{
    auto compareFuncRowMajor = [this](const CooEntry &left, const CooEntry &right) -> bool {
        return left.i < right.i || (left.i == right.i && left.j < right.j);
    };
    std::sort(entries_.begin(), entries_.end(), compareFuncRowMajor);
}

template<typename IndexType, typename FloatType>
void CooMatrix<IndexType, FloatType>::normalizeColumns()
{
    sortColumnMajor();
    
    IndexType currentColumn = 0;
    IndexType firstElementColumn = 0;
    FloatType columnSum = static_cast<FloatType>(0.0);
    
    for (IndexType i = 0; i < entries_.size(); ++i) {
        if (entries_[i].j == currentColumn) {
            columnSum += entries_[i].val;
        }
        else {
            for (IndexType normIdx = firstElementColumn; normIdx < i; ++normIdx) {
                entries_[normIdx].val /= columnSum;
            }
            currentColumn = entries_[i].j;
            columnSum = entries_[i].val;
            firstElementColumn = i;
        }
    }
    // do the same for the last Column
    for (IndexType normIdx = firstElementColumn; normIdx < entries_.size(); ++normIdx) {
        entries_[normIdx].val /= columnSum;
    }
}


template<typename IndexType, typename FloatType>
void CooMatrix<IndexType, FloatType>::normalizeRows()
{
    sortRowMajor();
    
    IndexType currentRow = 0;
    IndexType firstElementRow = 0;
    FloatType rowSum = static_cast<FloatType>(0.0);
    
    for (IndexType i = 0; i < entries_.size(); ++i) {
        if (entries_[i].i == currentRow) {
            rowSum += entries_[i].val;
        }
        else {
            for (IndexType normIdx = firstElementRow; normIdx < i; ++normIdx) {
                entries_[normIdx].val /= rowSum;
            }
            currentRow = entries_[i].i;
            rowSum = entries_[i].val;
            firstElementRow = i;
        }
    }
    // do the same for the last Column
    for (IndexType normIdx = firstElementRow; normIdx < entries_.size(); ++normIdx) {
        entries_[normIdx].val /= rowSum;
    }
}


// Problem: After one iteration, empty columns might result after removing the current ones.
template<typename IndexType, typename FloatType>
void CooMatrix<IndexType, FloatType>::removeEmptyColumns()
{
    isSymmetric_ = false;   //might destroy symmetry
    
    sortColumnMajor();
    
    IndexType lastRow{0};
    IndexType firstElementRow{0};
    FloatType rowSum = static_cast<FloatType>(0.0);
    
    std::vector<IndexType> emptyColumns;
    for (IndexType i = 0; i < entries_.size(); ++i) {
        if (i == 0 && entries_[i].j > 0) {
           for (std::size_t rowsToAdd = 0; rowsToAdd < entries_[i].j; ++rowsToAdd) {
                emptyColumns.push_back(rowsToAdd);
            }
        }
        else if (entries_[i].j > lastRow + 1) {
            const IndexType colsRemoved = emptyColumns.size();
            for (std::size_t colsToAdd = lastRow + 1; colsToAdd < entries_[i].j; ++colsToAdd) {
                emptyColumns.push_back(colsToAdd);
            }
        }
        lastRow = entries_[i].j;
        entries_[i].j -= emptyColumns.size();
    }
    // Add the remaining until the last row
    for (std::size_t colsToAdd = (entries_.size() > 0) ? lastRow + 1 : 0; colsToAdd < n_; ++colsToAdd) {
        emptyColumns.push_back(colsToAdd);
    }
    
    sortRowMajor();
    IndexType startingPoint{0};
    IndexType removedEntries{0};
    for (IndexType emptyColNumber = 0; emptyColNumber < emptyColumns.size(); ++emptyColNumber) {
        const IndexType idxToDelete = emptyColumns[emptyColNumber];
        
        for (IndexType curMatIdx = startingPoint; curMatIdx < entries_.size(); ++curMatIdx) {
            const IndexType curRow = entries_[curMatIdx].i;
            if (curRow == idxToDelete) {
                entries_[curMatIdx].i = m_;
                entries_[curMatIdx].j = n_;
                ++removedEntries;
            } else if (curRow > idxToDelete) {
                startingPoint = curMatIdx;
                break;
            } else {
                entries_[curMatIdx].i -= emptyColNumber;
            }
        }
    }
    for (IndexType curMatIdx = startingPoint; curMatIdx < entries_.size(); ++curMatIdx) {
        entries_[curMatIdx].i -= emptyColumns.size();
    }

    
    sortRowMajor();
    IndexType actualRemovedEntries{0};
    IndexType nnz = entries_.size();
    for (IndexType i = 0; i < nnz; ++i) {
        const auto lastElement = entries_.back();
        if (lastElement.i == m_ && lastElement.j == n_) {
            entries_.pop_back();
            ++actualRemovedEntries;
        }
    }
    if (actualRemovedEntries != removedEntries || nnz != entries_.size() + removedEntries)
        throw "BIG MISTAKE!!!";
    
    m_ -= emptyColumns.size();
    n_ -= emptyColumns.size();
    nnz_ = entries_.size();
}



template<typename IndexType, typename FloatType>
std::vector<IndexType> CooMatrix<IndexType, FloatType>::getEmptyColumnIndices() const
{
    std::vector<bool> isColEmpty(n_, true);
    for (std::size_t i = 0; i < entries_.size(); ++i) {
        isColEmpty[entries_[i].j] = false;
    }
    std::size_t numberEmptyCols = 0;
    for (std::size_t i = 0; i < n_; ++i) {
        if (isColEmpty[i]) {
            ++numberEmptyCols;
        }
    }
    std::vector<IndexType> emptyColumnIndices(numberEmptyCols);
    
    std::size_t nextEmptySlot = 0;
    for (std::size_t i = 0; i < n_; ++i) {
        if (isColEmpty[i]) {
            emptyColumnIndices[nextEmptySlot++] = i;
        }
    }
    std::sort(emptyColumnIndices.begin(), emptyColumnIndices.end());
    return emptyColumnIndices;
}


template<typename IndexType, typename FloatType>
std::vector<IndexType> CooMatrix<IndexType, FloatType>::getEmptyRowIndices() const
{
    std::vector<bool> isRowEmpty(n_, true);
    for (std::size_t i = 0; i < entries_.size(); ++i) {
        isRowEmpty[entries_[i].i] = false;
    }
    std::size_t numberEmptyRows = 0;
    for (std::size_t i = 0; i < n_; ++i) {
        if (isRowEmpty[i]) {
            ++numberEmptyRows;
        }
    }
    std::vector<IndexType> emptyRowIndices(numberEmptyRows);
    
    std::size_t nextEmptySlot = 0;
    for (std::size_t i = 0; i < n_; ++i) {
        if (isRowEmpty[i]) {
            emptyRowIndices[nextEmptySlot++] = i;
        }
    }
    std::sort(emptyRowIndices.begin(), emptyRowIndices.end());
    return emptyRowIndices;
}
