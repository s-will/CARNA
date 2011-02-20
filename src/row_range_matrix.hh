#ifndef ROW_RANGE_MATRIX_HH
#define ROW_RANGE_MATRIX_HH

// ----------------------------------------
//! matrix class with range for each row.
//
template <class elem_t>
class RowRangeMatrix {
public:
    typedef typename std::vector<elem_t>::size_type size_type;
    
    typedef std::pair<size_type,size_type> size_pair_type;

protected:
    std::vector<elem_t> mat_; //!< matrix elements go here

#ifndef NDEBUG
    //! store actual ranges, if we want range checks for debugging
    //! otherwise this information is not necessary after construction 
    std::vector<size_pair_type> ranges;
#endif
    
    std::vector<int> offset; //!< offset for each row
    
    size_type addr(size_type i, size_type j) const {
	assert(i<ranges.size());
	assert(ranges[i].first<=j && j<=ranges[i].second); 
	return offset[i]+j;
    }
    
public:    
    RowRangeMatrix() {
    }

    //! construct the matrix from range vector
    //! @param ranges_ a vector of column ranges for each matrix row
    RowRangeMatrix(const std::vector<size_pair_type> &ranges_) {
	resize(ranges_);
    }

    //! resize the matrix according to range vector
    //! @param ranges_ a vector of column ranges for each matrix row
    void
    resize(const std::vector<size_pair_type> &ranges_) {
	offset.resize(ranges_.size());
	
#ifndef NDEBUG
	// copy ranges_ in debugging mode
	ranges = ranges_;
#endif
	
	// compute offsets from ranges for fast address calculation
	size_type acc=0; // accumulate row widths
	for (size_type i=0; i<ranges_.size(); i++) {
	    offset[i] = acc - (int)ranges_[i].first; // requires integer subtraction!
	    // calc size of row i
	    size_type width=ranges_[i].second-ranges_[i].first+1;
	    // inc accumulator
	    acc += width;
	}
	
	// resize mat_, here acc is the total number of matrix entries
	mat_.resize(acc);
    }
    
    //! read only access to (i,j)
    const elem_t & operator() (size_type i,size_type j) const {
	return mat_[addr(i,j)];
    }
    
    //! read/write access to (i,j)
    elem_t & operator() (size_type i,size_type j) {
	return mat_[addr(i,j)];
    }

    //! read only access to (i,j)
    const elem_t get(size_type i,size_type j) const {
	return mat_[addr(i,j)];
    }

    //! write access to (i,j)
    void
    set(size_type i,size_type j, const elem_t &x) {
	mat_[addr(i,j)]=x;
    }
};

#endif
