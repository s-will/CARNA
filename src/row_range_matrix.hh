#ifndef ROW_RANGE_MATRIX_HH
#define ROW_RANGE_MATRIX_HH

#include <iostream>

// ----------------------------------------
//! matrix class with range for each row.
//
template <class elem_t>
class RowRangeMatrix {
public:
    //! size type
    typedef typename std::vector<elem_t>::size_type size_type;
    
    //! pair of sizes, entry in range vector
    typedef std::pair<size_type,size_type> size_pair_type;

protected:
    
    //! vector of matrix elements
    std::vector<elem_t> mat_; 

#ifndef NDEBUG
    //! stores actual ranges, if we want range checks for debugging
    //! otherwise this information is not necessary after construction 
    std::vector<size_pair_type> ranges;
#endif
    
    std::vector<int> offset; //!< offset for each row
    
    /** 
     * Address of matrix element
     * 
     * @param i row index
     * @param j column index
     * 
     * @return index of element in one dimensional vector
     */
    size_type addr(size_type i, size_type j) const {
	assert(i<ranges.size());
	assert(ranges[i].first<=j && j<=ranges[i].second); 
	return offset[i]+j;
    }
    
public:    
    RowRangeMatrix() {
    }

    /** 
     * construct the matrix from range vector
     * 
     * @param ranges_ vector of column ranges for each matrix row
     */
    RowRangeMatrix(const std::vector<size_pair_type> &ranges_) {
	resize(ranges_);
    }

    /** 
     * resize the matrix according to range vector
     * 
     * @param ranges_ vector of column ranges for each matrix row
     */void
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
    
    /** 
     *  read only access to matrix entry
     * 
     * @param i row index
     * @param j column index
     * 
     * @return constant reference to matrix entry
     */
    const elem_t & operator() (size_type i,size_type j) const {
	return mat_[addr(i,j)];
    }
    
    /** 
     *  read/write access to matrix entry
     * 
     * @param i row index
     * @param j column index
     * 
     * @return reference to matrix entry
     */
    elem_t & operator() (size_type i,size_type j) {
	return mat_[addr(i,j)];
    }

    /** 
     *  read only access to matrix entry
     * 
     * @param i row index
     * @param j column index
     * 
     * @return value of matrix entry
     */
    const elem_t get(size_type i,size_type j) const {
	return mat_[addr(i,j)];
    }

    /** 
     *  write access to matrix entry
     * 
     * @param i row index
     * @param j column index
     * @param x new value
     */
    void
    set(size_type i,size_type j, const elem_t &x) {
	mat_[addr(i,j)]=x;
    }
    
    template<class T>
    friend
    std::ostream &
    operator << (std::ostream &out, const RowRangeMatrix<T> &m);
};

template <class elem_t>
std::ostream &
operator << (std::ostream &out, const RowRangeMatrix<elem_t> &m) {
    #ifndef NDEBUG
    size_t i=0;
    for(typename std::vector<typename RowRangeMatrix<elem_t>::size_pair_type>::const_iterator it=m.ranges.begin();
	m.ranges.end()!=it;
	++it) {
	out << i <<" [ "<< it->first << "-" << it->second << "] : "; 
	for(size_t j=it->first; j<=it->second; ++j) {
	    out <<m.get(i,j)<<" ";
	}
	out << std::endl;
	i++;
    }
    #else
    std::cerr << "RowRangeMatrix::operator << : Turn on debugging."<<std::endl; 
    #endif
    return out;
}

#endif
