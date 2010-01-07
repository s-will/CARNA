#ifndef ANCHOR_CONSTRAINTS_HH
#define ANCHOR_CONSTRAINTS_HH

#include <string>
#include <vector>
#include <map>

#include <exception>

#include <iostream>


//! maintain the constraints on alignment (sequence) edges that
//! have to be satisfied during the alignment
class AnchorConstraints {
public:
    typedef size_t size_type;
    typedef std::pair<size_type,size_type> size_pair_t;

    class failure : public std::exception {
	std::string msg_;
    public:
	explicit failure (const std::string& msg): msg_(msg) {};
	virtual ~failure() throw();
	virtual const char* what() const throw();
    };
    
private:
    typedef std::map<std::string,size_type> name_tab_t;
    
    typedef std::vector<int> seq_t;

    typedef std::vector<std::string> name_seq_t;

    typedef size_pair_t range_t;
    typedef std::vector<range_t> range_seq_t;

    //! sequence of connected positions in seq B or seq A.
    //! a[i] = j (1<=j<=lenB) means there is an edge from A_i to B_j
    //! a[i] = 0 means there is no name for position i in A
    //! a[i] = -1 means, there is a name for A_i, but no match to B 
    seq_t a;
    
    //! sequence of connected positions in seq A for B.
    //! semantics analogously to a
    seq_t b;

    //! sequence of ranges for a, allows fast constraint checking.
    //! a position i can only be matched to any position in ar[i].first..ar[i].second
    //! without violating an anchor constraint (holds for arbitrary i: 1<=i<=lenA)
    range_seq_t ar; 
  
    //! map from position to names in a
    name_seq_t names_a;
    //! map from position to names in b
    name_seq_t names_b;
    
    //! length of the names
    size_type name_size_;
    
public:
    // ------------------------------------------------------------
    // constructors

    /**
       initialize the data structures
       given two string vectors
       
       The constraints (=alignment edges that have to be satisfied)
       are encoded as follows:
       equal symbols in the sequences for A and B form an edge
       
       In order to specify an arbitrary number of sequences,
       the strings can consist of several lines, then a symbol consists
       of all characters of the column. '.' and ' ' are neutral character
       
       Example:
       seqCA={"..123...."}
       seqCB={"...12.3...."}
       
       specifies the edges (3,4), (4,5), and (5,7)

       Example 2:
       seqCA={"..AAB....",
              "..121...."}
       seqCB={"...AA.B....",
              "...12.1...."}
       specifies the same constraints, allowing a larger name space for constraints.
    */
    
    AnchorConstraints(size_type lenA, 
		   const std::vector<std::string> &seqCA,
		   size_type lenB,
		   const std::vector<std::string> &seqCB);
    
    /**
       initialize the data structures
       given two strings
       
       for semantics see first constructor, where 
       the strings are split into vectors of strings at '#' or '\n'
    */
    AnchorConstraints(size_type lenA,
		   const std::string &seqCA,
		   size_type lenB,
		   const std::string &seqCB);
    
    // -----------------------------------------------------------
    // asking for constraint information
    
    //! is the edge (i,j) allowed?
    //! an edge is allowed, iff it is not in conflict with an anchor constraint
    bool
    allowed_edge(size_type i, size_type j) const {
        return
	  ar[i].first <= j 
	  && j <= ar[i].second;
    }

    //! matching position in b for position <i> in a
    int
    match_to_a(size_type i) const {
	return a[i];
    }
    
    //! matching position in b for position <i> in a
    int
    match_to_b(size_type i) const {
	return b[i];
    }
    
    //! is position i in sequence A aligned to any position in B
    bool
    aligned_in_a(size_type i) const {
	return match_to_a(i)>0;
    }
    
    //! is position j in sequence B aligned to any position in A
    bool
    aligned_in_b(size_type j) const {
	return match_to_b(j)>0;
    }
    
    //! get the name of position i in A
    std::string
    get_name_a(size_type i) const {return names_a[i];}
    
    //! get the name of position j in B
    std::string
    get_name_b(size_type j) const {return names_b[j];}

    //! returns length/size of the names 
    size_type
    name_size() const {return name_size_;};

    //! is the constraint declaration empty
    bool
    empty() const {return name_size()==0;}
    
    //! return the positions (i,j) of the rightmost anchor constraint
    size_pair_t rightmost_anchor() const {
      for (size_type i=a.size(); i>1; ) { // for i=lenA downto 1
	--i;
	if (a[i]>0) return size_pair_t(i,a[i]);
      }
      return size_pair_t(0,0);
    }

    //! return the positions (i,j) of the leftmost anchor constraint 
    size_pair_t leftmost_anchor() const {
      for (size_type i=0; i<=a.size(); i++) {
	if (a[i]>0) return size_pair_t(i,a[i]);
      }
      return size_pair_t(a.size()+1,b.size()+1);
    }

private:
    // ------------------------------------------------------------
    // construction helper
    
    //! Transform the input for the second constructor
    //! to the input for the first.
    //! Basically, splits string seq_str at '#'
    static
    void 
    transform_input(std::vector<std::string> &seqVec,
		    const std::string &seqStr);
    
    //! Translate input vector of strings <seq>
    //! to a map of position names to sequence indices
    //! (also tests input, return true for valid input)
    static
    void
    transform_input(name_tab_t &nameTab,
		    size_type len,
		    const std::vector<std::string> &seq);

    //! Initializes the data structures for efficiently answering
    //! constraint queries later on
    void
    init_tables(const name_tab_t &nameTabA,
		const name_tab_t &nameTabB);

    static
    void
    init_seq_table(seq_t & seq_tab,
		   name_seq_t & name_seq_tab,
		   const name_tab_t &nameTabA,
		   const name_tab_t &nameTabB);

    //! test, whether string/name consists of only don't care symbols
    //! (used for ignoring '.' and ' ' in constraint names)
    static
    bool
    only_dont_care(const std::string &s);
};


#endif // ANCHOR_CONSTRAINTS_HH
