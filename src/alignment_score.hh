#ifndef ALIGNMENT_SCORE_HH
#define ALIGNMENT_SCORE_HH

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/iter.hh>
#include <gecode/set.hh>
#include <gecode/search.hh>

#include "LocARNA/rna_data.hh"
#include "LocARNA/scoring.hh"


class RNAalignment;

class AlignmentScore : public Gecode::Propagator {
public:
    typedef Gecode::ViewArray<Gecode::Int::IntView> IntViewArray;
    // typedef Gecode::ViewArray<Gecode::Set::SetView> SetViewArray;
    typedef Gecode::ViewArray<Gecode::Int::BoolView> BoolViewArray;
    
    typedef size_t size_type;

private:
    
    const Sequence &seqA;
    const Sequence &seqB;
    const ArcMatches &arc_matches;
    const AlignerParams &params;
    const Scoring &scoring;
    IntViewArray MD;
    BoolViewArray M;
    Gecode::Int::IntView Score;
        
protected:
  
    /// Constructor for cloning \a p
    AlignmentScore(Gecode::Space& home, bool share, AlignmentScore& p);

    /// Constructor for posting \a p
    AlignmentScore(Gecode::Space& home,
		   const Sequence &seqA,
		   const Sequence &seqB, 
		   const ArcMatches &arcmatches,
		   const AlignerParams &params,
		   const Scoring &scoring,
		   IntViewArray &MD,
		   BoolViewArray &M,
		   Gecode::Int::IntView &Score
		   );

    virtual ~AlignmentScore();
    
protected:

    void
    print_vars() const;
    
    /*
    //! constrain an integer view to be less or undef
    //! @param home Home space
    //! @param x the integer view
    //! @param val value such that x<val or x=undef
    Gecode::ModEvent
    le_or_undef(Gecode::Space &home, Gecode::Int::IntView &xv, int val) {
	Gecode::Iter::Ranges::Singleton r(val,undef-1);
	return xv.minus_r(home,r,false);
    }

    //! exclude from a set view all values smaller or equal \a val
    //! @param xv set view variable
    //! @param val integer value
    Gecode::ModEvent
    exclude_lq(Gecode::Space &home, Gecode::Set::SetView &xv, int val) {
	Gecode::ModEvent ret=Gecode::ME_GEN_NONE;
	for(int i=xv.lubMin(); i<=val; i++) {
	    ret |= xv.exclude(home,i);
	}
	return ret;
    }

    //! exclude from a set view all values greater or equal \a val
    //! @param xv set view variable
    //! @param val integer value
    Gecode::ModEvent
    exclude_gq(Gecode::Space &home, Gecode::Set::SetView &xv, int val) {
	Gecode::ModEvent ret=Gecode::ME_GEN_NONE;
	for(int i=val; i<=xv.lubMax(); i++) {
	    ret |= xv.exclude(home,i);
	}
	return ret;
    }

    //! @param xv Integer View
    //! @returns maximal domain value of \a x that is not equal undef
    int max_non_undef( const Gecode::Int::IntView &xv ) {
	int m = xv.max();
	if (m!=undef) return m;
	return m-xv.regret_max(); 
    }
    */
    

    //! upper bound for the contribution of matching positions i
    //! and j of respective sequences R and S
    //! if tight!=NULL && tight!=false, return in tight, whether the bound is tight
    score_t
    ub_match(size_type i, size_type j,
	     const Matrix<bool> &considered_ams,
	     const Matrix<score_t> &match_scores,
	     bool *tight=NULL) const;

    //! calculate the score for an alingment given by trace vectors
    //! @param traceA trace vector for positions in sequence A
    //! @param traceB trace vector for positions in sequence B
    score_t
    evaluate_trace(const std::vector<size_type> &traceA,
		   const std::vector<size_type> &traceB,
		   const Matrix<bool> &considered_ams,
		   const Matrix<score_t> &match_scores		   
		   ) const;

    score_t
    evaluate_tracematch(const std::vector<size_type> &traceA,
			const std::vector<size_type> &traceB,
			const Matrix<bool> &considered_ams,
			const Matrix<score_t> &match_scores,
			size_type i,size_type j) const;


    //! test whether a match or deletion is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the match between i and j is of deletion of i between j,j+1 is allowed
    //! assume that all the consistency checking with other variables MD,M has been done 
    bool
    match_or_deletion_allowed(size_type i,size_type j) const {
	return MD[i].in((int)j);
    }


    //! test whether a match is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the match between i and j is allowed
    //! assume that all the consistency checking with other variables MD,M has been done 
    bool
    match_allowed(size_type i,size_type j) const {
	return MD[i].in((int)j) && M[i].in(1);
    }

    //! test whether a match is guaranteed (i.e. forced) by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the match between i and j is guaranteed
    //! assume that all the consistency checking with other variables MD,M has been done. 
    //! match_forced(i,j) implies match_allowed(i,j)
    bool
    match_forced(size_type i,size_type j) const {
	return MD[i].assigned() && (size_t)MD[i].val()==j && !M[i].in(0);
    }
    
    //! test whether an insertion is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the insertion j between i and i+1 is allowed
    bool
    insertion_allowed(size_type i,size_type j) const {
	// insertion has to be inferred from the MD and M vars in this and following line
	// note that a match in line i+1 requires an extra non-insertion position in line i
	const size_t n=seqA.length();
	
	if (i<n) {
	    size_t max=MD[i+1].max();
	    if (!M[i+1].in(0)) { // guaranteed match
		max--;
	    }
	    return MD[i].min()<(int)j && j<=max;
	} else {
	    return MD[i].min()<(int)j;
	}
    }
    
    //! test whether a deletion is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the deletion of of i between j and j+1 is allowed
    bool
    deletion_allowed(size_type i,size_type j) const {
	return MD[i].in((int)j) && M[i].in(0);
    }

    //! used for testing whether propagator can be deleted
    //! @returns whether all vars are fixed
    bool
    all_vars_fixed() const;

    void
    forward_algorithm(Gecode::Space& home,
		      Matrix<infty_score_t> &Fwd,
		      const Matrix<score_t> &UBM);

    void
    backtrace_forward(Gecode::Space &home, 
		      const Matrix<infty_score_t> &Fwd,
		      const Matrix<score_t> &UBM,
		      std::vector<size_type> &traceA,
		      std::vector<size_type> &traceB
		      );

    void
    backward_algorithm(Gecode::Space& home, 
		       Matrix<infty_score_t> &Bwd,
		       const Matrix<score_t> &UBM);
    
    //! set the MD and M variables in the range [start..end] to
    //! the values described by the trace
    //! @param start Begin of range
    //! @param end End of range
    //! @param traceA The trace
    //! @returns Modification event
    //! @pre the range must describe a run
    Gecode::ModEvent
    fix_vars_to_trace(Gecode::Space &home,
		      size_t start,
		      size_t end,
		      const std::vector<size_type> &traceA,
		      const std::vector<size_type> &traceB);

    Gecode::ModEvent
    fix_tight_runs(Gecode::Space &home,
		   const std::vector<size_type> &traceA,
		   const std::vector<size_type> &traceB,
		   const Matrix<bool> &tight);

    Gecode::ModEvent 
    prune(Gecode::Space& home, 
	  const Matrix<infty_score_t> &Fwd,
	  const Matrix<infty_score_t> &Bwd,
	  const Matrix<score_t> &UBM);

    //! determines the indices of arcs that are forced to occur in any
    //! arc match.  return result in the output parameters forcedA and
    //! forcedB
    template<class AdjList>
    void
    determine_forced_arcs(const AdjList &adjlA,
			  const AdjList &adjlB,
			  std::vector<bool> &forcedA, 
			  std::vector<bool> &forcedB,
			  bool right,
			  bool *tight) const;
	
    //! bound on all arcmatches to the right (left) from i,j
    //! if tight!=null and tight!=false return in tight, whether the bound is tight
    template<class AdjList>
    score_t
    bound_arcmatches(const AdjList &adjlA, 
		     const AdjList &adjlB,
		     const Matrix<bool> &considered_ams,
		     bool right,
		     bool *tight) const;
    
    // computes choice for the current space
    void
    choice(RNAalignment &s,
	   const Matrix<infty_score_t> &Fwd,
	   const Matrix<infty_score_t> &Bwd,
	   const std::vector<size_type> &traceA,
	   const std::vector<size_type> &traceB,
	   const Matrix<score_t> &UBM
	   ) const;
    
    void
    prune_decided_arc_matches(Matrix<bool> &considered_ams, Matrix<score_t> &match_scores);

 
public:
    //! post a binary neighbor constraint
    static Gecode::ExecStatus post(Gecode::Space& home,
				   const Sequence &seqA,
				   const Sequence &seqB,
				   const ArcMatches &arc_matches,
				   const AlignerParams &params,
				   const Scoring &scoring,
				   Gecode::IntVarArray &MD,
				   Gecode::BoolVarArray &M,
				   Gecode::IntVar &Score
				   );
	    
    /// Copy propagator during cloning
    virtual Gecode::Actor* copy(Gecode::Space& home, bool share);
    
    /// Cost function
    virtual Gecode::PropCost cost(const Gecode::Space& home, const Gecode::ModEventDelta& med) const;
    
    /// Perform propagation
    virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta&);
};


#endif // ALIGNMENT_SCORE_HH
