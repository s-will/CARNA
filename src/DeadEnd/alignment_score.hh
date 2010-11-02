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
    
    typedef size_t size_type;

private:
    
    const Sequence &seqA;
    const Sequence &seqB;
    const ArcMatches &arc_matches;
    const AlignerParams &params;
    const Scoring &scoring;
    IntViewArray M;
    IntViewArray T;
    Gecode::Int::IntView Score;
    
    const int undef;
    
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
		   IntViewArray &M,
		   IntViewArray &T,
		   Gecode::Int::IntView &Score
		   );

protected:
    
    //! constrain an integer view to be less or undef
    //! @param home Home space
    //! @param x the integer view
    //! @param val value such that x<val or x=undef
    Gecode::ModEvent
    le_or_undef(Gecode::Space &home, Gecode::Int::IntView &xv, int val) {
	Gecode::Iter::Ranges::Singleton r(val,undef-1);
	return xv.minus_r(home,r,false);
    }

    //! @param xv Integer View
    //! @returns maximal domain value of \a x that is not equal undef
    int max_non_undef( const Gecode::Int::IntView &xv ) {
	int m = xv.max();
	if (m!=undef) return m;
	return m-xv.regret_max(); 
    }

    

    //! upper bound for the contribution of matching positions i
    //! and j of respective sequences R and S
    score_t
    ub_match(size_type i, size_type j) const;

    //! calculate the score for an alingment given by trace vectors
    //! @param traceA trace vector for positions in sequence A
    //! @param traceB trace vector for positions in sequence B
    score_t
    evaluate_trace(const std::vector<size_type> &traceA,const std::vector<size_type> &traceB) const;

    score_t
    evaluate_tracematch(const std::vector<size_type> &traceA,const std::vector<size_type> &traceB,
			size_type i,size_type j) const;


    //! test whether a match is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the match between i and j is allowed
    //! assume that all the consistency checking with other variables M,T has been done 
    bool
    match_allowed(size_type i,size_type j) const {
	return M[i].in((int)j) && T[i].in((int)j);
    }

    //! test whether a match is guaranteed (i.e. forced) by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the match between i and j is guaranteed
    //! assume that all the consistency checking with other variables M,T has been done 
    bool
    match_forced(size_type i,size_type j) const {
	return T[i].in((int)j) && M[i].assigned() && (size_t)M[i].val()==j;
    }

    
    //! test whether an insertion is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the insertion j between i and i+1 is allowed
    bool
    insertion_allowed(size_type i,size_type j) const {
	return T[i].in((int)j);
    }
    

    //! test whether a deletion is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the deletion of of i between j and j+1 is allowed
    bool
    deletion_allowed(size_type i,size_type j) const {
	return T[i].in((int)j) && M[i].in((int)undef);
    }

    //! used for testing whether propagator can be deleted
    //! @returns whether all vars are fixed
    bool
    all_vars_fixed() const;


    //! first phase of propagation;
    //! performs bound independent propagation on variables M,T
    Gecode::ModEvent
    simple_consistency(Gecode::Space& home);

    void
    forward_algorithm(Gecode::Space& home, Matrix<infty_score_t> &Fwd);

    void
    backtrace_forward(Gecode::Space &home, const Matrix<infty_score_t> &Fwd,
		      std::vector<size_type> &traceA,
		      std::vector<size_type> &traceB
		      );

    void
    backward_algorithm(Gecode::Space& home, Matrix<infty_score_t> &Bwd);
	
    Gecode::ModEvent 
    prune(Gecode::Space& home, 
			  const Matrix<infty_score_t> &Fwd,
			  const Matrix<infty_score_t> &Bwd);


    //! bound on all arcmatches to the right (left) from i,j
    template<class AdjList>
    score_t
    bound_arcmatches(size_t i, size_t j, AdjList adjlA, AdjList adjlB, bool right) const;
    
    // computes choice for the current space
    void
    choice(RNAalignment &s,
	   const Matrix<infty_score_t> &Fwd,
	   const Matrix<infty_score_t> &Bwd,
	   const std::vector<size_type> &traceA,
	   const std::vector<size_type> &traceB) const;
 
public:
    //! post a binary neighbor constraint
    static Gecode::ExecStatus post(Gecode::Space& home,
				   const Sequence &seqA,
				   const Sequence &seqB,
				   const ArcMatches &arc_matches,
				   const AlignerParams &params,
				   const Scoring &scoring,
				   Gecode::IntVarArray &M,
				   Gecode::IntVarArray &T,
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
