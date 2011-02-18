#ifndef ALIGNMENT_SCORE_HH
#define ALIGNMENT_SCORE_HH


// include config.h
#include "../config.h"

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/iter.hh>
#include <gecode/set.hh>
#include <gecode/search.hh>

#include "LocARNA/rna_data.hh"
#include "LocARNA/scoring.hh"



class RNAalignment;


/**
   Propagator that relates a pairwise aligmnent described by variables
   MD and M to its sequence-structure alignment score.
 */
class AlignmentScore : public Gecode::Propagator {
public:
    typedef Gecode::ViewArray<Gecode::Int::IntView> IntViewArray;
    // typedef Gecode::ViewArray<Gecode::Set::SetView> SetViewArray;
    typedef Gecode::ViewArray<Gecode::Int::BoolView> BoolViewArray;
    
    typedef size_t size_type;

    typedef std::vector<bool> BoolVec;
    typedef std::vector<size_type> SizeVec;

    typedef LocARNA::Matrix<LocARNA::score_t> ScoreMatrix;
    typedef LocARNA::Matrix<LocARNA::infty_score_t> InftyScoreMatrix;
    typedef LocARNA::Matrix<bool> BoolMatrix;
    
private:
    
    const LocARNA::Sequence &seqA;
    const LocARNA::Sequence &seqB;
    const LocARNA::ArcMatches &arc_matches;
    const LocARNA::AlignerParams &params;
    const LocARNA::Scoring &scoring;

    //! MD[i] is the position in seqB to that position i in seqA is
    //! either matched or after that position i is deleted.  Note that
    //! positions of insertions are therefore implicit. The semantics
    //! of MD and M is encapsulated by the "allowed" methods.
    IntViewArray MD;
    //! M[i] flags whether position i is deleted after j=MD[i] or matched i~MD[i]
    BoolViewArray M;
    Gecode::Int::IntView Score;
        
protected:
  
    /// Constructor for cloning \a p
    AlignmentScore(Gecode::Space& home, bool share, AlignmentScore& p);

    /// Constructor for posting \a p
    AlignmentScore(Gecode::Space& home,
		   const LocARNA::Sequence &seqA,
		   const LocARNA::Sequence &seqB, 
		   const LocARNA::ArcMatches &arcmatches,
		   const LocARNA::AlignerParams &params,
		   const LocARNA::Scoring &scoring,
		   IntViewArray &MD,
		   BoolViewArray &M,
		   Gecode::Int::IntView &Score
		   );

    virtual ~AlignmentScore();
    
protected:

    void
    print_vars() const;
    
    //! upper bound for the contribution of matching positions i
    //! and j of respective sequences R and S
    //! if tight!=NULL && tight!=false, return in tight, whether the bound is tight
    LocARNA::score_t
    ub_match(size_type i, size_type j,
	     const ScoreMatrix &considered_ams,
	     const ScoreMatrix &match_scores,
	     bool *tight=NULL) const;

    //! calculate the score for an alingment given by trace vectors
    //! @param traceA trace vector for positions in sequence A
    //! @param traceB trace vector for positions in sequence B
    LocARNA::score_t
    evaluate_trace(const SizeVec &traceA,
		   const SizeVec &traceB,
		   const ScoreMatrix &considered_ams,
		   const ScoreMatrix &match_scores		   
		   ) const;

    LocARNA::score_t
    evaluate_tracematch(const SizeVec &traceA,
			const SizeVec &traceB,
			const ScoreMatrix &considered_ams,
			const ScoreMatrix &match_scores,
			size_type i,size_type j) const;


    
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

    //! test whether a match or deletion is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the match between i and j is of deletion of i between j,j+1 is allowed
    //! assume that all the consistency checking with other variables MD,M has been done 
    bool
    match_or_deletion_allowed(size_type i,size_type j) const {
	return MD[i].in((int)j);
    }

    //! minimum column in matrix row.
    //! @params i matrix row
    //! @returns the minimal j where (i,j) is potentially on a trace
    //! through the alignment matrices due to the local MD, M
    //! variables.
    size_t min_col(size_t i) const {
	assert(i<=(size_t)MD.size());
	return MD[i].min();
    }
    
    //! maximum column in matrix row.
    //! @params i matrix row
    //! @returns the maximal j where (i,j) is potentially on a trace
    //! through the alignment matrices due to the local MD, M
    //! variables.
    size_t max_col(size_t i) const {
	assert(i<=(size_t)MD.size());
	
	if (i==(size_t)MD.size()-1) {
	    return seqB.length();
	}
	
	int max_c = MD[i+1].max();
	
	if (!M[i+1].in(0)) --max_c;
	
	return max_c;
    }
    
    // test whether match/insertion/deletion at (i,j) is allowed due
    // to the variables M and MD. If allowed, the methods guarantee
    // that origin and target of the trace arrow are valid, i.e. min_col(i)<=j<=max_col(i).
    
    //! test whether a match is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the match between i and j is allowed
    bool
    match_arrow_allowed(size_type i,size_type j) const {
	assert(1<=i && i<=seqA.length());
	assert(1<=j && j<=seqB.length());
	
	return MD[i-1].min()<(int)j && MD[i].in((int)j) && M[i].in(1);
    }

    //! test whether an insertion is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the insertion j between i and i+1 is allowed
    bool
    insertion_arrow_allowed(size_type i,size_type j) const {
	assert(0<=i && i<=seqA.length());
	assert(1<=j && j<=seqB.length());
	return min_col(i)<j && j<=max_col(i);
    }
    
    //! test whether a deletion is allowed by the constraint store
    //! @param i position in sequence 1
    //! @param j position in sequence 2
    //! @returns whether the deletion of i between j and j+1 is allowed
    bool
    deletion_arrow_allowed(size_type i,size_type j) const {
	assert(1<=i && i<=seqA.length());
	assert(0<=j && j<=seqB.length());
	
	return min_col(i-1)<=j && MD[i].in((int)j) && M[i].in(0);
    }
    
    //! used for testing whether propagator can be deleted
    //! @returns whether all vars are fixed
    bool
    all_vars_fixed() const;
    
    
    //! run forward algorithm for computing prefix alignment scores
    //! in Fwd, given matrix of upper bounds for matches. Support affine gap cost model.
    //! @param home home space
    //! @param Fwd matrix of size n+1 x m+1, output parameter
    //! @param FwdA matrix of size n+1 x m+1, output parameter
    //! @param FwdB matrix of size n+1 x m+1, output parameter
    //! @param precomputed UBM upper bounds for all matches i~j
    //!
    //! Result is returned in Fwd, FwdA, and FwdB. Fwd(i,j) is the
    //! best prefix alignment score for prefixes A_1..i and
    //! B_1..j. FwdA(i,j) same as Fwd(i,j) with restriction that
    //! A_i is deleted. FwdB(i,j) same as Fwd(i,j)
    //! with restriction that B_j is inserted.
    void
    forward_algorithm(Gecode::Space& home,
			     InftyScoreMatrix &Fwd,
			     InftyScoreMatrix &FwdA,
			     InftyScoreMatrix &FwdB,
			     const ScoreMatrix &UBM);

    void
    backtrace_forward(Gecode::Space &home, 
		      const InftyScoreMatrix &Fwd,
		      const InftyScoreMatrix &FwdA,
		      const InftyScoreMatrix &FwdB,
		      const ScoreMatrix &UBM,
		      SizeVec &traceA,
		      SizeVec &traceB
		      );

    void
    backward_algorithm(Gecode::Space& home, 
			      InftyScoreMatrix &Bwd,
			      InftyScoreMatrix &BwdA,
			      InftyScoreMatrix &BwdB,
			      const ScoreMatrix &UBM);
    
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
		      const SizeVec &traceA,
		      const SizeVec &traceB);

    Gecode::ModEvent
    fix_tight_runs(Gecode::Space &home,
		   const SizeVec &traceA,
		   const SizeVec &traceB,
		   const BoolMatrix &tight);


    Gecode::ModEvent 
    prune(Gecode::Space& home, 
	  const InftyScoreMatrix &Fwd,
	  const InftyScoreMatrix &FwdA,
	  //const InftyScoreMatrix &FwdB,
	  const InftyScoreMatrix &Bwd,
	  const InftyScoreMatrix &BwdA,
	  //const InftyScoreMatrix &BwdB,
	  const ScoreMatrix &UBM);

    //! determines the indices of arcs that are forced to occur in any
    //! arc match.  return result in the output parameters forcedA and
    //! forcedB
    template<class AdjList>
    void
    determine_forced_arcs(const AdjList &adjlA,
			  const AdjList &adjlB,
			  BoolVec &forcedA, 
			  BoolVec &forcedB,
			  bool right,
			  bool *tight) const;
	
    //! bound on all arcmatches to the right (left) from i,j
    //! if tight!=null and tight!=false return in tight, whether the bound is tight
    template<class AdjList>
    LocARNA::score_t
    bound_arcmatches(const AdjList &adjlA, 
		     const AdjList &adjlB,
		     const ScoreMatrix &considered_ams,
		     bool right,
		     bool *tight) const;
    
    // computes choice for the current space
    void
    choice(RNAalignment &s,
	   const InftyScoreMatrix &Fwd,
	   const InftyScoreMatrix &Bwd,
	   const SizeVec &traceA,
	   const SizeVec &traceB,
	   const ScoreMatrix &UBM,
	   const ScoreMatrix &match_scores
	   ) const;
    
    void
    prune_decided_arc_matches(ScoreMatrix &considered_ams, ScoreMatrix &match_scores);

 
public:
    //! post constraint
    static Gecode::ExecStatus post(Gecode::Space& home,
				   const LocARNA::Sequence &seqA,
				   const LocARNA::Sequence &seqB,
				   const LocARNA::ArcMatches &arc_matches,
				   const LocARNA::AlignerParams &params,
				   const LocARNA::Scoring &scoring,
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
