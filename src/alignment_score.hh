#ifndef ALIGNMENT_SCORE_HH
#define ALIGNMENT_SCORE_HH


// include config.h
#include "../config.h"

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/iter.hh>
#include <gecode/set.hh>
#include <gecode/search.hh>

#include "LocARNA/params.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/scoring.hh"

#include "row_range_matrix.hh"


// TODO: OPTIMIZE SPEED/SPACE of propagation

// Storing n^2 matrices is not necessary, since most of the time, the
// matrices are very sparse. One can cheaply store only the defined
// values for each matrix row. (introduce new matrix class for this purpose!)
// DONE => use class RowRangeMatrix
// This optimization seems to be very effective, which
// emphasizes the importance of optimizing space and memory management.


// Gecode supports optimized memory management. This is currently not
// used for the dynamic programming matrices, which is likely to have
// a strong negative performance impact.  Currently, the matrices are
// allocated and freed for each call of propagate!  NOTE that we do
// not define them as class memebers, since we don't want to copy the
// matrices with the propagator.  Due to affine gap cost this problem
// becomes even more severe.

// Affine gap cost computation can be space optimized: apparently
// pruning does make use of Fwd/Bwd and FwdA/BwdA but does not require
// FwdB/BwdB.  For this reason, the full matrices Fwd,Bwd,FwdA,BwdA
// have to be stored, but storing FwdB/BwdB is actually not
// necessary. Since we don't need full matrices for computation, we
// could save the space (and allocation/deallocation)!

// Finally, one could combine backward computation and pruning for
// another great speed up!


#include <LocARNA/basepairs.hh>

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
    typedef std::vector<unsigned int> SizeVec;

    typedef LocARNA::Matrix<LocARNA::score_t> ScoreMatrix;
    typedef LocARNA::Matrix<LocARNA::infty_score_t> InftyScoreMatrix;
    typedef RowRangeMatrix<LocARNA::infty_score_t> InftyScoreRRMatrix;
    typedef LocARNA::Matrix<bool> BoolMatrix;

    typedef LocARNA::BasePairs__Arc Arc; //!< arc
    
private:
    
    const LocARNA::Sequence &seqA;
    const LocARNA::Sequence &seqB;
    const LocARNA::ArcMatches &arc_matches;
    const LocARNA::AlignerParams &params;
    const LocARNA::Scoring &scoring;

    /**
     * MD[i] is the position in seqB to that position i in seqA is
     * either matched or after that position i is deleted.  Note that
     * positions of insertions are therefore implicit. The semantics
     * of MD and M is encapsulated by the "allowed" methods.
    */
    IntViewArray MD;
    //! M[i] flags whether position i is deleted after j=MD[i] or matched i~MD[i]
    BoolViewArray M;
    Gecode::Int::IntView Score;
        
protected:
  
    //! Constructor for cloning \a p
    AlignmentScore(Gecode::Space& home, bool share, AlignmentScore& p);

    /**
     * Constructor for posting \a p
     * @param *this home space
     * @param seqA sequence A
     * @param seqB sequence B
     * @param arcmatches locarna object defining arc matches
     * @param aligner_params locarna object defining parameters for the alignment
     * @param scoring locarna object defining the scoring scheme
     * @param MD constraint variables defining "match or deletion" position in seqB for each position of seqA
     * @param M constraint variables defining whether MD refers to match (or deletion)
     * @param Score constraint variable for alignment score
    */
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

    /** 
     * \brief virtual destructor
     * 
     */
    virtual ~AlignmentScore();
    
protected:
    
    /** 
     * \brief print constraint variables MD,M,Score for debugging
     * 
     * @param ostream output stream
     */
    void
    print_vars(std::ostream &out) const;

    /** 
     * \brief print trace (DEBUGGING)
     * 
     * @param out output stream 
     * @param traceA trace vector a
     * @param traceB trace vector b
     * @param trace_score score of the trace
     */
    void
    print_trace(std::ostream &out, const SizeVec &traceA,
		const SizeVec &traceB, LocARNA::score_t trace_score) const;
    
    /**
     * upper bound for the contribution of a match 
     *
     * @param i index in first sequence R
     * @param j index in second sequence S
     * @param considered_ams 
     * @param match_scores 
     * @param[in,out] tight return in tight, whether the bound is
     * tight (if tight!=NULL && tight!=false)
     * 
     * @note the bound is tight iff no arc matches contribute to the bound,
     * i.e. iff the bound is 0
     *
     * @return upper bound for contribution of matching i and j 
     */
    LocARNA::score_t
    ub_match(size_type i, size_type j,
	     const ScoreMatrix &considered_ams,
	     const ScoreMatrix &match_scores,
	     bool *tight=NULL) const;

    /**
     * calculate the score for an alingment given by trace vectors
     * @param traceA trace vector for positions in sequence A
     * @param traceB trace vector for positions in sequence B
    */
    LocARNA::score_t
    evaluate_trace(const SizeVec &traceA,
		   const SizeVec &traceB,
		   const ScoreMatrix &considered_ams,
		   const ScoreMatrix &match_scores		   
		   ) const;

    /** 
     * evaluate one match in the trace
     * 
     * @param traceA 
     * @param traceB 
     * @param considered_ams 
     * @param match_scores 
     * @param i 
     * @param j 
     * 
     * @return score contribution for matching i and j in trace 
     */
    LocARNA::score_t
    evaluate_tracematch(const SizeVec &traceA,
			const SizeVec &traceB,
			const ScoreMatrix &considered_ams,
			const ScoreMatrix &match_scores,
			size_type i,size_type j) const;
    
    /**
     * test whether a match is guaranteed (i.e. forced) by the constraint store
     * @param i position in sequence 1
     * @param j position in sequence 2
     * @returns whether the match between i and j is guaranteed
     * assume that all the consistency checking with other variables MD,M has been done. 
     * match_forced(i,j) implies match_allowed(i,j)
    */
    bool
    match_forced(size_type i,size_type j) const {
	return MD[i].assigned() && (size_t)MD[i].val()==j && !M[i].in(0);
    }

    /**
     * test whether a match or deletion is allowed by the constraint store
     * @param i position in sequence 1
     * @param j position in sequence 2
     * @returns whether the match between i and j is of deletion of i between j,j+1 is allowed
     * assume that all the consistency checking with other variables MD,M has been done 
    */
    bool
    match_or_deletion_allowed(size_type i,size_type j) const {
	return MD[i].in((int)j);
    }

    /**
     * minimum column in matrix row.
     * @params i matrix row
     * @returns the minimal j where (i,j) is potentially on a trace
     * through the alignment matrices due to the local MD, M
     * variables.
    */
    size_t min_col(size_t i) const {
	assert(i<=(size_t)MD.size());
	return MD[i].min();
    }
    
    /**
     * maximum column in matrix row.
     * @params i matrix row
     * @returns the maximal j where (i,j) is potentially on a trace
     * through the alignment matrices due to the local MD, M
     * variables.
    */
    size_t max_col(size_t i) const {
	assert(i<=(size_t)MD.size());
	
	if (i==(size_t)MD.size()-1) {
	    return seqB.length();
	}
	
	int max_c = MD[i+1].max();
	
	if (!M[i+1].in(0)) --max_c;
	
	return max_c;
    }

    /** 
     * test whether trace through a matrix entry is allowed
     * 
     * @param i index in first sequence R
     * @param j index in second sequence S
     * 
     * @return whether trace through (i,j) is allowed 
     */
    bool trace_allowed(size_t i, size_t j) const {
	return min_col(i)<=j && j<=max_col(i);
    }
    
    // test whether match/insertion/deletion at (i,j) is allowed due
    // to the variables M and MD. If allowed, the methods guarantee
    // that origin and target of the trace arrow are valid, i.e. min_col(i)<=j<=max_col(i).
    
    /**
     * test whether a match is allowed by the constraint store
     * @param i position in sequence 1
     * @param j position in sequence 2
     * @returns whether the match between i and j is allowed
    */
    bool
    match_arrow_allowed(size_type i,size_type j) const {
	assert(1<=i && i<=seqA.length());
	assert(1<=j && j<=seqB.length());
	
	return MD[i-1].min()<(int)j && MD[i].in((int)j) && M[i].in(1);
    }

    /**
     * test whether an insertion is allowed by the constraint store
     * @param i position in sequence 1
     * @param j position in sequence 2
     * @returns whether the insertion j between i and i+1 is allowed
    */
    bool
    insertion_arrow_allowed(size_type i,size_type j) const {
	assert(0<=i && i<=seqA.length());
	assert(1<=j && j<=seqB.length());
	return min_col(i)<j && j<=max_col(i);
    }
    
    /**
     * test whether a deletion is allowed by the constraint store
     * @param i position in sequence 1
     * @param j position in sequence 2
     * @returns whether the deletion of i between j and j+1 is allowed
    */
    bool
    deletion_arrow_allowed(size_type i,size_type j) const {
	assert(1<=i && i<=seqA.length());
	assert(0<=j && j<=seqB.length());
	
	return min_col(i-1)<=j && MD[i].in((int)j) && M[i].in(0);
    }
    
    /**
     * used for testing whether propagator can be deleted
     * @returns whether all vars are fixed
    */
    bool
    all_vars_fixed() const;
    
    
    /**
     * run forward algorithm for computing prefix alignment scores
     * in Fwd, given matrix of upper bounds for matches. Support affine gap cost model.
     * @param home home space
     * @param Fwd matrix of size n+1 x m+1, output parameter
     * @param FwdA matrix of size n+1 x m+1, output parameter
     * @param FwdB matrix of size n+1 x m+1, output parameter
     * @param UBM precomputed upper bounds for all matches i~j
     *
     * Result is returned in Fwd, FwdA, and FwdB. Fwd(i,j) is the
     * best prefix alignment score for prefixes A_1..i and
     * B_1..j. FwdA(i,j) same as Fwd(i,j) with restriction that
     * A_i is deleted. FwdB(i,j) same as Fwd(i,j)
     * with restriction that B_j is inserted.
    */
    void
    forward_algorithm(Gecode::Space& home,
			     InftyScoreRRMatrix &Fwd,
			     InftyScoreRRMatrix &FwdA,
			     InftyScoreRRMatrix &FwdB,
			     const ScoreMatrix &UBM);

    /** 
     * Perform forward algorithm
     * 
     * @param home 
     * @param Fwd 
     * @param FwdA 
     * @param FwdB 
     * @param UBM precomputed upper bounds for all matches i~j
     * @param traceA 
     * @param traceB 
     */
    void
    backtrace_forward(Gecode::Space &home, 
		      const InftyScoreRRMatrix &Fwd,
		      const InftyScoreRRMatrix &FwdA,
		      const InftyScoreRRMatrix &FwdB,
		      const ScoreMatrix &UBM,
		      SizeVec &traceA,
		      SizeVec &traceB
		      );

    /** 
     * Perform backward algorithm
     * 
     * @param home 
     * @param Bwd 
     * @param BwdA 
     * @param BwdB 
     * @param UBM precomputed upper bounds for all matches i~j 
     */
    void
    backward_algorithm(Gecode::Space& home, 
			      InftyScoreRRMatrix &Bwd,
			      InftyScoreRRMatrix &BwdA,
			      InftyScoreRRMatrix &BwdB,
			      const ScoreMatrix &UBM);
    
    /**
     * set the MD and M variables in the range [start..end] to
     * the values described by the trace
     * @param start Begin of range
     * @param end End of range
     * @param traceA The trace
     * @returns Modification event
     * @pre the range must describe a run
    */
    Gecode::ModEvent
    fix_vars_to_trace(Gecode::Space &home,
		      size_t start,
		      size_t end,
		      const SizeVec &traceA,
		      const SizeVec &traceB);

    /**
     * Determine the runs in the trace that are 'tight' and therefore
     * have to occur in any solution as in the trace. Then, fix the
     * constraint variables accordingly.
     * 
     * @param home the home space
     * @param traceA the trace vector for A
     * @param traceB the trace vector for B
     * @param tight matrix of precomputed tightness of upper bounds
     * for all matches i~j
     * 
     * @return gecode mod event 
     */
    Gecode::ModEvent
    fix_tight_runs(Gecode::Space &home,
		   const SizeVec &traceA,
		   const SizeVec &traceB,
		   const BoolMatrix &tight);


    /**
     * prune the domains of the constraint variables M and MD
     * according to the bounds derivable from forward and backward
     * matrices.
     * 
     * @param home 
     * @param Fwd 
     * @param FwdA 
     * @param FwdB 
     * @param Bwd 
     * @param BwdA 
     * @param BwdB 
     * @param UBM precomputed upper bounds for all matches i~j
     * 
     * @return 
     */
    Gecode::ModEvent 
    prune(Gecode::Space& home, 
	  const InftyScoreRRMatrix &Fwd,
	  const InftyScoreRRMatrix &FwdA,
	  //const InftyScoreRRMatrix &FwdB,
	  const InftyScoreRRMatrix &Bwd,
	  const InftyScoreRRMatrix &BwdA,
	  //const InftyScoreRRMatrix &BwdB,
	  const ScoreMatrix &UBM);

    /**
     * determines the indices of arcs that are forced to occur in any
     * arc match.  return result in the output parameters forcedA and
     * forcedB
     * 
     * @param adjlA 
     * @param adjlB 
     * @param forcedA 
     * @param forcedB 
     * @param right 
     * @param tight 
     */
    template<class AdjList>
    void
    determine_forced_arcs(const AdjList &adjlA,
			  const AdjList &adjlB,
			  BoolVec &forcedA, 
			  BoolVec &forcedB,
			  bool right,
			  bool *tight) const;
	
    /** 
     * bound on all arcmatches to the right (left) from i,j
     * if tight!=null and tight!=false return in tight, whether the bound is tight
     * 
     * @param adjlA adjacency list of a position in A 
     * @param adjlB adjacency list of a position in B
     * @param considered_ams the matrix of (scores of) considered arc matches
     * @param right whether the lists adjlA and adjlB are right adjacency lists
     * 
     * @note the adjacency lists have to be either both right or both
     * left adjacency lists; as usual in the locarna lib, the lists
     * have to be sorted
     *
     * @return score upper bound of contributions due to the arc matches 
     *
     * @note the upper bound takes into account that the other ends of
     * simultaneously scored arc matches can not cross each
     * other. (The method determines the maximum clique of compatible
     * arc matches by DP.)
     *
     * @note arc matches with negative scores do never contribute
     */
    template<class AdjList>
    LocARNA::score_t
    bound_arcmatches(const AdjList &adjlA, 
		     const AdjList &adjlB,
		     const ScoreMatrix &considered_ams,
		     bool right) const;
    
    /** 
     * \brief compute choice for the current space
     * 
     * @param s the space
     * @param Fwd the forward matrix
     * @param Bwd the backward matrix
     * @param traceA trace vector A
     * @param traceB trace vector B
     * @param trace_score score of the trace
     * @param UBM precomputed upper bounds for all matches i~j
     * @param match_scores matrix of match scores
     *
     * Compute the choice in the propagator, such that
     * the heuristic can be guided by the bounds without the
     * need for costly re-computation in the brancher.
     */
    void
    choice(RNAalignment &s,
	   const InftyScoreRRMatrix &Fwd,
	   const InftyScoreRRMatrix &Bwd,
	   const SizeVec &traceA,
	   const SizeVec &traceB,
	   LocARNA::score_t trace_score,
	   const ScoreMatrix &UBM,
	   const ScoreMatrix &match_scores
	   ) const;

    /** 
     * \brief remove already decided arc matches from further consideration
     *
     * Determine arc matches where the match of the left or right ends
     * is determined and remove them from further consideration. In
     * such cases, the similarity of the arc match can be transfered
     * to the match similarity of the undetermined ends
     * 
     * @note this method initializes matrices of match scores and
     * considered arc matches, the match score matrix is modified to
     * additionally reflect the contributions of decided arc matches.
     * These arc matches are then removed from the "list"
     * considered_ams.
     *
     * @param[out] considered_ams matrix of scores for arc matches
     * @param[out] match_scores   matrix of corrected scores for single matchs i~j
     *
     * @note the matrix considered_ams contains entries (a,b) for
     * matching the arcs with indices a and b. If the arc is still
     * "considered", the entry is the score of the arcmatch and 0
     * otherwise. "considered" means that the score contribution
     * is not reflected in match_scores.
     */
    void
    prune_decided_arc_matches(ScoreMatrix &considered_ams, ScoreMatrix &match_scores);

public:
    
    /**
     * post the alignment score constraint
     *
     * @param home home space
     * @param seqA sequence A
     * @param seqB sequence B
     * @param arcmatches locarna object defining arc matches
     * @param aligner_params locarna object defining parameters for the alignment
     * @param scoring locarna object defining the scoring scheme
     * @param MD constraint variables defining "match or deletion" position in seqB for each position of seqA
     * @param M constraint variables defining whether MD refers to match (or deletion)
     * @param Score constraint variable for alignment score
    */
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
    
    /** 
     * Copy propagator during cloning
     * 
     * @param home home space 
     * @param share Gecode share flag
     * 
     * @return pointer to actor
     */
    virtual Gecode::Actor* copy(Gecode::Space& home, bool share);
    
    /** 
     * Cost function
     * 
     * @param home home space 
     * @param med delta for mod event
     * 
     * @return cost of propagation
     */
    virtual Gecode::PropCost cost(const Gecode::Space& home, const Gecode::ModEventDelta& med) const;
    
    /** 
     * Perform propagation
     * 
     * @param home the home space 
     * 
     * @return Gecode execution status
     */
    virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta&);
};


#endif // ALIGNMENT_SCORE_HH
