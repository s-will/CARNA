#ifndef ALIGNMENT_SCORE_HH
#define ALIGNMENT_SCORE_HH

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/set.hh>
#include <gecode/search.hh>

#include "LocARNA/rna_data.hh"
#include "LocARNA/scoring.hh"

class AlignmentScore : public Gecode::Propagator {
public:
    typedef Gecode::ViewArray<Gecode::Int::IntView> IntViewArray;
    typedef Gecode::ViewArray<Gecode::Set::SetView> SetViewArray;
    
    typedef size_t size_type;

private:
    
    const Sequence &seqA;
    const Sequence &seqB;
    const ArcMatches &arc_matches;
    const AlignerParams &params;
    const Scoring &scoring;
    IntViewArray M;
    IntViewArray G;
    SetViewArray H;
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
		   IntViewArray &G,
		   SetViewArray &H,
		   Gecode::Int::IntView &Score
		   );


    //! upper bound for the contribution of matching positions i
    //! and j of respective sequences R and S
    infty_score_t
    ub_match(size_type i, size_type j) const;

    //! test whether a match is allowed by the constraint store
    //! @params i position in sequence 1
    //! @params j position in sequence 2
    //! @returns whether the match between i and j is allowed
    bool
    match_allowed(size_type i,size_type j) const {
	return M[i].in((int)j); // is j a member of M[i]
    }

    //! test whether an insertion is allowed by the constraint store
    //! @params i position in sequence 1
    //! @params j position in sequence 2
    //! @returns whether the insertion j between i and i+1 is allowed
    bool
    insertion_allowed(size_type i,size_type j) const {
	return  ! H[i].notContains((int)j); // is j potentially a member of the set G[i]
    }

    //! test whether a deletion is allowed by the constraint store
    //! @params i position in sequence 1
    //! @params j position in sequence 2
    //! @returns whether the deletion of of i between j and j+1 is allowed
    bool
    deletion_allowed(size_type i,size_type j) const {
	return G[i].in((int)j); // is j a member of G[i]
    }

    //! used for testing whether propagator can be deleted
    //! @returns whether all vars are fixed
    bool
    all_vars_fixed() const;


public:
    //! post a binary neighbor constraint
    static Gecode::ExecStatus post(Gecode::Space& home,
				   const Sequence &seqA,
				   const Sequence &seqB,
				   const ArcMatches &arc_matches,
				   const AlignerParams &params,
				   const Scoring &scoring,
				   Gecode::IntVarArray &M,
				   Gecode::IntVarArray &G,
				   Gecode::SetVarArray &H,
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
