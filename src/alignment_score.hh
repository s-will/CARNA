#ifndef ALIGNMENT_SCORE_HH
#define ALIGNMENT_SCORE_HH

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>

#include "LocARNA/rna_data.hh"
#include "LocARNA/scoring.hh"

class AlignerParams;

class AlignmentScore : public Gecode::Propagator {
private:

  const RnaData &rna_data_S;
  const RnaData &rna_data_R;
  const AlignerParams &params;
  IntViewArray &M;
  IntViewArray &G;
  SetViewArray &H;
  Gecode::Int::IntView &Score;

public:
  typedef Gecode::ViewArray<Gecode::Int::IntView> IntViewArray;
  typedef Gecode::ViewArray<Gecode::Set::SetView> SetViewArray;

protected:
    
  
  
  /// Constructor for cloning \a p
  AlignmentScore(Gecode::Space& home, bool share, AlignmentScore& p);

  /// Constructor for posting \a p
  AlignmentScore(Gecode::Space& home,
		 const RnaData &rna_data_S,
		 const RnaData &rna_data_R,
		 const AlignerParams &params,
		 IntViewArray &M,
		 IntViewArray &G,
		 SetViewArray &H,
		 Gecode::Int::IntView &Score
		 );


  //! upper bound for the contribution of matching positions i
  //! and j of respective sequences R and S
  score_t
  ub_match(size_type i, size_type j) const;


  //! test whether a match is allowed by the constraint store
  //! @params i position in sequence 1
  //! @params j position in sequence 2
  //! @returns whether the match between i and j is allowed
  bool
  AlignmentScore::match_allowed(size_type i,size_type j) const {
    return M[i].in(j); // is j a member of M[i]
  }

  //! test whether an insertion is allowed by the constraint store
  //! @params i position in sequence 1
  //! @params j position in sequence 2
  //! @returns whether the insertion j between i and i+1 is allowed
  bool
  AlignmentScore::insertion_allowed(size_type i,size_type j) const {
    return G[i].in(j); // is j a member of G[i]
  }

  //! test whether a deletion is allowed by the constraint store
  //! @params i position in sequence 1
  //! @params j position in sequence 2
  //! @returns whether the deletion of of i between j and j+1 is allowed
  bool
  AlignmentScore::deletion_allowed(size_type i,size_type j) const {
    return H[i].contains(j); // is j potentially a member of the set G[i]
  }

public:
  //! post a binary neighbor constraint
  static Gecode::ExecStatus post(Gecode::Space& home,
				 const RnaData &rna_data_S,
				 const RnaData &rna_data_R,
				 const AlignerParams &params,
				 IntViewArray &M,
				 IntViewArray &G,
				 SetViewArray &H,
				 Gecode::Int::IntView &Score
				 );
	    
  /// Copy propagator during cloning
  virtual Gecode::Actor* copy(Gecode::Space& home, bool share);
    
  /// Cost function
  virtual Gecode::PropCost cost(void) const;
    
  /// Perform propagation
  virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta&);
};


#endif // ALIGNMENT_SCORE_HH
