#ifndef ALIGNMENT_SCORE_HH
#define ALIGNMENT_SCORE_HH

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>

#include "LocARNA/rna_data.hh"

class AlignmentScore : public Gecode::Propagator {
private:

protected:
    
  /// Constructor for cloning \a p
  AlignmentScore(Gecode::Space& home, bool share, AlignmentScore& p);

  /// Constructor for posting \a p
  AlignmentScore(Gecode::Space& home,
		 const RnaData &rna_data_S,
		 const RnaData &rna_data_R,
		 const AlignerParams &params,
		 Gecode::Int::IntViewArray M,
		 Gecode::Int::IntViewArray G,
		 Gecode::Int::SetViewArray H,
		 Gecode::Int::InvVarView Score
		 );


  //! upper bound for the contribution of matching positions i
  //! and j of respective sequences R and S
  score_t
  ub_match(size_type i, size_type j) const;

public:
  //! post a binary neighbor constraint
  static Gecode::ExecStatus post(Gecode::Space& home,
				 const RnaData &rna_data_S,
				 const RnaData &rna_data_R,
				 const AlignerParams &params,
				 Gecode::Int::IntViewArray M,
				 Gecode::Int::IntViewArray G,
				 Gecode::Int::SetViewArray H,
				 Gecode::Int::InvVarView Score
				 );
	    
  /// Copy propagator during cloning
  virtual Gecode::Actor* copy(Gecode::Space& home, bool share);
    
  /// Cost function
  virtual Gecode::PropCost cost(void) const;
    
  /// Perform propagation
  virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta&);
};


#endif // ALIGNMENT_SCORE_HH
