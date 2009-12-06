#include "AlignmentScore.hh"
#include "../LocARNA/rna_data.hh"

//#include <iostream>


AlignmentScore
::AlignmentScore(Gecode::Space& home,
		 const RNAData &rna_data_S,
		 const RNAData &rna_data_R,
		 const AlignerParams &params,
		 Gecode::Int::IntViewArray M,
		 Gecode::Int::IntViewArray G,
		 Gecode::Int::SetViewArray H,
		 Gecode::Int::InvVarView Score
		 )
  : Propagator(home)
{
  dots;
}

AlignmentScore::AlignmentScore(Gecode::Space& home,
			       bool share,
			       AlignmentScore& p) :
  Propagator(home,share)
{
  dots;
}

Gecode::ExecStatus
AlignmentScore::
post(Gecode::Space& home,
     const RNAData &rna_data_S,
     const RNAData &rna_data_R,
     const AlignerParams &params,
     Gecode::Int::IntViewArray M,
     Gecode::Int::IntViewArray G,
     Gecode::Int::SetViewArray H,
     Gecode::Int::InvVarView Score
     )     
{
  if ( simple_fail_test ) {
    return Gecode::ES_FAILED;
  } else {
    new (home) AlignmentScore(home,
			      rna_data_S,
			      rna_data_R,
			      params,
			      M,G,H,
			      Score);
  }
  return Gecode::ES_OK;
}

Gecode::Actor*
AlignmentScore::copy(Gecode::Space& home, bool share) {
  return new (home) AlignmentScore(home,share,*this);
}

Gecode::PropCost
AlignmentScore::cost(void) const {
  if (x0.assigned() || x1.assigned()) {
    return Gecode::PropCost::unary(Gecode::PropCost::LO);
  }
  return Gecode::PropCost::binary(Gecode::PropCost::LO);
}


score_t
AlignmentScore::ub_match(size_type i, size_type j) const {
  
  // compute upper bound for the contribution of matching positions i
  // and j of respective sequences R and S
  // by ...
  
  
}


Gecode::ExecStatus
AlignmentScore::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {

  Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
  
  // ----------------------------------------
  // Forward algorithm
  //
  Matrix<score_t> Fwd; 
  
  
  
  // ----------------------------------------
  // Backward algorithm
  //
  Matrix<score_t> Bwd; 
  
  
  
  GECODE_ME_CHECK(ret);
    
  if ( test_for_subsumed ) {	// no further propagation possible
    return Gecode::ES_SUBSUMED(*this, home);
  }
    
  // return Gecode::ES_FIX;
}
