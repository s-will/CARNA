#include "alignment_score.hh"

//#include <iostream>


AlignmentScore::AlignmentScore(Gecode::Space& home,
			       const RNAData &rna_data_S_,
			       const RNAData &rna_data_R_,
			       const AlignerParams &params_,
			       Gecode::Int::IntViewArray &M_,
			       Gecode::Int::IntViewArray &G_,
			       Gecode::Int::SetViewArray &H_,
			       Gecode::Int::InvVarView &Score_
			       )
  : Propagator(home),
    rna_data_S(rna_data_S_),
    rna_data_R(rna_data_R_),
    params(params_),
    M(M_),
    G(G_),
    H(H_),
    Score(Score_)
{
}

AlignmentScore::AlignmentScore(Gecode::Space& home,
			       bool share,
			       AlignmentScore& p) :
  Propagator(home,share,p),
  rna_data_R(p.rna_data_R),
  rna_data_S(p.rna_data_S),
  params(p.params)
{
  M.update(p.M);
  G.update(p.G);
  H.update(p.H);
  Score.update(p.Score);
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
  // by summing over all possible base pair matchs
  //
  // when selecting a certain structure: maximize over possible base pair matchs
  // when not selecting one single structure: make use of the fact that
  // there is at most one alignment edge per base ==> each arc can be matched at most once
  
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ACHTUNG: we need infinity arithmetic (as in LocARNA)
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Gecode::ExecStatus
AlignmentScore::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {

  // ----------------------------------------
  // Propagation strategy:
  //
  // propagation is performed in three steps:
  // 1) forward algo 2) backward algo 3) compute edge upper bounds and prune domains
  // in 1) and 2) we fill the respective matrices Fwd and Bwd. The combination
  // of matrix entries Fwd(i-1,j-1)+ub_mach(i,j)+Bwd(i,j) yields the upper bound for match i~j.
  //
  // ----------------------------------------


  Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
  
  const &seqR=rna_data_R.get_sequence();
  const &seqS=rna_data_S.get_sequence();
  
  const int n=seqR.length();
  const int m=seqS.length();
  
  const score_t gamma = params.gap_extension;

  // ----------------------------------------
  // Forward algorithm
  //
  
  Matrix<score_t> Fwd(n+1,m+1); 
  // Fwd(i,j) := max score of a alignment R[1..i], S[1..j]
  
  // initialize
  Fwd(0,0)=0;
  for(size_type i=1; i<=n; i++) {
    if (deletion_allowed(i,0)) {
      Fwd(i,0) = Fwd(i-1,0) + gamma;
    }
  }
  for(size_type j=1; j<=m; i++) {
    Fwd(0,j) = Fwd(0,j-1) + gamma;
  }
  
  // recurse
  for(size_type i=1; i<=n; i++) {
    for(size_type j=1; j<=m; i++) {
      Fwd(i,j) = neg_infty;
      if (match_allowed(i,j)) {
	Fwd(i,j) = Fwd(i-1,j-1) + ub_match(i,j);
      }
      if (deletion_allowed(i,j)) {
	Fwd(i,j) = std::max(Fwd(i,j),Fwd(i-1,j)+gamma);
      }
      if (insertion_allowed(i,j)) {
	Fwd(i,j) = std::max(Fwd(i,j),Fwd(i,j-1)+gamma);
      }
  }
  
  // ----------------------------------------
  // Backward algorithm
  //
  
  Matrix<score_t> Bwd(n+1,m+1); 
  // Bwd(i,j) := max score of a alignment R[i+1..n], S[j+1..m]
  
  // initialize
  Bwd(n,m)=0;
  for(size_type i=n; i>0;) { // for i=n-1 downto 0
    --i;
    Bwd(i,m) = Bwd(i-1,m) + gamma;
  }
  for(size_type j=m; j>0;) { // for j=m-1 downto 0
    --j;
    Bwd(n,j) = Bwd(n,j-1) + gamma;
  }
  
  // recurse
  for(size_type i=n; i>0;) { // for i=n-1 downto 0
    --i;
    for(size_type j=m; j>0;) { // for j=m-1 downto 0
      --j;
      Bwd(i,j) = Bwd(i+1,j+1) + ub_match(i+1,j+1); // here match i+1 ~ j+1 (!)
      Bwd(i,j) = std::max(Bwd(i,j),
			  std::max(Bwd(i+1,j),
				   Bwd(i,j+1))+gamma);
    }
  }
  
  // ----------------------------------------
  // Combine matrices and prune
  //
  
  for(size_type i=1; i<=n; i++) {
    for(size_type j=1; j<=m; i++) {
      
      // prune M variable for match i~j
      score_t ubm = Fwd(i-1,j-1)+ub_match(i,j)+Bwd(i,j);
      if (ubm < Score.min()) {
	M[i].nq(j);
      }

      // prune G variable for gap i~-
      
    }
  }
  
  
  GECODE_ME_CHECK(ret);
    
  if ( test_for_subsumed ) {	// no further propagation possible
    return Gecode::ES_SUBSUMED(*this, home);
  }
    
  // return Gecode::ES_FIX;
}
