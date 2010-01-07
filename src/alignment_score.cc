#include "alignment_score.hh"
#include "LocARNA/rna_data.hh"

//#include <iostream>


AlignmentScore::AlignmentScore(Gecode::Space& home,
			       const RnaData &rna_data_S_,
			       const RnaData &rna_data_R_,
			       const AlignerParams &params_,
			       const Scoring &scoring_,
			       IntViewArray &M_,
			       IntViewArray &G_,
			       SetViewArray &H_,
			       Gecode::Int::IntView &Score_
			       )
    : Propagator(home),
      rna_data_S(rna_data_S_),
      rna_data_R(rna_data_R_),
      params(params_),
      scoring(scoring_),
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
    params(p.params),
    scoring(p.scoring)
{
    M.update(home,share,p.M);
    G.update(home,share,p.G);
    H.update(home,share,p.H);
    Score.update(home,share,p.Score);
}

Gecode::ExecStatus
AlignmentScore::
post(Gecode::Space& home,
     const RnaData &rna_data_S,
     const RnaData &rna_data_R,
     const AlignerParams &params,
     const Scoring &scoring,
     Gecode::IntVarArray &M,
     Gecode::IntVarArray &G,
     Gecode::SetVarArray &H,
     Gecode::IntVar &Score
     )
{
    if ( false /*simple_fail_test*/ ) {
	return Gecode::ES_FAILED;
    } else {
	// translate vars/var arrays to views/view arrays
      
	Gecode::Int::IntView ScoreView = Score;
      
	Gecode::VarArgArray<Gecode::IntVar> MArg = M;
	IntViewArray MView = IntViewArray(home,MArg);

	Gecode::VarArgArray<Gecode::IntVar> GArg = G;
	IntViewArray GView = IntViewArray(home, GArg);

	Gecode::VarArgArray<Gecode::SetVar> HArg = H;
	SetViewArray HView = SetViewArray(home, HArg);
      
	new (home) AlignmentScore(home,
				  rna_data_S,
				  rna_data_R,
				  params,
				  scoring,
				  MView,GView,HView,
				  ScoreView);
    }
    return Gecode::ES_OK;
}

Gecode::Actor*
AlignmentScore::copy(Gecode::Space& home, bool share) {
    return new (home) AlignmentScore(home,share,*this);
}

Gecode::PropCost
AlignmentScore::cost(const Gecode::Space& home, const Gecode::ModEventDelta& med) const {
    return Gecode::PropCost::linear(Gecode::PropCost::HI,
				    rna_data_S.get_sequence().length()
				    *rna_data_R.get_sequence().length()
				    );
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



Gecode::ExecStatus
AlignmentScore::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {

    // ----------------------------------------
    // Propagation strategy:
    //
    // This is the naive implementation of the propagator!
    // The idea is strictly following "Keep It Simple and Stupid" (KISS-design)
    //
    // propagation is performed in three steps:
    // 1) forward algo 2) backward algo 3) compute edge upper bounds and prune domains
    // in 1) and 2) we fill the respective matrices Fwd and Bwd. The combination
    // of matrix entries Fwd(i-1,j-1)+ub_mach(i,j)+Bwd(i,j) yields the upper bound for match i~j.
    //
    // ----------------------------------------


    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
  
    const Sequence &seqR=rna_data_R.get_sequence();
    const Sequence &seqS=rna_data_S.get_sequence();
  
    const int n=seqR.length();
    const int m=seqS.length();
  

    // ----------------------------------------
    // Forward algorithm
    //
  
    Matrix<infty_score_t> Fwd(n+1,m+1); // ForWarD matrix
  
    Fwd.fill(infty_score_t::neg_infty);
  
    // Fwd(i,j) := max score of a alignment R[1..i], S[1..j]
  
    // initialize
    Fwd(0,0)=(infty_score_t)0;
    for(size_type i=1; i<=n; i++) {
	if (deletion_allowed(i,0)) {
	    Fwd(i,0) = Fwd(i-1,0) + scoring.gapA(i,0);
	}
    }
  
    for(size_type j=1; j<=m; j++) {
	if (insertion_allowed(0,j)) {
	    Fwd(0,j) = Fwd(0,j-1)  + scoring.gapB(0,j);
	}
    }
  
    // recurse
    for(size_type i=1; i<=n; i++) {
	for(size_type j=1; j<=m; i++) {
	    if (match_allowed(i,j)) {
		Fwd(i,j) = Fwd(i-1,j-1) + ub_match(i,j);
	    }
	    if (deletion_allowed(i,j)) {
		Fwd(i,j) = std::max(Fwd(i,j),Fwd(i-1,j)+scoring.gapA(i,j));
	    }
	    if (insertion_allowed(i,j)) {
		Fwd(i,j) = std::max(Fwd(i,j),Fwd(i,j-1)+scoring.gapB(i,j));
	    }
	}
  
	// ----------------------------------------
	// Backward algorithm
	//
  
	Matrix<infty_score_t> Bwd(n+1,m+1); 
	// Bwd(i,j) := max score of a alignment R[i+1..n], S[j+1..m]
  
	// init with "invalid" value negative infinity
	Bwd.fill(infty_score_t::neg_infty);

	// initialize
	Bwd(n,m)=(infty_score_t)0;
	for(size_type i=n; i>0;) { // for i=n-1 downto 0
	    --i;
	    if (deletion_allowed(i+1,m)) {
		Bwd(i,m) = Bwd(i-1,m) + scoring.gapA(i+1,m);
	    }
	}
	for(size_type j=m; j>0;) { // for j=m-1 downto 0
	    --j;
	    if (insertion_allowed(n,j+1)) {
		Bwd(n,j) = Bwd(n,j-1) + scoring.gapB(n,j+1);
	    }
	}
  
	// recurse
	for(size_type i=n; i>0;) { // for i=n-1 downto 0
	    --i;
	    for(size_type j=m; j>0;) { // for j=m-1 downto 0
		--j;
		if (match_allowed(i+1,j+1)) {
		    Bwd(i,j) = Bwd(i+1,j+1) + ub_match(i+1,j+1); // match i+1 ~ j+1 (!)
		}
		if (deletion_allowed(i+1,j+1)) {
		    Bwd(i,j) = std::max(Bwd(i,j), Bwd(i+1,j)+scoring.gapA(i+1,j));
		}
		if (insertion_allowed(i+1,j+1)) {
		    Bwd(i,j) = std::max(Bwd(i,j), Bwd(i,j+1)+scoring.gapB(i,j+1));
		}
	    }
  
	    // ----------------------------------------
	    // Combine matrices and prune
	    //
  
	    for(size_type i=1; i<=n; i++) {
		for(size_type j=1; j<=m; i++) {      
		    // prune M variable for match i~j
		    infty_score_t ubm = Fwd(i-1,j-1)+ub_match(i,j)+Bwd(i,j);
		    if (ubm < Score.min()) {
			ret |= M[i].nq(j);
		    }

		    // prune G variable for gap i~- (deletion of i between j and j+1)
		    infty_score_t ubd = Fwd(i-1,j) + gamma + Bwd(i,j);
		    if (ubd < Score.min()) {
			ret |= G[i].nq(j);
		    }
      
		    // prune H variable for gap -~j (insertion of j between i and i+1)
		    infty_score_t ubi = Fwd(i,j-1) + gamma + Bwd(i,j);
		    if (ubi < Score.min()) {
			ret |= H[i].exclude(j);
		    }
		}
	    }
  
	    GECODE_ME_CHECK(ret);
  
	    // can the propagator be subsumed before all vars are fixed and is it efficiently detectable?
  
	    // test whether all vars fixed, then subsume
	    if ( all_vars_fix ) {
		return ES_SUBSUMED(*this,dispose(home)); 
	    }
  
	    return Gecode::ES_FIX;
	}
