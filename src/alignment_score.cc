#include "alignment_score.hh"
#include "LocARNA/arc_matches.hh"

#include "WinDisplay.cpp"
#include "RNAalignment.hh"

//#include <iostream>

// DETAIL: for the upper bound in ub_match, arc match scores to the
// left and to the right are added. In total, thus each arcmatch score is added twice.
// Dividing the score by two may cause rounding problems. Thus, better multiply all other
// scores by two (not done yet).
// NOTE: it does not work to add only scores either to the left or 
// to the right (depending on whether doing forward or backward computation!)


AlignmentScore::AlignmentScore(Gecode::Space& home,
			       const Sequence &seqA_,
			       const Sequence &seqB_, 
			       const ArcMatches &arc_matches_,
			       const AlignerParams &params_,
			       const Scoring &scoring_,
			       IntViewArray &M_,
			       IntViewArray &G_,
			       SetViewArray &H_,
			       Gecode::Int::IntView &Score_
			       )
    : Propagator(home),
      seqA(seqA_),
      seqB(seqB_),
      arc_matches(arc_matches_),
      params(params_),
      scoring(scoring_),
      M(M_),
      G(G_),
      H(H_),
      Score(Score_),
      undef(seqB.length()+1)
{
    // subscribe the propagator to any variable change 
    M.subscribe(home,*this,Gecode::Int::PC_INT_DOM);
    G.subscribe(home,*this,Gecode::Int::PC_INT_DOM);
    H.subscribe(home,*this,Gecode::Set::PC_SET_ANY);
    Score.subscribe(home,*this,Gecode::Int::PC_INT_BND);    
}

AlignmentScore::AlignmentScore(Gecode::Space& home,
			       bool share,
			       AlignmentScore& p) :
    Propagator(home,share,p),
    seqA(p.seqA),
    seqB(p.seqB),
    arc_matches(p.arc_matches),
    params(p.params),
    scoring(p.scoring),
    undef(p.undef)
{
    M.update(home,share,p.M);
    G.update(home,share,p.G);
    H.update(home,share,p.H);
    Score.update(home,share,p.Score);
}

Gecode::ExecStatus
AlignmentScore::
post(Gecode::Space& home,
     const Sequence &seqA,
     const Sequence &seqB,
     const ArcMatches &arc_matches,
     const AlignerParams &params,
     const Scoring &scoring,
     Gecode::IntVarArray &M,
     Gecode::IntVarArray &G,
     Gecode::SetVarArray &H,
     Gecode::IntVar &Score
     )
{

    // the post method converts Vars to Views and constructs a new propagator
    // in the home-space

    // translate vars/var arrays to views/view arrays
	
    Gecode::Int::IntView ScoreView = Score;
      
    Gecode::VarArgArray<Gecode::IntVar> MArg = M;
    IntViewArray MView = IntViewArray(home,MArg);
    
    Gecode::VarArgArray<Gecode::IntVar> GArg = G;
    IntViewArray GView = IntViewArray(home, GArg);
    
    Gecode::VarArgArray<Gecode::SetVar> HArg = H;
    SetViewArray HView = SetViewArray(home, HArg);
    
    new (home) AlignmentScore(home,
			      seqA,
			      seqB,
			      arc_matches,
			      params,
			      scoring,
			      MView,GView,HView,
			      ScoreView);
    
    return Gecode::ES_OK;
}

Gecode::Actor*
AlignmentScore::copy(Gecode::Space& home, bool share) {
    return new (home) AlignmentScore(home,share,*this);
}

Gecode::PropCost
AlignmentScore::cost(const Gecode::Space& home, const Gecode::ModEventDelta& med) const {
    return Gecode::PropCost::linear(Gecode::PropCost::LO,
				    (unsigned int) (seqA.length() * seqB.length())
				    );
}


score_t
AlignmentScore::ub_match(size_type i, size_type j) const {
    //std::cout << "AlignmentScore::ub_match"<<std::endl;
  
    // compute upper bound for the contribution of matching positions i
    // and j of respective sequences A and B
    // by summing over all possible base pair matchs
    //
    // when selecting a certain structure: maximize over possible base pair matchs
    // when not selecting one single structure: make use of the fact that
    // there is at most one alignment edge per base ==> each arc can be matched at most once
    
    score_t bound=2*scoring.basematch(i,j);
    
    // NOTE: there are different possibilities to enumerate the possible arc-matches.
    // traversing the arc-matches is good for theoretical complexity,
    // since we assume this is limited by a constant
    
    
    // for all pairs of arcs in A and B that have right ends i and j, respectively
    //
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
	arc_matches.common_right_end_list(i,j).end() != it; ++it ) {
	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	if ( match_allowed(am.arcA().left(),am.arcB().left()) ) { 
	    // if the left ends of arcs arcA and arcB can match 
	    bound += scoring.arcmatch(am);
	}
    }
    
    // same for common left ends
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(i,j).begin();
	arc_matches.common_left_end_list(i,j).end() != it; ++it ) {	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	if ( match_allowed(am.arcA().right(),am.arcB().right()) ) { 
	    // if the right ends of arcs arcA and arcB can match 
	    bound += scoring.arcmatch(am);
	}
    }
    
    return (bound/2);
}


bool
AlignmentScore::all_vars_fixed() const {
    // test all variable views and return false as soon as one is not fixed/assigned
    // otherwise return true
    
    //M
    for (size_type i=1; i<M.size(); ++i) {
	if (!M[i].assigned()) return false;
    }

    //G
    for (size_type i=1; i<G.size(); ++i) {
	if (!G[i].assigned()) return false;
    }
    
    //H
    for (size_type i=0; i<H.size(); ++i) {
	if (!H[i].assigned()) return false;
    }

    return true;
}


score_t
AlignmentScore::evaluate_trace(std::vector<size_type> &traceA,
			       std::vector<size_type> &traceB) const {
    
    const int n=seqA.length();
    const int m=seqB.length();

    score_t score=0;

    size_type i=1;
    size_type j=1;
    
    while (i<n && j<m) {
	if (traceA[i]==0) {
	    score += scoring.gapA(i,j);
	    ++i;
	} else if (traceA[j]==0) {
	    score += scoring.gapB(i,j);
	    ++j;
	} else {
	    // match between traceB[j]==i and traceA[i]==j
	    
	    score_t matchscore=scoring.basematch(i,j);
	    
	    // for all pairs of arcs in A and B that have right ends i and j, respectively
	    //
	    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
		arc_matches.common_right_end_list(i,j).end() != it; ++it ) {
		
		const ArcMatch &am = arc_matches.arcmatch(*it);
		
		if ( traceA[am.arcA().left()] == am.arcB().left() ) { 
		    // if the left ends of arcs arcA and arcB match 
		    matchscore += scoring.arcmatch(am);
		}
	    }
	    
	    // same for common left ends
	    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(i,j).begin();
		arc_matches.common_left_end_list(i,j).end() != it; ++it ) {	
		const ArcMatch &am = arc_matches.arcmatch(*it);
		
		if ( traceA[am.arcA().right()] == am.arcB().right() ) { 
		    // if the right ends of arcs arcA and arcB  match 
		    matchscore += scoring.arcmatch(am);
		}
	    }
	    
	    score += matchscore/2;

	    ++i;
	    ++j;
	}
    }
    
    return score;
}


Gecode::ExecStatus
AlignmentScore::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
    
    //
    // TODO: clean up code by having subroutines for the single steps of the propagation, i.e.
    // consistency checking, forward, backward, combination and pruning
    //
    
    std::cout << "AlignmentScore::propagate " <<std::endl;
    
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
    
    const int n=seqA.length();
    const int m=seqB.length();
    
    std::cout << "Before consistency"<<std::endl;
    std::cout << "Matches:    " << M << std::endl;
    std::cout << "Deletions:  " << G << std::endl;
    std::cout << "Insertions: " << H << std::endl;
    std::cout << "Score:      " << Score << std::endl;

    // ----------------------------------------
    // consistency checking all variables
    // such that during forward and backward algo,
    // we can do simple checks for allowed match, deletion and insertion
    // and such that inconsistent values are removes early
    // 
    
    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
    
    // propagate in two sweeps    

    // handling of insertions is still wrong
    
    // fails for Matches:    {10, 3, 4, {4,10}, {4,10}, 5, 8, 9, 10, 10}

    int min_j=0;
    int min_j2=-1;
    for (int i=1; i<M.size(); i++) {
	// invariant:
	// min_j is the minimal j that is consumed left of i,
	// i.e. either it is matched, deleted after j or inserted
	//
	// min_j2 is the minimal j that is matched, such that G[i] has to be strictly larger
	
	ret |= M[i].gr(home,min_j);
	ret |= G[i].gr(home,min_j-1);
	ret |= G[i].gr(home,min_j2);
	ret |= exclude_lq(home,H[i],min_j);

	if (!M[i].in(undef)) {
	    ret |= G[i].eq(home,undef);
	    if (M[i].assigned()) {
		ret |= H[i].exclude(home,M[i].val());
	    }
	    min_j = std::max( min_j, M[i].min() );
	    min_j2 = M[i].min();
	}
	
	if (!G[i].in(undef)) {
	    ret |= M[i].eq(home,undef);
	    min_j = std::max( min_j, G[i].min() );	
	}
	
	if (H[i].glbSize()!=0) { // if empty set is not possible value of H[i]
	    min_j = std::max( min_j, H[i].glbMin() );
	}
    }

    
    int max_j=m+1;
    for (int i=n; i>1;) {
	--i;
	
	ret |= le_or_undef(home,M[i],max_j);
	ret |= le_or_undef(home,G[i],max_j);
	ret |= exclude_gq(home,H[i],max_j);
	
	if (!M[i].in(undef)) {
	    ret |= G[i].eq(home,undef);
	    if (M[i].assigned()) {
		ret |= H[i].exclude(home,M[i].val());
	    }
	    max_j = std::min( max_j, M[i].max());
	}
	
	if (!G[i].in(undef)) {
	    ret |= M[i].eq(home,undef);
	    max_j = std::min( max_j, G[i].max() + 1 );
	}
	
	if (H[i].glbSize()!=0) { // empty set is not possible value of H[i]
	    max_j = std::min( max_j, H[i].glbMax() );
	}
	
    }

    ret |= exclude_gq(home,H[0],max_j);
    
    // each i is either matched or deleted
    for (int i=1; i<=n; i++) {
	if (M[i].assigned() && M[i].val()==undef) {
	    ret |= G[i].nq(home,undef);
	}
	if (G[i].assigned() && G[i].val()==undef) {
	    ret |= M[i].nq(home,undef);
	}
    }
    
    // each j is either matched or deleted
     
    // determine values j that are not consumed right of i
    int maxGLB[n+1]; // maximal value that has to be included in H[i] at a matched position
    maxGLB[n]=m;
    for (int i=n+1; i>0;) {
	--i;
	if ( !M[i].in(undef) && G[i].in(undef) ) { // in case we have a match
	    maxGLB[i-1] = M[i].min()-1; // what's certainly in for smaller i
        } else if ( M[i].in(undef) && !G[i].in(undef))  { // case where we have a gap
	    maxGLB[i-1] = maxGLB[i]; // carry over deletions
	} else { // when we don't know
	    maxGLB[i-1]=0;
	}
    }
     
    for (int i=0; i<=n; i++) {
	if ( !M[i].in(undef) && G[i].in(undef) ) {
	    for (int k=M[i].min()+1; k<=maxGLB[i]; k++) {
		ret |= H[i].include(home,k);
	    }
	}
	if ( !G[i].in(undef) ) {
	    ret |= H[i].cardMax(home,0); // symmetry breaking
	}
    }
    
    for (int k=1; k<=maxGLB[0]; k++) {
	ret |= H[0].include(home,k);
    }
    
    GECODE_ME_CHECK(ret);

    // missing deletions and insertions and cross compatibility with matches
    
    // DONE consistency check

    ret=Gecode::ME_GEN_NONE;

    std::cout << "Made consistent:"<<std::endl;
    std::cout << "Matches:    " << M << std::endl;
    std::cout << "Deletions:  " << G << std::endl;
    std::cout << "Insertions: " << H << std::endl;
    std::cout << "Score:      " << Score << std::endl;

    
    
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
	for(size_type j=1; j<=m; j++) {
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
    }
    
    std::cout << "Fwd"<<std::endl<<Fwd <<std::endl;
    
    // the forward algorithm yields an upper bound
    // ==> constrain the score variable
    if (Fwd(n,m).is_finite()) {
	GECODE_ME_CHECK(Score.lq(home,(int)Fwd(n,m).finite_value()));
    } else {
	std::cout << "infinite upper bound"<<std::endl;
	return Gecode::ES_FAILED;
    }

    // ----------------------------------------
    // do backtracking of the Fwd matrix in order to obtain a lower bound of alignment score
    //

    //
    // OBSERVATION:
    // doing the trace back helps in the very beginning for having a lower bound at all
    // later trace back usually won't yield better lower bounds
    // i.e. in general its probably best to use some cheapo heuristic for getting
    // a first good alignment and starting with some lower bound. Traceback in the propagator
    // can consequently be omitted.
    //


    // trace vectors. idea one entry per position gives position in the other sequence to that
    // the pos is matched or 0 (in case of gap).
    std::vector< size_type > traceA;
    std::vector< size_type > traceB;
    
    traceA.resize(n+1);
    traceB.resize(m+1);
    
    {
	int i=n;
	int j=m;
	
	while (i>0 && j>0) {
	    if ( Fwd(i,j) == Fwd(i-1,j-1) + ub_match(i,j) ) {
		traceA[i]=j;
		traceB[j]=i;
		--i; 
		--j;
	    } else if ( Fwd(i,j) == Fwd(i-1,j) + scoring.gapA(i,j) ) {
		traceA[i]=0;
		--i;
	    } else if ( Fwd(i,j) == Fwd(i,j-1) + scoring.gapB(i,j) ) {
		traceB[j]=0;
		--j;
	    }
	}
	
	while (i>0) {
	    traceA[i]=0;
	    --i;
	}
	while (j>0) {
	    traceB[j]=0;
	    --j;
	}
    }
    
    // print trace (DEBUGGING)
    std::cout << "TRACE"<<std::endl;
    for (size_type i=1; i<=n; ++i ) {
	std::cout << traceA[i] << " ";
    }
    std::cout << std::endl;
    for (size_type j=1; j<=m; j++) {
	std::cout << traceB[j] << " ";
    }
    std::cout << std::endl;
    

    // copy trace vectors to the space of the propagator
    static_cast<RNAalignment&>(home).traceA=traceA;
    static_cast<RNAalignment&>(home).traceB=traceB;

    // evaluate alignment for obtaining lower bound
    score_t lower_bound = evaluate_trace(traceA,traceB);
    
    std::cout << "Score: " << lower_bound << std::endl;

    
    // constrain Score by lower bound
    Score.gq(home,(int)lower_bound);
    
    // end of traceback and lower bounding
    // ----------------------------------------
    
    
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
	    Bwd(i,m) = Bwd(i+1,m) + scoring.gapA(i+1,m);
	}
    }
    for(size_type j=m; j>0;) { // for j=m-1 downto 0
	--j;
	if (insertion_allowed(n,j+1)) {
	    Bwd(n,j) = Bwd(n,j+1) + scoring.gapB(n,j+1);
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
	    if (deletion_allowed(i+1,j)) {
		Bwd(i,j) = std::max(Bwd(i,j), Bwd(i+1,j)+scoring.gapA(i+1,j));
	    }
	    if (insertion_allowed(i,j+1)) {
		Bwd(i,j) = std::max(Bwd(i,j), Bwd(i,j+1)+scoring.gapB(i,j+1));
	    }
	}
    }

    std::cout << "Bwd"<<std::endl<<Bwd <<std::endl;
    
    // ----------------------------------------
    // Combine matrices and prune
    //
        
    for(size_type i=1; i<=n; i++) {
	for(size_type j=1; j<=m; j++) {      
	    infty_score_t score_min = (infty_score_t) Score.min();
	    
	    // prune M variable for match i~j
	    if (match_allowed(i,j)) {
		infty_score_t ubm = Fwd(i-1,j-1)+ub_match(i,j)+Bwd(i,j);
		
		if ( (!ubm.is_finite()) || (ubm < score_min)) {
		    ret |= M[i].nq(home,(int)j);
		}
	    }
		
	    // prune G variable for gap i~- (deletion of i between j and j+1)
	    if (deletion_allowed(i,j)) {
		infty_score_t ubd = Fwd(i-1,j) + scoring.gapA(i+1,j) + Bwd(i,j);
		if ( (!ubd.is_finite()) || (ubd < score_min)) {
		    ret |= G[i].nq(home,(int)j);
		}
	    }
	    
	    // prune H variable for gap -~j (insertion of j between i and i+1)
	    if (insertion_allowed(i,j)) {
		infty_score_t ubi = Fwd(i,j-1) + scoring.gapB(i,j+1) + Bwd(i,j);
		if ( (!ubi.is_finite()) || (ubi < score_min)) {
		    ret |= H[i].exclude(home,(int)j);
		}
	    }
	}
    }

    GECODE_ME_CHECK(ret);
    
    // test whether all vars fixed, then subsume (can we subsume earlier?)
    if ( all_vars_fixed() ) {
	Score.eq(home,Score.max());
	return ES_SUBSUMED(*this,dispose(home)); 
    }
    
    return Gecode::me_modified(ret)?Gecode::ES_NOFIX:Gecode::ES_FIX;
}

