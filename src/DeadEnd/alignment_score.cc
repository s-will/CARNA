
#include "alignment_score.hh"
#include "LocARNA/arc_matches.hh"

#include "WinDisplay.cpp"
#include "RNAalignment.hh"

#include <limits>

#include <assert.h>

// const bool debug_out=false;
const bool debug_out=true;

//#include <iostream>

// DETAIL: for the upper bound in ub_match, arc match scores to the
// left and to the right are added. In total, thus each arcmatch score is added twice.
// Dividing the score by two may cause rounding problems. Thus, we multiply all other
// scores by two.
//
// NOTE: it does not work to add only scores either to the left or 
// to the right (depending on whether doing forward or backward computation!)


AlignmentScore::AlignmentScore(Gecode::Space& home,
			       const Sequence &seqA_,
			       const Sequence &seqB_, 
			       const ArcMatches &arc_matches_,
			       const AlignerParams &params_,
			       const Scoring &scoring_,
			       IntViewArray &M_,
			       IntViewArray &T_,
			       Gecode::Int::IntView &Score_
			       )
    : Propagator(home),
      seqA(seqA_),
      seqB(seqB_),
      arc_matches(arc_matches_),
      params(params_),
      scoring(scoring_),
      M(M_),
      T(T_),
      Score(Score_),
      undef(seqB.length()+1)
{
    // subscribe the propagator to any variable change 
    M.subscribe(home,*this,Gecode::Int::PC_INT_DOM);
    T.subscribe(home,*this,Gecode::Int::PC_INT_DOM);
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
    T.update(home,share,p.T);
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
     Gecode::IntVarArray &T,
     Gecode::IntVar &Score
     )
{

    // the post method converts Vars to Views and constructs a new propagator
    // in the home-space

    // translate vars/var arrays to views/view arrays
	
    Gecode::Int::IntView ScoreView = Score;
      
    Gecode::VarArgArray<Gecode::IntVar> MArg = M;
    IntViewArray MView = IntViewArray(home,MArg);
    
    Gecode::VarArgArray<Gecode::IntVar> TArg = T;
    IntViewArray TView = IntViewArray(home, TArg);
    
    new (home) AlignmentScore(home,
			      seqA,
			      seqB,
			      arc_matches,
			      params,
			      scoring,
			      MView,TView,
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


/*
score_t
AlignmentScore::ub_match_simple(size_type i, size_type j) const {
    //std::cout << "AlignmentScore::ub_match"<<std::endl;
  
    // compute upper bound for the contribution of matching positions i
    // and j of respective sequences A and B
    // by summing over all possible base pair matchs
    //
    
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
	    score_t amscore = scoring.arcmatch(am);
	    if (amscore>=0 || match_forced(am.arcA().left(),am.arcB().left()) ) {
		bound+=amscore;
	    }
	}
    }
    
    // same for common left ends
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(i,j).begin();
	arc_matches.common_left_end_list(i,j).end() != it; ++it ) {	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	if ( match_allowed(am.arcA().right(),am.arcB().right()) ) { 
	    // if the right ends of arcs arcA and arcB can match 
	    score_t amscore = scoring.arcmatch(am);
	    if (amscore>=0 || match_forced(am.arcA().right(),am.arcB().right()) ) {
		bound+=amscore;
	    }
	}
    }
    
    return (bound);
}

// simple ub_match: 46 16 68598 68613 
*/


template<class AdjList>
score_t
AlignmentScore::bound_arcmatches(size_t i, size_t j, AdjList adjlA, AdjList adjlB, bool right) const {
    // solve the problem of finding a maximal clique in the compatibility graph of arcs to the left(right)
    // by dynamic programming
    //
    
    // note that the arcs in the right(left) adjacency lists are sorted

    
    // for DP it suffices to store a row of the matrix
    std::vector<infty_score_t> bmat;
    infty_score_t bmat_match; // this is an extra cell for entry ii-1,jj-1 in the matrix b, which is necessary, since we want to overwrite the stored row of bmat
    
    size_type ii=0;
    size_type jj=0;
    
    // for all arcsA and arcsB determine whether its ok to delete them
    // or they are part in a forced arc match
    std::vector<bool> forcedA;
    std::vector<bool> forcedB;
    forcedA.resize(adjlA.size()+1);
    forcedB.resize(adjlB.size()+1);
    fill(forcedA.begin(),forcedA.end(),false);
    fill(forcedB.begin(),forcedB.end(),false);
    
    // for all pairs of arcs in A and B that have right ends i and j, respectively
    //
    ii=0;
    for (typename AdjList::const_iterator arcA=adjlA.begin();
	 arcA!=adjlA.end() ; ++arcA) {
	++ii;
	jj=0; // start with arc number start_jj+1 in sequence B
	for (typename AdjList::const_iterator arcB=adjlB.begin();
	     arcB!=adjlB.end() ; ++arcB) {
	    ++jj;
	    
	    size_t k=(right)?arcA->left():arcA->right();
	    size_t l=(right)?arcB->left():arcB->right();

	    bool forced = match_forced(k,l);
	    
	    forcedA[ii] = forcedA[ii] | forced;
	    forcedB[jj] = forcedB[jj] | forced;
	}
    }
    
    // init bmat; first row
    bmat.resize(adjlB.size()+1);
    bmat[0]=(infty_score_t)0;
    jj=0;
    for (typename AdjList::const_iterator arcB=adjlB.begin();
	 arcB!=adjlB.end() ; ++arcB) {
	++jj;
	//init with score for deleting all arcsB <=jj
	bmat[jj]=bmat[jj-1];
	if (forcedB[jj]) bmat[jj]=infty_score_t::neg_infty;
    }
    
    ii=0;
    // for all pairs of arcs in A and B that have right ends i and j, respectively
    //
    for (typename AdjList::const_iterator arcA=adjlA.begin();
	 arcA!=adjlA.end() ; ++arcA) {
	++ii;
	
	bmat_match=bmat[0];
	if (forcedA[ii]) {
	    bmat[0]=infty_score_t::neg_infty;
	}
	
	jj=0; // start with arc number start_jj+1 in sequence B
	for (typename AdjList::const_iterator arcB=adjlB.begin();
	     arcB!=adjlB.end() ; ++arcB) {
	    ++jj;
	    
	    assert(jj<bmat.size());
	    
	    infty_score_t entry=infty_score_t::neg_infty;

	    size_t k=(right)?arcA->left():arcA->right();
	    size_t l=(right)?arcB->left():arcB->right();

	    if ( match_allowed(k,l) ) {
		score_t amscore = scoring.arcmatch(*arcA,*arcB);
		entry = bmat_match + amscore;
	    }
	    if (!forcedA[ii]) {
		entry = std::max(entry,bmat[jj]); // do not match arcA
	    }
	    if (!forcedB[jj]) {
		entry = std::max(entry,bmat[jj-1]); // do not match arcB
	    }
	    
	    bmat_match=bmat[jj]; // remember entry bmat[jj], which is still for the last row before overwriting it
	    bmat[jj]=entry;
	}
    }
    
    // maximal bound for arcs matches to the left is in bmat[jj] 
    //std::cout <<adjlA.size()<<" "<<adjlB.size()<<" "<<ii<<" "<<jj<<" "<< bmat[jj]<<std::endl;
    return bmat[jj].finite_value();
}


score_t
AlignmentScore::ub_match(size_type i, size_type j) const {
    //std::cout << "AlignmentScore::ub_match "<<i<< " "<<j<<std::endl;
  
    // compute upper bound for the contribution of matching positions i
    // and j of respective sequences A and B
    // by summing over all possible base pair matchs
    //
    // when selecting a certain structure: maximize over possible base pair matchs
    
    score_t bound=2*scoring.basematch(i,j);    
    
    //std::cout << bound<<std::endl;

    bound += bound_arcmatches
	(i,j,
	 arc_matches.get_base_pairsA().right_adjlist(i),
	 arc_matches.get_base_pairsB().right_adjlist(j),
	 true);
    
    bound += bound_arcmatches
	(i,j,
	 arc_matches.get_base_pairsA().left_adjlist(i),
	 arc_matches.get_base_pairsB().left_adjlist(j),
	 false);
    
    //std::cout << "end AlignmentScore::ub_match "<<i<< " "<<j<<" "<<bound<<std::endl;
    
    return bound;
}




bool
AlignmentScore::all_vars_fixed() const {
    // test all variable views and return false as soon as one is not fixed/assigned
    // otherwise return true
    
    //M
    for (size_type i=1; i<(size_t)M.size(); ++i) {
	if (!M[i].assigned()) return false;
    }

    if (!Score.assigned()) return false;
    
    return true;
}


score_t 
AlignmentScore::evaluate_tracematch(const std::vector<size_type> &traceA,
				    const std::vector<size_type> &traceB,
				    size_type i,
				    size_type j) const{
    
    score_t matchscore=2*scoring.basematch(i,j);
    
    if (debug_out) std::cout <<"[ "<<matchscore<<" ";
    // for all pairs of arcs in A and B that have right ends i and j, respectively
    //
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
	arc_matches.common_right_end_list(i,j).end() != it; ++it ) {
	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	if ( traceA[am.arcA().left()] == am.arcB().left() ) { 
	    // if the left ends of arcs arcA and arcB match 
	    if (debug_out) std::cout<<am.arcA().left()<<","<<am.arcB().left()<<":"<<scoring.arcmatch(am)<<" ";
	    matchscore += scoring.arcmatch(am);
	}
    }

    if (debug_out) std::cout <<" ; ";
    
    // same for common left ends
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(i,j).begin();
	arc_matches.common_left_end_list(i,j).end() != it; ++it ) {	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	if ( traceA[am.arcA().right()] == am.arcB().right() ) { 
	    // if the right ends of arcs arcA and arcB  match 
	    if (debug_out) std::cout<<am.arcA().right()<<","<<am.arcB().right()<<":"<<scoring.arcmatch(am)<<" ";
	    matchscore += scoring.arcmatch(am);
	}
    }

    if (debug_out) std::cout <<"]";

    
    return matchscore;    
}

score_t
AlignmentScore::evaluate_trace(const std::vector<size_type> &traceA,
			       const std::vector<size_type> &traceB) const {
    
    const size_t n=seqA.length();
    const size_t m=seqB.length();

    score_t score=0;

    size_type i=1;
    size_type j=1;
    
    if (debug_out) std::cout << "TRACE SCORE "<<n<<" "<<m<<std::endl;
    while (i<=n || j<=m) {
	if (i<=n && traceA[i]==0) {
	    score += 2*scoring.gapA(i,j-1);
	    if (debug_out) std::cout << "del "<<i<<std::endl;
	    ++i;
	} else if (j<=m && traceB[j]==0) {
	    score += 2*scoring.gapB(i-1,j);
	    if (debug_out) std::cout << "ins "<<i<<std::endl;
	    ++j;
	} else {
	    // match between traceB[j]==i and traceA[i]==j

	    score_t matchscore = evaluate_tracematch(traceA,traceB,i,j);
	    
	    score += matchscore;

	    if (debug_out) std::cout << "match "<<i<<" "<<j<<" "<<matchscore<<std::endl;
	    
	    ++i;
	    ++j;
	}
    }
    
    return score;
}

Gecode::ModEvent
AlignmentScore::simple_consistency(Gecode::Space& home) {
    const int n=seqA.length();
    const int m=seqB.length();

    if (debug_out) {
	std::cout << "Before consistency"<<std::endl;
	std::cout << "Matches:    " << M << std::endl;
	std::cout << "Trace:      " << T << std::endl;
	std::cout << "Score:      " << Score << std::endl;
    }

    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
    
    if (!T[0].in(0) || !T[n].in(m)) ret=Gecode::ME_GEN_FAILED;

    if (debug_out) {    
	std::cout << "After consistency"<<std::endl;
	std::cout << "Matches:    " << M << std::endl;
	std::cout << "Trace:      " << T << std::endl;
	std::cout << "Score:      " << Score << std::endl;
    }

    return ret;
}


void
AlignmentScore::forward_algorithm(Gecode::Space& home, Matrix<infty_score_t> &Fwd) {
    const size_t n=seqA.length();
    const size_t m=seqB.length();

     // ----------------------------------------
    // Forward algorithm
    //
    
    Fwd.fill(infty_score_t::neg_infty);
  
    // Fwd(i,j) := max score of a alignment R[1..i], S[1..j]
  
    // initialize
    Fwd(0,0)=(infty_score_t)0;
    for(size_type i=1; i<=n; i++) {
	if (deletion_allowed(i,0)) {
	    Fwd(i,0) = Fwd(i-1,0) + 2*scoring.gapA(i,0);
	}
    }
    
    for(size_type j=1; j<=m; j++) {
	if (insertion_allowed(0,j)) {
	    Fwd(0,j) = Fwd(0,j-1)  + 2*scoring.gapB(0,j);
	}
    }
  
    // recurse
    for(size_type i=1; i<=n; i++) {
	for(size_type j=1; j<=m; j++) {
	    if (match_allowed(i,j)) {
		Fwd(i,j) = Fwd(i-1,j-1) + ub_match(i,j);
	    }
	    if (deletion_allowed(i,j)) {
		Fwd(i,j) = std::max(Fwd(i,j),Fwd(i-1,j)+2*scoring.gapA(i,j));
	    }
	    if (insertion_allowed(i,j)) {
		Fwd(i,j) = std::max(Fwd(i,j),Fwd(i,j-1)+2*scoring.gapB(i,j));
	    }
	}
    }
}

void
AlignmentScore::
backtrace_forward(Gecode::Space &home, const Matrix<infty_score_t> &Fwd,
		  vector<size_type> &traceA,
		  vector<size_type> &traceB
		  ) {
    const int n=seqA.length();
    const int m=seqB.length();

    // ----------------------------------------
    // do backtracking of the Fwd matrix in order to obtain a lower
    // bound of alignment score and guide the heuristic
    
    traceA.resize(n+1);
    traceB.resize(m+1);
    
    {
	int i=n;
	int j=m;
	
	while (i>0 && j>0) {
	    if ( match_allowed(i,j) && Fwd(i,j) == Fwd(i-1,j-1) + ub_match(i,j) ) {
		traceA[i]=j;
		traceB[j]=i;
		--i; 
		--j;
	    } else if ( deletion_allowed(i,j) && Fwd(i,j) == Fwd(i-1,j) + 2*scoring.gapA(i,j) ) {
		traceA[i]=0;
		--i;
	    } else if ( insertion_allowed(i,j) && Fwd(i,j) == Fwd(i,j-1) + 2*scoring.gapB(i,j) ) {
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
}

void
AlignmentScore::backward_algorithm(Gecode::Space& home, Matrix<infty_score_t> &Bwd) {
    const int n=seqA.length();
    const int m=seqB.length();

    // ----------------------------------------
    // Backward algorithm
    //
    
    // Bwd(i,j) := max score of a alignment R[i+1..n], S[j+1..m]
    
    // init with "invalid" value negative infinity
    Bwd.fill(infty_score_t::neg_infty);
    
    // initialize
    Bwd(n,m)=(infty_score_t)0;
    for(size_type i=n; i>0;) { // for i=n-1 downto 0
	--i;
	if (deletion_allowed(i+1,m)) {
	    Bwd(i,m) = Bwd(i+1,m) + 2*scoring.gapA(i+1,m);
	}
    }
    for(size_type j=m; j>0;) { // for j=m-1 downto 0
	--j;
	if (insertion_allowed(n,j+1)) {
	    Bwd(n,j) = Bwd(n,j+1) + 2*scoring.gapB(n,j+1);
	}
    }
    
    // recurse
    for(size_type i=n; i>0;) { // for i=n-1 downto 0
	--i;
	for(size_type j=m; j>0;) { // for j=m-1 downto 0
	    --j;
	    if (match_allowed(i+1,j+1)) {   // match i+1 ~ j+1 (!)
		Bwd(i,j) = Bwd(i+1,j+1) + ub_match(i+1,j+1);
	    }
	    if (deletion_allowed(i+1,j)) { // delete i+1 after j
		Bwd(i,j) = std::max(Bwd(i,j), Bwd(i+1,j)+2*scoring.gapA(i+1,j));
	    }
	    if (insertion_allowed(i,j+1)) { // insert j+1 after i
		Bwd(i,j) = std::max(Bwd(i,j), Bwd(i,j+1)+2*scoring.gapB(i,j+1));
	    }
	}
    }
}

Gecode::ModEvent
AlignmentScore::prune(Gecode::Space& home, 
		      const Matrix<infty_score_t> &Fwd,
		      const Matrix<infty_score_t> &Bwd) {
    // ----------------------------------------
    // Combine matrices and prune
    //

    const size_t n=seqA.length();
    const size_t m=seqB.length();
    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
    
    for(size_type i=1; i<=n; i++) {
	for(size_type j=1; j<=m; j++) {
	    // prune M variable for match i~j
	    if (match_allowed(i,j)) {
		infty_score_t ub = Fwd(i-1,j-1)+ub_match(i,j)+Bwd(i,j);
		
		if ( (!ub.is_finite()) || (ub < (infty_score_t) Score.min())) {
		    ret |= M[i].nq(home,(int)j);
		}
	    }
	}
    }

    for(size_type i=1; i<=n; i++) {
	for(size_type j=0; j<=m; j++) {      
	    
	    // prune T variable
	    infty_score_t ub = Fwd(i,j) + Bwd(i,j);
	    if ( (!ub.is_finite()) || (ub < (infty_score_t) Score.min())) {
		ret |= T[i].nq(home,(int)j);
	    }
	}
    }

    return ret;
}

void
AlignmentScore::choice(RNAalignment &s,
		       const Matrix<infty_score_t> &Fwd,
		       const Matrix<infty_score_t> &Bwd,
		       const std::vector<size_type> &traceA,
		       const std::vector<size_type> &traceB) const {
    
    const size_t n=seqA.length();
    const size_t m=seqB.length();

    // select T[pos] with largest domain
    size_t maxs=0;
    size_t pos=0;
    for (size_t i=0; i<=n; i++) {
	if (T[i].size()>maxs) {
	    pos=i;
	    maxs=T[i].size();
	}
    }
    assert(maxs!=0);
    
    size_t val=traceA[pos];
    if (pos==0 || val==0) val=undef;
    
    // copy choice to the space
    s.pos=pos;
    s.val=val;
    
    // determine values for splitting choice
    
    score_t minb=numeric_limits<score_t>::max();
    score_t maxb=numeric_limits<score_t>::min();
    
    std::cout << "CHOICE AT "<<pos <<std::endl;
    for (size_t j=0; j<=m; j++) {
	infty_score_t b;
	if (s.T[pos].in(j) && (b=Fwd(pos,j)+Bwd(pos,j)).is_finite()) {
	    minb=min(minb,b.finite_value());
	    maxb=max(maxb,b.finite_value());
	}
	//std::cout <<j<<": "<< b << std::endl;
    }
    
    //std::cout << "minb: "<<minb<<" maxb: "<<maxb<<" T[pos]: "<<s.T[pos]<<std::endl;

    size_t minval=m+1;
    size_t maxval=0;
    
    for (size_t j=0; j<=m; j++) {
	infty_score_t b;
	if (s.T[pos].in(j) && (b=Fwd(pos,j)+Bwd(pos,j)).is_finite() ) {
	    //std::cout << j << ":" << b<<" ";
	    if (b.finite_value()>=minb+0.5*(maxb-minb)) { 
		minval=min(minval,j);
		maxval=max(maxval,j);
	    }
	}
    }
    
    if (minval==maxval) {
	std::cout << " IN: " << minval;
    } else {
	std::cout << " IN: " << minval << "-" << maxval;
    }
    std::cout << std::endl;

    if (minval==(size_t)s.T[pos].min() && maxval==(size_t)s.T[pos].max()) {
	maxval=(maxval-minval)/2+minval;
    }
    
    
    s.minval=minval;
    s.maxval=maxval;
}


Gecode::ExecStatus
AlignmentScore::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
    
    if (debug_out) {
	std::cout << "AlignmentScore::propagate " <<std::endl;
    }

    // ----------------------------------------
    // Propagation strategy:
    //
    // This is the naive implementation of the propagator!
    // The idea is strictly following "Keep It Simple and Stupid" (KISS)
    //
    // propagation is performed in three steps:
    // 1) forward algo 2) backward algo 3) compute edge upper bounds
    // and prune domains in 1) and 2) we fill the respective matrices
    // Fwd and Bwd. The combination of matrix entries
    // Fwd(i-1,j-1)+ub_match(i,j)+Bwd(i,j) yields the upper bound for
    // match i~j.
    //
    // ----------------------------------------

    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;

    // -------------------- RUN SIMPLE CONSISTENCY
    ret |= simple_consistency(home);
 
    if (Gecode::me_failed(ret)) {
	if (debug_out) std::cout << "simple consistency failed"<<std::endl;
	return Gecode::ES_FAILED;
    }
    
    const size_t n=seqA.length();
    const size_t m=seqB.length();

    // -------------------- RUN FORWARD ALGORITHM
    Matrix<infty_score_t> Fwd(n+1,m+1); // ForWarD matrix
    
    forward_algorithm(home,Fwd);
    
    if (debug_out) {std::cout << "Fwd"<<std::endl<<Fwd <<std::endl;}
    
    // the forward algorithm yields an upper bound
    // ==> constrain the score variable
    if (Fwd(n,m).is_finite()) {
	if (Gecode::me_failed(Score.lq(home,(int)Fwd(n,m).finite_value()))) {
	    if (debug_out) std::cout << "Bounding after Fwd failed."<<std::endl;
	    return Gecode::ES_FAILED;
	}
    } else {
	if (debug_out) std::cout << "Fwd has infinite bound: fail."<<std::endl;
	return Gecode::ES_FAILED;
    }

    // trace vectors. idea one entry per position gives position in
    // the other sequence to that the pos is matched or 0 (in case of
    // gap).
    std::vector< size_type > traceA;
    std::vector< size_type > traceB;

    // -------------------- BACKTRACE FWD
    backtrace_forward(home,Fwd,traceA,traceB);
    
    
    if (debug_out) {
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
    }

    // -------------------- Lower bounding

    // evaluate alignment for obtaining lower bound
    score_t trace_score = evaluate_trace(traceA,traceB);
    
    if (debug_out) {
	std::cout << "TRACE Score: " << trace_score << std::endl;
    }
    
    // constrain Score by trace_score, which is a lower bound
    
    if (Gecode::me_failed(Score.gq(home,(int)trace_score))) {
     	if (debug_out) std::cout << "Setting of new lower bound failed."<<std::endl;
      	return Gecode::ES_FAILED;
    }
    
    // -------------------- RUN BACKWARD ALGORITHM
    Matrix<infty_score_t> Bwd(n+1,m+1); 

    backward_algorithm(home,Bwd);
    //if (debug_out) {std::cout << "Bwd"<<std::endl<<Bwd <<std::endl;}


    if (debug_out) {
	// ------------------------------------------------------------
	// DEBUGGING: check trace
	for (size_t i=1; i<=n; i++) {
	    std::cout << "TRACE " << i <<" " << traceA[i] <<" ";
	    if (traceA[i]!=0) {
		size_t j = traceA[i];
		score_t matchscore=evaluate_tracematch(traceA,traceB,i,j);
		
		cout
		    << match_allowed(i,j)<<" "
		    <<Fwd(i-1,j-1)<<"+"<< matchscore<<"["<<ub_match(i,j)<<"]"<<"+"<<Bwd(i,j)<<" = "
		    <<Fwd(i-1,j-1)+matchscore+Bwd(i,j);
	    }
	    std::cout << std::endl;
	}
    }
    
 
    // -------------------- PRUNE VARIABLES
    ret |= prune(home,Fwd,Bwd);


    
    // -------------------- select CHOICE for the space
    choice(static_cast<RNAalignment&>(home),Fwd,Bwd,traceA,traceB);

    if (debug_out) {
	std::cout << "After Pruning:"<<std::endl;
	std::cout << "Matches:    " << M << std::endl;
	std::cout << "Trace:      " << T << std::endl;
	std::cout << "Score:      " << Score << std::endl;
    }
    
    if (Gecode::me_failed(ret)) {
	if (debug_out) std::cout << "Fail after pruning."<<std::endl;
	return Gecode::ES_FAILED;
    }
    
    // test whether all vars fixed, then subsume (can we subsume earlier?)
    if ( all_vars_fixed() ) {
	Score.eq(home,Score.max());
	return ES_SUBSUMED(*this,dispose(home)); 
    }
    
    return Gecode::me_modified(ret)?Gecode::ES_NOFIX:Gecode::ES_FIX;
}
