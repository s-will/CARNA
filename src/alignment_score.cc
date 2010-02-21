
#include "alignment_score.hh"
#include "LocARNA/arc_matches.hh"

#include "WinDisplay.hh"
#include "RNAalignment.hh"

#include <limits>

#include <assert.h>


/*
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IDEA OF TIGHTNESS IS NOT WORKING
  
  TO MAKE IT WORK:
  Score arc matchs only for non-guaranteed matches
  For all guaranteed matches transfer their arc match score to the bsae match score of
  the other arc match end
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */



const bool debug_out=false;
//const bool debug_out=true;

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
			       IntViewArray &MD_,
			       BoolViewArray &M_,
			       Gecode::Int::IntView &Score_
			       )
    : Propagator(home),
      seqA(seqA_),
      seqB(seqB_),
      arc_matches(arc_matches_),
      params(params_),
      scoring(scoring_),
      MD(MD_),
      M(M_),
      Score(Score_)
{
    // subscribe the propagator to any variable change 
    MD.subscribe(home,*this,Gecode::Int::PC_INT_DOM);
    M.subscribe(home,*this,Gecode::Int::PC_INT_DOM);
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
    scoring(p.scoring)
{
    MD.update(home,share,p.MD);
    M.update(home,share,p.M);
    Score.update(home,share,p.Score);
}

AlignmentScore::~AlignmentScore() {}

Gecode::ExecStatus
AlignmentScore::
post(Gecode::Space& home,
     const Sequence &seqA,
     const Sequence &seqB,
     const ArcMatches &arc_matches,
     const AlignerParams &params,
     const Scoring &scoring,
     Gecode::IntVarArray &MD,
     Gecode::BoolVarArray &M,
     Gecode::IntVar &Score
     )
{

    // the post method converts Vars to Views and constructs a new
    // propagator in the home-space

    // translate vars/var arrays to views/view arrays
	
    Gecode::Int::IntView ScoreView = Score;
      
    Gecode::VarArgArray<Gecode::IntVar> MDArg = MD;
    IntViewArray MDView = IntViewArray(home,MDArg);
    
    Gecode::VarArgArray<Gecode::BoolVar> MArg = M;
    BoolViewArray MView = BoolViewArray(home, MArg);
        
    new (home) AlignmentScore(home,
			      seqA,
			      seqB,
			      arc_matches,
			      params,
			      scoring,
			      MDView,MView,
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

template<class AdjList>
void
AlignmentScore::determine_forced_arcs(const AdjList &adjlA,
				      const AdjList &adjlB,
				      std::vector<bool> &forcedA, 
				      std::vector<bool> &forcedB,
				      bool right,
				      bool *tight) const {
    // for all pairs of arcs in A and B that have right ends i and j,
    // respectively
    //
    size_t i=1;
    for (typename AdjList::const_iterator arcA=adjlA.begin();
	 arcA!=adjlA.end() ; ++arcA, ++i) {
	
	size_t j=1;
	for (typename AdjList::const_iterator arcB=adjlB.begin();
	     arcB!=adjlB.end() ; ++arcB, ++j) {
	    
	    size_t k=(right)?arcA->left():arcA->right();
	    size_t l=(right)?arcB->left():arcB->right();
	    
	    bool forced = match_forced(k,l);

	    if ( (!forced) && match_allowed(k,l) ) { // if the match
						     // is allowed,
						     // but not forced
		*tight = false; // the bound will not be tight
	    }
	    
	    forcedA[i] = forcedA[i] | forced;
	    forcedB[j] = forcedB[j] | forced;
	}
    }
}

template<class AdjList>
score_t
AlignmentScore::bound_arcmatches(const AdjList &adjlA,
				 const AdjList &adjlB,
				 const Matrix<bool> &considered_ams,
				 bool right,
				 bool *tight_ptr) const {
    
    // short cut for empty adjacency lists
    if (adjlA.size()==0 || adjlB.size()==0) {
	return 0;
    }

    // solve the problem of finding a maximal clique in the
    // compatibility graph of arcs to the left(right) by dynamic
    // programming
    //
    
    // note that the arcs in the right(left) adjacency lists are sorted
    
    
    // for DP it suffices to store a row of the matrix
    std::vector<infty_score_t> MAT(adjlB.size()+1);
    infty_score_t MAT_match; // this is an extra cell for entry
			      // i-1,j-1 in the matrix b, which is
			      // necessary, since we want to overwrite
			      // the stored row of MAT
    
    size_type i;
    size_type j;
    
    // init MAT; first row
    MAT[0]=(infty_score_t)0;
    for (size_t j=1; j<=adjlB.size(); ++j) {
	//init with score for deleting all arcsB <=j
	MAT[j]=MAT[j-1];
    }
    
    // for all pairs of arcs in A and B that have right ends i and j, respectively
    //
    i=1;
    for (typename AdjList::const_iterator arcA=adjlA.begin();
	 arcA!=adjlA.end() ; ++arcA, ++i) {
	
	MAT_match=MAT[0];
	
	j=1;
	for (typename AdjList::const_iterator arcB=adjlB.begin();
	     arcB!=adjlB.end() ; ++arcB, ++j) {

	    //size_t k=(right)?arcA->left():arcA->right();
	    //size_t l=(right)?arcB->left():arcB->right();
	    
	    
	    infty_score_t entry=infty_score_t::neg_infty;
	    if ( considered_ams.get(arcA->idx(), arcB->idx()) ) { // note: for considered am the match k,l is allowed
		entry = MAT_match + scoring.arcmatch(*arcA,*arcB);
	    }
	   
	    entry = std::max(entry,MAT[j]); // do not match arcA
	    entry = std::max(entry,MAT[j-1]); // do not match arcB
	    	    
	    MAT_match=MAT[j]; // remember entry MAT[j], which is still
			      // a matrix entry of the last row,
			      // before overwriting it
	    MAT[j]=entry;
	}
    }

    score_t bound = MAT[adjlB.size()].finite_value();

    // when is the bound TIGHT? A: if there are no considered arc matches
    // contributing to the bound.
    // Since arc matches with negative score do never contribute,
    // the bound is tight iff the bound is 0!!!
    bool tight=(bound==0);

    // RETURN tight
    if (tight_ptr!=NULL && *tight_ptr==true) {
	*tight_ptr = tight;
    }
    
    // RETURN bound
    // maximal bound for arcs matches to the left is in MAT[adjlB.size()] 
    //std::cout <<adjlA.size()<<" "<<adjlB.size()<<" "<< MAT[adjlB.size()]<<std::endl;
    return bound;
}


score_t
AlignmentScore::ub_match(size_type i, size_type j,
			 const Matrix<bool> &considered_ams,
			 const Matrix<score_t> &match_scores,
			 bool *tight) const {
    //std::cout << "AlignmentScore::ub_match "<<i<< " "<<j<<std::endl;
  
    // compute upper bound for the contribution of matching positions i
    // and j of respective sequences A and B
    // by summing over all possible base pair matchs
    //
    // when selecting a certain structure: maximize over possible base pair matchs
    
    score_t bound=match_scores(i,j);
    
    //std::cout << bound<<std::endl;

    bound += 
	bound_arcmatches(arc_matches.get_base_pairsA().right_adjlist(i),
			 arc_matches.get_base_pairsB().right_adjlist(j),
			 considered_ams,
			 true,tight);
    
    bound += 
	bound_arcmatches(arc_matches.get_base_pairsA().left_adjlist(i),
			 arc_matches.get_base_pairsB().left_adjlist(j),
			 considered_ams,
			 false,tight);
    
    //std::cout << "end AlignmentScore::ub_match "<<i<< " "<<j<<" "<<bound<<std::endl;
    
    return bound;
}


bool
AlignmentScore::all_vars_fixed() const {
    // test all variable views and return false as soon as one is not fixed/assigned
    // otherwise return true
    
    //MD
    for (size_type i=0; i<(size_t)MD.size(); ++i) {
	if (!MD[i].assigned()) return false;
    }

    //M
    for (size_type i=0; i<(size_t)M.size(); ++i) {
	if (!M[i].assigned()) return false;
    }

    return Score.assigned();
}


score_t 
AlignmentScore::evaluate_tracematch(const std::vector<size_type> &traceA,
				    const std::vector<size_type> &traceB,
				    const Matrix<bool> &considered_ams,
				    const Matrix<score_t> &match_scores,
				    size_type i,
				    size_type j) const{
    
    score_t matchscore=match_scores(i,j);

    // for all pairs of arcs in A and B that have right ends i and j, respectively
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
	arc_matches.common_right_end_list(i,j).end() != it; ++it ) {
	
	const ArcMatch &am = arc_matches.arcmatch(*it);

	if (!considered_ams.get(am.arcA().idx(),am.arcB().idx())) continue;
	
	if ( traceA[am.arcA().left()] == am.arcB().left() ) { 
	    // if the left ends of arcs arcA and arcB match 
	    matchscore += scoring.arcmatch(am);
	}
    }

    // same for common left ends
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(i,j).begin();
	arc_matches.common_left_end_list(i,j).end() != it; ++it ) {	
	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	if (!considered_ams.get(am.arcA().idx(),am.arcB().idx())) continue;
	
	if ( traceA[am.arcA().right()] == am.arcB().right() ) { 
	    // if the right ends of arcs arcA and arcB  match 
	    matchscore += scoring.arcmatch(am);
	}
    }

    return matchscore;
}

score_t
AlignmentScore::evaluate_trace(const std::vector<size_type> &traceA,
			       const std::vector<size_type> &traceB,
			       const Matrix<bool> &considered_ams,
			       const Matrix<score_t> &match_scores
			       ) const {
    
    const size_t n=seqA.length();
    const size_t m=seqB.length();

    score_t score=0;

    size_type i=1;
    size_type j=1;
    
    // while not all positions i and j consumed
    while (i<=n || j<=m) {
	// invariant: score is the score for aligning A_1..i-1 to B_1..j-1 according to the trace 
	if (i<=n && traceA[i]==0) { // deletion of i after j-1
	    score += 2*scoring.gapA(i,j-1);
	    ++i;
	} else if (j<=m && traceB[j]==0) { // insertion of j after i-1
	    score += 2*scoring.gapB(i-1,j);
	    ++j;
	} else { // match between i and j
	    assert(j==traceA[i]);
	    assert(i==traceB[j]);
	    
	    score += evaluate_tracematch(traceA,traceB,considered_ams,match_scores,i,j);
	    ++i;
	    ++j;
	}
    }
    
    return score;
}

void
AlignmentScore::print_vars() const {
    std::cout << "Matches/Deletions:    ";
    for (size_t i=0; i<=seqA.length(); i++) {
	std::cout <<i<<(M[i].assigned()?(M[i].val()==0?"g":"~"):"?")<<MD[i]<<", "; 
    }
    
    //std::cout << "Matches/Deletions:    " << MD << std::endl;
    //std::cout << "Match Flags:          " << M << std::endl;
    std::cout << "Score:      " << Score << std::endl;
}



void
AlignmentScore::forward_algorithm(Gecode::Space& home,
				  Matrix<infty_score_t> &Fwd,
				  const Matrix<score_t> &UBM
				  ) {
    const size_t n=seqA.length();
    const size_t m=seqB.length();

     // ----------------------------------------
    // Forward algorithm
    //
    
    //Fwd.fill(infty_score_t::neg_infty);
    
    // Definition( Fwd-matrix )
    // Fwd(i,j) := max score of a alignment R[1..i], S[1..j]
    
    ////////////////////
    // initialize
    Fwd(0,0)=(infty_score_t)0;
    for(size_type i=1; i<=n; i++) {
	if (deletion_allowed(i,0)) {
	    Fwd(i,0) = Fwd(i-1,0) + 2*scoring.gapA(i,0);
	} else {
	    Fwd(i,0) = infty_score_t::neg_infty;
	}
    }
    
    for(size_type j=1; j<=m; j++) {
	if (insertion_allowed(0,j)) {
	    Fwd(0,j) = Fwd(0,j-1)  + 2*scoring.gapB(0,j);
	} else {
	    Fwd(0,j) = infty_score_t::neg_infty;
	}
    }
    
    // special initialization:
    // in each matrix row i the entry Fwd(i,MD[i].min()-1) is invalid,
    // but could be accessed in the recursion
    for(size_type i=1; i<=n; i++) {
	if (MD[i].min()>1) {
	    Fwd(i,MD[i].min()-1) = infty_score_t::neg_infty;
	}
    }
    
    Fwd(n,m) = infty_score_t::neg_infty;
    // end init
    ////////////////////
    // recurse
    for(size_type i=1; i<=n; i++) {
	size_t maxj=m;
	if (i<n) {
	    maxj = MD[i+1].max();
	    if (!M[i+1].in(0)) {
		maxj--;
	    }
	}
	
	for(size_type j=max(1,MD[i].min()); j<=maxj; j++) {
	    Fwd(i,j) = infty_score_t::neg_infty;
	    if (match_allowed(i,j)) {
		Fwd(i,j) = Fwd(i-1,j-1) + UBM(i,j);
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
		  const Matrix<score_t> &UBM,
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
	    if ( match_allowed(i,j) && Fwd(i,j) == Fwd(i-1,j-1) + UBM(i,j) ) {
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
AlignmentScore::backward_algorithm(Gecode::Space& home, Matrix<infty_score_t> &Bwd,
				   const Matrix<score_t> &UBM) {
    const size_t n=seqA.length();
    const size_t m=seqB.length();

    // ----------------------------------------
    // Backward algorithm
    //
    
    // Bwd(i,j) := max score of a alignment R[i+1..n], S[j+1..m]
    
    // init with "invalid" value negative infinity
    // Bwd.fill(infty_score_t::neg_infty);
    
    // initialize
    Bwd(n,m)=(infty_score_t)0;
    for(size_type i=n; i>0;) { // for i=n-1 downto 0
	--i;
	if (deletion_allowed(i+1,m)) {
	    Bwd(i,m) = Bwd(i+1,m) + 2*scoring.gapA(i+1,m);
	} else {
	    Bwd(i,m) = infty_score_t::neg_infty;
	}
    }
    for(size_type j=m; j>0;) { // for j=m-1 downto 0
	--j;
	if (insertion_allowed(n,j+1)) {
	    Bwd(n,j) = Bwd(n,j+1) + 2*scoring.gapB(n,j+1);
	} else {
	    Bwd(n,j) = infty_score_t::neg_infty;
	}
    }
    
    Bwd(0,0) = infty_score_t::neg_infty;
    
    // special initialization
    for(size_type i=n; i>0;) { // for i=n-1 downto 0
	--i;

	size_t maxj = MD[i+1].max();
	if (!M[i+1].in(0)) {
	    maxj--;
	}
	
	if (maxj+1<m) {
	    Bwd(i,maxj+1) = infty_score_t::neg_infty;
	}
    }
    
    // recurse
    for(size_type i=n; i>0;) { // for i=n-1 downto 0
	--i;
	
	// we need to cover allowed matches/deletions/insertions in the row i+1!
	size_t minj = MD[i].min();
	
	size_t maxj = MD[i+1].max();
	if (!M[i+1].in(0)) {
	    maxj--;
	}
	maxj = min(maxj,m-1);

	
	for(size_type j=maxj+1; j>minj;) { // for j=m-1 downto 0
	    --j;
	    Bwd(i,j) = infty_score_t::neg_infty;
	    if (match_allowed(i+1,j+1)) {   // match i+1 ~ j+1 (!)
		Bwd(i,j) = Bwd(i+1,j+1) + UBM(i+1,j+1);
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
		      const Matrix<infty_score_t> &Bwd,
		      const Matrix<score_t> &UBM
		      ) {
    // NOTE: using the precomputed upper bounds for matching is weaker than
    // new computation via ub_match (which may give stronger bounds due to prior pruning)
    // Pruning will however anyway cause a re-run of propagate()
    // ==> use precomputed UBM since it's significantly faster!
    
    
    // ----------------------------------------
    // Combine matrices and prune
    //

    const size_t n=seqA.length();
    // const size_t m=seqB.length();
    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;

    ////////////////////////////////////////
    // prune MD and M
    for(size_type i=1; i<=n; i++) {
	for(size_type j=(size_t)MD[i].min(); j<=(size_t)MD[i].max(); j++) {
	    // prune MD variable for match i~j
	    if (match_or_deletion_allowed(i,j)) {
		infty_score_t ub = Fwd(i,j)+Bwd(i,j);
		if ( (!ub.is_finite()) || (ub < (infty_score_t) Score.min())) {
		    ret |= MD[i].nq(home,(int)j);
		}
	    }
	}
    }
    
    ////////////////////////////////////////
    // prune on M

    // prune impossible deletions
    for(size_type i=1; i<=n; i++) {
	bool del=false; // is deletion possible?
	for(size_type j=(size_t)MD[i].min(); !del && j<=(size_t)MD[i].max(); j++) {
	    infty_score_t ubd = Fwd(i-1,j)+2*scoring.gapA(i,j)+Bwd(i,j);
	    if ( ubd.is_finite() && (ubd >= (infty_score_t) Score.min())) {
		del=true;
	    }
	}
	if (!del) {
	    ret |= M[i].nq(home,0); // deletion of i not possible
	}
    }
    
    // prune impossible matches
    for(size_type i=1; i<=n; i++) {
	bool match=false; // is match possible?
	for(size_type j=std::max((size_t)1,(size_t)MD[i].min()); !match && j<=(size_t)MD[i].max(); j++) {
	    infty_score_t ubm = Fwd(i-1,j-1)+UBM(i,j)+Bwd(i,j);
	    if ( ubm.is_finite() && (ubm >= (infty_score_t) Score.min())) {
		match = true; // match is possible
	    }
	}
	if (!match) {
	    M[i].nq(home,1); // match of i not possible
	}
    }
    

    return ret;
}

void
AlignmentScore::choice(RNAalignment &s,
		       const Matrix<infty_score_t> &Fwd,
		       const Matrix<infty_score_t> &Bwd,
		       const std::vector<size_type> &traceA,
		       const std::vector<size_type> &traceB,
		       const Matrix<score_t> &UBM) const {
    
    // NOTE: even more severe than in prune, UBM does not contain the best upper bounds anymore
    // to save time, we use UBM nevertheless instead of re-computing the bounds 
    
    const size_t n=seqA.length();
    const size_t m=seqB.length();
    
    // determine position with largest weight
    // use tie-breaking
    
    if (debug_out) std::cout <<"Determine choice"<<std::endl;
    //if (debug_out) print_vars();
    
    vector<score_t> weights;
    score_t maxweight=numeric_limits<score_t>::min();
    
    weights.resize(n+1);
    for (size_t i=0; i<=n; i++) {
	
	// weight by maximal base pair bound
	weights[i]=numeric_limits<score_t>::min();
	
	if (!s.MD[i].assigned()) {
	    for (size_t j=s.MD[i].min(); j<=(size_t)s.MD[i].max(); j++) {
		if (s.MD[i].in((int)j)) {
		    weights[i]=max(weights[i],
				   UBM(i,j)-2*scoring.basematch(i,j) //ub_match(i,j,false)
				   + 
				   (score_t)s.MD[i].size());
		}
	    }
	}
    }
    
    size_t pos=0;
    
    for (size_t i=0; i<=n; i++) {
	if (weights[i]>maxweight) {
	    maxweight=weights[i];
	    pos=i;
	}
    }
    
    if(maxweight == numeric_limits<score_t>::min()) {	
	return;
    }


    // determine longest run of vars that have weight maxweight
    
    int last_notchosen=-1;
    size_t best_left_end=0;
    size_t best_run_len=0;
    
    for (size_t i=0; i<=n; i++) {
	if (weights[i]<maxweight) {
	    size_t cur_run_len = i-last_notchosen-1;
	    if (cur_run_len>best_run_len) {
		best_run_len=cur_run_len;
		best_left_end=last_notchosen+1;
	    }
	    last_notchosen=i;
	}
    }
    
    size_t cur_run_len = s.MD.size()-last_notchosen-1;
    if (cur_run_len>best_run_len) {
	best_run_len=cur_run_len;
	best_left_end=last_notchosen+1;
    }
    
    assert(best_run_len>0); // otherwise status() is incorrect
    
    pos=best_left_end + (best_run_len/2); // split the longest run
    
    assert(0<=pos && pos<=n);


    // print trace (DEBUGGING)
    //if (debug_out) {std::cout << "TRACE"<<std::endl;
    //	for (size_type i=1; i<=n; ++i ) {
    //	    std::cout << traceA[i] << " ";
    //	}
    //	std::cout<<std::endl;
    //}
    
    size_t val = (pos==0)?0:traceA[pos];
    
    if (val==0) { // trace suggests to delete pos
	// look in traceA, where i is deleted
	for (size_t i=pos; i>0;) {
	    --i;
	    if (traceA[i]!=0) {
		val=traceA[i];
		break;
	    }
	}
	
	if (debug_out) std::cout << "CHOICE AT "<<pos <<" DELETION AFTER "<<val<<" "<<MD[pos] <<std::endl;
	
	if (debug_out) print_vars();
    }
    

    // copy choice to the space
    s.pos=pos;
    s.val=val;
    
    if (debug_out) std::cout << "Choose pos:"<<pos<<" val:"<<val<<std::endl;

    // determine values for splitting choice
    
    score_t minb=numeric_limits<score_t>::max();
    score_t maxb=numeric_limits<score_t>::min();

    if (debug_out) std::cout << "CHOICE AT "<<pos <<std::endl;
    for (size_t j=1; j<=m; j++) {
	infty_score_t b;
	if (s.MD[pos].in(j) && (b=Fwd(pos,j)+Bwd(pos,j)).is_finite()) {
		minb=min(minb,b.finite_value());
		maxb=max(maxb,b.finite_value());
	}
    }

    size_t minval=m+1;
    size_t maxval=0;
    
    for (size_t j=1; j<=m; j++) {
	infty_score_t b;
	if (s.MD[pos].in(j) && (b=Fwd(pos,j)+Bwd(pos,j)).is_finite() ) {
	    if (debug_out) std::cout << j << ":" << b<<" ";
	    if (b.finite_value()>=minb+0.80*(maxb-minb)) { 
		minval=min(minval,j);
		maxval=max(maxval,j);
	    }
	}
    }	    
    
    
    if (minval==(size_t)s.MD[pos].min() && maxval==(size_t)s.MD[pos].max()) {
	
	size_t medval = (maxval-minval)/2+minval;
	
	// choose according to trace
	
	if (val<=medval) { 
	    maxval=medval;
	} else {
	    minval=medval+1;
	}
    }
        
    if (debug_out) { 
	if (minval==maxval) {
	    std::cout << " IN: " << minval;
	} else {
	    std::cout << " IN: " << minval << "-" << maxval;
	}
	std::cout << std::endl;
    }
    
    // copy choice to the space
    s.minval=minval;
    s.maxval=maxval;
    
}


Gecode::ModEvent
AlignmentScore::fix_vars_to_trace(Gecode::Space &home,
				  size_t start,
				  size_t end,
				  const std::vector<size_type> &traceA,
				  const std::vector<size_type> &traceB
				  ) {
    ////////////////////
    // some assertions
    assert(start>=0);
    assert(end<=seqA.length());
    if (start>0) {
	assert(MD[start-1].assigned());
    }
    if (end<seqA.length()) {
	assert(MD[end+1].assigned());
    }
    ////////////////////
    
    
    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
    int last=0;
    if (start>0) {
	last=MD[start-1].val();
    }
    
    for (size_t i=start; i<=end; i++) {
	if (traceA[i]!=0) { // match of i
	    // assert( ub_match(i,traceA[i]) == evaluate_tracematch(traceA,traceB,i,traceA[i]) );
	    ret |= MD[i].eq(home,(int)traceA[i]);
	    ret |= M[i].eq(home,1);
	    last=traceA[i];
	} else { // deletion of i
	    ret |= MD[i].eq(home,last);
	    ret |= M[i].eq(home,0);
	}
    }
    return ret;
}


Gecode::ModEvent
AlignmentScore::fix_tight_runs(Gecode::Space &home,
			       const std::vector<size_type> &traceA,
			       const std::vector<size_type> &traceB,
			       const Matrix<bool> &tight) {
    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
    
    size_t n=seqA.length();
    
    // determine tight runs
    
    int last_assigned=-1;
    bool run_is_tight=true;
    
    for (size_t i=0; i<=n; i++) {
	if (MD[i].assigned() && M[i].assigned()) {// end of a run
	    
	    // if the whole run is tight
	    if (i>(size_t)(last_assigned+1) && run_is_tight) {
		
		// we can fix the whole run
		/*
		print_vars();
		for (size_t k=0; k<=n; k++) {
		    if (traceA[k]!=0) {
			bool t=tight.get(k,traceA[k]);
			std::cout <<k<<"~"<<traceA[k]<<":"<<t<<" ";
		    }
		}
		std::cout << std::endl;
		std::cout << "fix "<<(last_assigned+1)<<"-"<<(i-1)<<std::endl;
		*/
		ret |= fix_vars_to_trace(home,last_assigned+1,i-1,traceA,traceB);
	       
	    }
	    
	    last_assigned=i;
	    run_is_tight=true;
	    
	} else { // within a run
	    if (traceA[i]!=0) {
		run_is_tight &= tight.get(i,traceA[i]);
	    }
	}
    }
    
    
     // if the whole run from last_assigned to the end is tight
    if ((size_t)(last_assigned+1)<n+1 && run_is_tight) {
	// we can fix the whole run
	//std::cout << "fix "<<(last_assigned+1)<<"-"<<"end"<<std::endl;

	ret |= fix_vars_to_trace(home,last_assigned+1,n,traceA,traceB);
	
    }
    //print_vars();
    return ret;
}

Gecode::ExecStatus
AlignmentScore::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
    
    if (debug_out) {
	std::cout << "AlignmentScore::propagate " <<std::endl;
    }

    if (debug_out) {
	std::cout << "Begin propagation"<<std::endl;
	print_vars();
    }

    // ----------------------------------------
    // Propagation strategy:
    //
    // This is the naive implementation of the propagator
    // with only some efficiency improvements.
    // 
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
    
    const size_t n=seqA.length();
    const size_t m=seqB.length();
    
    
    ////////////////////////////////////////////////////////////
    // determine the undecided arc matches
    // 1.) arc matches that are not valid can be removed
    // 2.) arc matches where at least one end is certainly
    //     matched can be removed. Their score is added to the 
    //     unassigned end or, if both ends are matched, 
    //     arbitrarily to the left end

    
    const BasePairs &bpsA = arc_matches.get_base_pairsA();
    const BasePairs &bpsB = arc_matches.get_base_pairsB();

    // boolean Matrix giving for every arc match if we still consider
    // the arc match.
    Matrix<bool> considered_ams;
    considered_ams.resize(bpsA.num_bps(), bpsB.num_bps());
    

    Matrix<score_t> match_scores(n+1,m+1);
    // initialize match_score matrix
    for(size_type i=1; i<=n; i++) {
	for(size_type j=MD[i].min(); j<=(size_t)MD[i].max(); j++) {
	    match_scores(i,j)=2*scoring.basematch(i,j);
	}
    }
    
    
    for(size_t i=0; i<bpsA.num_bps(); ++i) {
	for(size_t j=0; j<bpsB.num_bps(); ++j) {
	    
	    const Arc &arcA=bpsA.arc(i);
	    const Arc &arcB=bpsB.arc(j);
	    
	    // check whether both ends can match
	    if (match_allowed(arcA.left(),arcB.left())
		&&
		match_allowed(arcA.right(),arcB.right())) {
		
		// if arc match allowed but none of the ends is forced
		// then the arc match is considered
		considered_ams.set(i,j,true); 
		
		// if one end is forced to match, then add twice the arc match score 
		// to the match score of the other end
		// Note: if both ends are forced, we add arc match scores arbitrarily to both ends
		//
		if (match_forced(arcA.right(),arcB.right())
		    &&
		    match_forced(arcA.left(),arcB.left())) {
		    match_scores(arcA.left(),arcB.left()) += scoring.arcmatch(arcA,arcB);
		    match_scores(arcA.right(),arcB.right()) += scoring.arcmatch(arcA,arcB);
		    considered_ams.set(i,j,false); // don't consider arc match  anymore
		} else if (match_forced(arcA.right(),arcB.right())) {
		    match_scores(arcA.left(),arcB.left()) += 2 * scoring.arcmatch(arcA,arcB);
		    considered_ams.set(i,j,false); // don't consider arc match  anymore
		} else if (match_forced(arcA.left(),arcB.left())) {
		    match_scores(arcA.right(),arcB.right()) += 2 * scoring.arcmatch(arcA,arcB);
		    considered_ams.set(i,j,false); // don't consider arc match  anymore
		}		
	    } else {
		// invalid arc match
		considered_ams.set(i,j,false);
	    } 
	}
    }
    // end of determining considered arcs and transferring scores
    ////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////
    // precompute upper bounds and determine tightness
    //
    Matrix<bool> tight;
    tight.resize(n+1,m+1);
    
    // precompute upper bounds for single matches i~j
    // for use in forward and backward algorithm
    Matrix<score_t> UBM(n+1,m+1);
    for (size_t i=1; i<=n; i++) {
	for (size_t j=MD[i].min(); j<=(size_t)MD[i].max(); j++) {
	    
	    bool tight_ij=true;
	    UBM(i,j)=ub_match(i,j,considered_ams,match_scores,&tight_ij);
	    tight.set(i,j,tight_ij);
	    
	}
	//std::cout << i << ": " << tight_i << std::endl;
    }
    ////////////////////////////////////////////////////////////



    // -------------------- RUN FORWARD ALGORITHM
    Matrix<infty_score_t> Fwd(n+1,m+1); // ForWarD matrix
    
    forward_algorithm(home,Fwd,UBM);
    
    //if (debug_out) {std::cout << "Fwd"<<std::endl<<Fwd <<std::endl;}
    
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
    backtrace_forward(home,Fwd,UBM,traceA,traceB);
    
    // test whether for all vars in a run the bound is tight
    // the variables in a "tight run" can be fixed to trace
    
    
    // -------------------- FIX TIGHT RUNS
    ret |= fix_tight_runs(home,traceA,traceB,tight);
    
    // -------------------- Lower bounding

    // evaluate alignment for obtaining lower bound
    score_t trace_score = evaluate_trace(traceA,traceB,considered_ams,match_scores);
    
    
    if (debug_out) {
	// print trace (DEBUGGING)
	std::cout << "TRACE"<<std::endl;
	std::cout << "Score: " << trace_score << std::endl;
	for (size_type i=1; i<=n; ++i ) {
	    std::cout << i<<"~"<<traceA[i] << " ";
	}
	std::cout << std::endl;
	for (size_type j=1; j<=m; j++) {
	    std::cout << traceB[j]<<"~"<<j << " ";
	}
	std::cout << std::endl;
    }


    // constrain Score by trace_score, which is a lower bound
    
    if (Gecode::me_failed(Score.gq(home,(int)trace_score))) {
     	if (debug_out) std::cout << "Setting of new lower bound failed."<<std::endl;
      	return Gecode::ES_FAILED;
    }
    
    


    // -------------------- RUN BACKWARD ALGORITHM
    Matrix<infty_score_t> Bwd(n+1,m+1); 

    backward_algorithm(home,Bwd,UBM);
    //if (debug_out) {std::cout << "Bwd"<<std::endl<<Bwd <<std::endl;}


    if (false && debug_out) {
	// ------------------------------------------------------------
	// DEBUGGING: check trace
	for (size_t i=1; i<=n; i++) {
	    std::cout << "TRACE " << i <<" " << traceA[i] <<" ";
	    if (traceA[i]!=0) {
		size_t j = traceA[i];
		score_t matchscore=evaluate_tracematch(traceA,traceB,considered_ams,match_scores,i,j);
		
		cout
		    << match_allowed(i,j)<<" "
		    <<Fwd(i-1,j-1)<<"+"<< matchscore<<"["<<ub_match(i,j,considered_ams,match_scores)<<"]"<<"+"<<Bwd(i,j)<<" = "
		    <<Fwd(i-1,j-1)+matchscore+Bwd(i,j);
	    }
	    std::cout << std::endl;
	}
    }
    
    
    // -------------------- PRUNE VARIABLES
    ret |= prune(home,Fwd,Bwd,UBM);
        
    
    if (debug_out) {
	std::cout << "After Pruning:"<<std::endl;
	print_vars();
    }
    
    if (Gecode::me_failed(ret)) {
	if (debug_out) std::cout << "Fail after pruning."<<std::endl;
	return Gecode::ES_FAILED;
    }
    
    // test whether all vars fixed, then subsume (can we subsume earlier?)
    if ( all_vars_fixed() ) {
	return ES_SUBSUMED(*this,dispose(home));
    }

    // -------------------- select CHOICE for the space
    if (!Gecode::me_modified(ret)) {
	// don't call choice if propagate will be called again anyway (due to modifications)
	choice(static_cast<RNAalignment&>(home),Fwd,Bwd,traceA,traceB,UBM);
    }

    return Gecode::me_modified(ret)?Gecode::ES_NOFIX:Gecode::ES_FIX;
}



////////////////////////////////////////////////////////////
// OBSOLETE CODE

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

