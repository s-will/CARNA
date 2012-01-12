#include "alignment_score.hh"
#include "LocARNA/arc_matches.hh"

//#include "WinDisplay.hh"
#include "RNAalignment.hh"

#include <limits>

#include <assert.h>

using namespace LocARNA;
using namespace std;

const bool debug_out=false;
// const bool debug_out=true;


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
				      BoolVec &forcedA, 
				      BoolVec &forcedB,
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

	    if ( (!forced) && match_arrow_allowed(k,l) ) { // if the match
		// is allowed, but not forced
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
				 const ScoreMatrix &considered_ams,
				 bool right) const {
    
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

	    size_t k=(right)?arcA->left():arcA->right();
	    size_t l=(right)?arcB->left():arcB->right();
	    
	    
	    infty_score_t entry=infty_score_t::neg_infty;
	    if ( match_arrow_allowed(k,l) && considered_ams.get(arcA->idx(), arcB->idx())!=0 ) {
		entry = MAT_match + considered_ams.get(arcA->idx(), arcB->idx());
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
    
    // RETURN bound
    // maximal bound for arcs matches to the left is in MAT[adjlB.size()] 
    //std::cout <<adjlA.size()<<" "<<adjlB.size()<<" "<< MAT[adjlB.size()]<<std::endl;
    return bound;
}


score_t
AlignmentScore::ub_match(size_type i, size_type j,
			 const ScoreMatrix &considered_ams,
			 const ScoreMatrix &match_scores,
			 bool *tight_ptr) const {
    //std::cout << "AlignmentScore::ub_match "<<i<< " "<<j<<std::endl;
  
    // compute upper bound for the contribution of matching positions i
    // and j of respective sequences A and B
    // by summing over all possible base pair matchs
    //
    // when selecting a certain structure: maximize over possible base pair matchs
    
    score_t bound=0;
    
    //std::cout << bound<<std::endl;
        
    bound += 
	bound_arcmatches(arc_matches.get_base_pairsA().right_adjlist(i),
			 arc_matches.get_base_pairsB().right_adjlist(j),
			 considered_ams,
			 true);
    
    bound += 
	bound_arcmatches(arc_matches.get_base_pairsA().left_adjlist(i),
			 arc_matches.get_base_pairsB().left_adjlist(j),
			 considered_ams,
			 false);
    
    // when is the bound TIGHT? A: if there are no considered arc matches
    // contributing to the bound.
    // Since arc matches with negative score do never contribute,
    // the bound is tight iff the bound is 0!!!
    bool tight=(bound==0);

    // return tightness
    if (tight_ptr!=NULL && *tight_ptr==true) {
	*tight_ptr = tight;
    }
    
    bound += match_scores(i,j);

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
AlignmentScore::evaluate_tracematch(const SizeVec &traceA,
				    const SizeVec &traceB,
				    const ScoreMatrix &considered_ams,
				    const ScoreMatrix &match_scores,
				    size_type i,
				    size_type j) const{
    
    score_t matchscore=match_scores(i,j);

    // for all pairs of arcs in A and B that have right ends i and j, respectively
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_right_end_list(i,j).begin();
	arc_matches.common_right_end_list(i,j).end() != it; ++it ) {
	
	const ArcMatch &am = arc_matches.arcmatch(*it);

	
	if ( traceA[am.arcA().left()] == am.arcB().left() ) { 
	    // if the left ends of arcs arcA and arcB match 
	    
	    score_t amscore=considered_ams.get(am.arcA().idx(),am.arcB().idx());
	    matchscore += amscore;
	}
    }

    // same for common left ends
    for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(i,j).begin();
	arc_matches.common_left_end_list(i,j).end() != it; ++it ) {	
	
	const ArcMatch &am = arc_matches.arcmatch(*it);
	
	if ( traceA[am.arcA().right()] == am.arcB().right() ) { 
	    // if the right ends of arcs arcA and arcB  match 

	    score_t amscore=considered_ams.get(am.arcA().idx(),am.arcB().idx());
	    matchscore += amscore;
	}
    }

    return matchscore;
}

score_t
AlignmentScore::evaluate_trace(const SizeVec &traceA,
			       const SizeVec &traceB,
			       const ScoreMatrix &considered_ams,
			       const ScoreMatrix &match_scores
			       ) const {
    
    const size_t n=seqA.length();
    const size_t m=seqB.length();

    score_t score=0;

    size_type i=1;
    size_type j=1;
    
    bool gap_openA=false; // is a gap open in sequence A?
    bool gap_openB=false; // is a gap open in sequence B?

    // while not all positions i and j consumed
    while (i<=n || j<=m) {
	// invariant: score is the score for aligning A_1..i-1 to B_1..j-1 according to the trace 
	if (i<=n && traceA[i]==0) { // deletion of i after j-1
	    if (!gap_openA) score+=2*scoring.indel_opening();
	    gap_openA=true;
	    gap_openB=false;
	    score += 2*scoring.gapA(i,j-1);
	    ++i;
	} else if (j<=m && traceB[j]==0) { // insertion of j after i-1
	    if (!gap_openB) score+=2*scoring.indel_opening();
	    gap_openB=true;
	    gap_openA=false;
	    score += 2*scoring.gapB(i-1,j);
	    ++j;
	} else { // match between i and j
	    assert(j==traceA[i]);
	    assert(i==traceB[j]);
	    
	    gap_openA=false;
	    gap_openB=false;

	    score += evaluate_tracematch(traceA,traceB,considered_ams,match_scores,i,j);
	    ++i;
	    ++j;
	}
    }
    
    return score;
}

void
AlignmentScore::print_vars(std::ostream &out) const {
    out << "Matches/Deletions:    ";
    for (size_t i=0; i<=seqA.length(); i++) {
	out <<i<<(M[i].assigned()?(M[i].val()==0?"g":"~"):"?")<<MD[i]<<", "; 
    }
    
    //out << "Matches/Deletions:    " << MD << std::endl;
    //out << "Match Flags:          " << M << std::endl;
    out << "Score:      " << Score << std::endl;
}



void
AlignmentScore::
forward_algorithm(Gecode::Space& home,
			 InftyScoreRRMatrix &Fwd,
			 InftyScoreRRMatrix &FwdA,
			 InftyScoreRRMatrix &FwdB,
			 const ScoreMatrix &UBM
			 ) {

    const size_t n=seqA.length();
    const size_t m=seqB.length();

    // ----------------------------------------
    // Forward algorithm
    //
    
    // Definition( Fwd-matrix )
    // for j \in MD[i]:
    // Fwd(i,j)  := max score of an alignment R[1..i], S[1..j]
    // for j \in MD[i]
    // FwdA(i,j) := max score of an alignment R[1..i], S[1..j] where R[i] is deleted
    // for j \in MD[i]
    // FwdB(i,j) := max score of an alignment R[1..i], S[1..j] where S[j] is inserted
    
    
    ////////////////////
    // initialize
    Fwd(0,0)=(infty_score_t)0;
    FwdA(0,0) = (infty_score_t) (2*scoring.indel_opening()); // infty_score_t::neg_infty;
    FwdB(0,0) = (infty_score_t) (2*scoring.indel_opening()); // infty_score_t::neg_infty;
    for(size_type i=1; i<=n; i++) {
	if (!trace_allowed(i,0)) continue; // assuming monotonicity:
					   // one could as well break
	
	if (deletion_arrow_allowed(i,0)) {
	    FwdA(i,0) = FwdA(i-1,0) + 2*scoring.gapA(i,0);
	} else {
	    FwdA(i,0) = infty_score_t::neg_infty;
	}
	FwdB(i,0)=infty_score_t::neg_infty;
	Fwd(i,0)  = FwdA(i,0);
    }
    
    for(size_type j=1; j<=max_col(0); j++) {
	if (!trace_allowed(0,j)) continue; // assuming monotonicity:
					   // one could as well break
	
	FwdA(0,j)=infty_score_t::neg_infty;
	FwdB(0,j) = FwdB(0,j-1)  + 2*scoring.gapB(0,j);
	Fwd(0,j) = FwdB(0,j);
    }
    
    Fwd(n,m) = infty_score_t::neg_infty;
    // end init
    ////////////////////
    // recurse
    for(size_type i=1; i<=n; i++) {
	size_t minj = std::max((size_t)1,min_col(i));
	size_t maxj = max_col(i);
		
	for(size_type j=minj; j<=maxj; j++) {
	    //note: trace_allowed(i,j) holds (by j loop range)

	    Fwd(i,j) = infty_score_t::neg_infty;
	    FwdA(i,j) = infty_score_t::neg_infty;
	    FwdB(i,j) = infty_score_t::neg_infty;
	    if (match_arrow_allowed(i,j)) {
		Fwd(i,j) = Fwd(i-1,j-1) + UBM(i,j);
	    }
	    if (deletion_arrow_allowed(i,j)) {
		FwdA(i,j) = std::max(FwdA(i,j),FwdA(i-1,j)+2*scoring.gapA(i,j));
		FwdA(i,j) = std::max(FwdA(i,j),Fwd(i-1,j) +2*scoring.gapA(i,j)+2*scoring.indel_opening());
		
		Fwd(i,j)  = std::max(Fwd(i,j),FwdA(i,j));
	    }
	    if (insertion_arrow_allowed(i,j)) {
		FwdB(i,j) = std::max(FwdB(i,j),FwdB(i,j-1)+2*scoring.gapB(i,j));
		FwdB(i,j) = std::max(FwdB(i,j),Fwd(i,j-1) +2*scoring.gapB(i,j)+2*scoring.indel_opening());
		
		Fwd(i,j) = std::max(Fwd(i,j),FwdB(i,j));
	    }
	}
    }
}

void
AlignmentScore::
backtrace_forward(Gecode::Space &home, 
		  const InftyScoreRRMatrix &Fwd,
		  const InftyScoreRRMatrix &FwdA,
		  const InftyScoreRRMatrix &FwdB,
		  const ScoreMatrix &UBM,
		  SizeVec &traceA,
		  SizeVec &traceB
		  ) {
    const int n=seqA.length();
    const int m=seqB.length();

    // ----------------------------------------
    // do backtracking of the Fwd matrix in order to obtain a lower
    // bound of alignment score and guide the heuristic
    
    traceA.resize(n+1);
    traceB.resize(m+1);
    
    {
	enum { FWD, FWD_A, FWD_B } state;
	state = FWD;
	int i=n;
	int j=m;
	
	while (i>0 && j>0) {
	    // std::cerr << i << " " << j << " " << state << std::endl;
	    switch(state) {
	    case FWD:
		if (match_arrow_allowed(i,j) && Fwd(i,j) == Fwd(i-1,j-1) + UBM(i,j) ) {
		    traceA[i]=j;
		    traceB[j]=i;
		    --i;
		    --j;
		} else if ( deletion_arrow_allowed(i,j) && Fwd(i,j)==FwdA(i,j) ) {
		    state = FWD_A;
		} else if ( insertion_arrow_allowed(i,j) && Fwd(i,j)==FwdB(i,j) ) {
		    state = FWD_B;
		} else {
		    std::cerr << "Traceback error FWD\n"<<std::endl;
		    exit(-1);
		}
		break;
	    case FWD_A:
		if ( FwdA(i,j) == FwdA(i-1,j) + 2*scoring.gapA(i,j) ) {
		    traceA[i]=0;
		    --i;
		} else if ( FwdA(i,j) == Fwd(i-1,j) + 2*scoring.gapA(i,j) + 2*scoring.indel_opening() ) {
		    traceA[i]=0;
		    --i;
		    state=FWD;
		} else {
		    std::cerr << "Traceback error FWD_A\n"<<std::endl;
		    exit(-1);
		}
		break;
	    case FWD_B:
		if ( FwdB(i,j) == FwdB(i,j-1) + 2*scoring.gapB(i,j) ) {
		    traceB[j]=0;
		    --j;
		} else if ( FwdB(i,j) == Fwd(i,j-1) + 2*scoring.gapB(i,j) + 2*scoring.indel_opening() ) {
		    traceB[j]=0;
		    --j;
		    state=FWD;
		} else {
		    std::cerr << "Traceback error FWD_B\n"<<std::endl;
		    exit(-1);
		}
		break;
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
AlignmentScore::backward_algorithm(Gecode::Space& home, 
				   InftyScoreRRMatrix &Bwd,
				   InftyScoreRRMatrix &BwdA,
				   InftyScoreRRMatrix &BwdB,
				   const ScoreMatrix &UBM) {
    const size_t n=seqA.length();
    const size_t m=seqB.length();

    // ----------------------------------------
    // Backward algorithm
    //
    
    // Definition( Fwd-matrix )
    // for j \in MD[i]:
    // Bwd(i,j)  := max score of an alignment R[i+1..n], S[j+1..m]
    // for j \in MD[i]
    // BwdA(i,j) := max score of an alignment R[i+1..n], S[j+1..m] where R[i+1] is deleted
    // for j \in MD[i]
    // BwdB(i,j) := max score of an alignment R[i+1..n], S[j+1..m] where S[j+1] is inserted
    
    // initialize
    Bwd(n,m)=(infty_score_t)0;
    BwdA(n,m) = (infty_score_t)(2*scoring.indel_opening()); // infty_score_t::neg_infty;
    BwdB(n,m) = (infty_score_t)(2*scoring.indel_opening()); // infty_score_t::neg_infty;
  
    for(size_type i=n; i>0;) { // for i=n-1 downto 0
	--i;

	if (!trace_allowed(i,m)) continue; // assuming monotonicity:
					   // one could as well break

	if (deletion_arrow_allowed(i+1,m)) {
	    BwdA(i,m) = BwdA(i+1,m) + 2*scoring.gapA(i+1,m);
	} else {
	    Bwd(i,m) = infty_score_t::neg_infty;
	}
	BwdB(i,m)=infty_score_t::neg_infty;
	Bwd(i,m)=BwdA(i,m);
    }
    for(size_type j=max_col(n); j>min_col(n);) { // for j=m-1 downto 0
	--j;

	if (!trace_allowed(n,j)) continue; // assuming monotonicity:
					   // one could as well break

	if (insertion_arrow_allowed(n,j+1)) {
	    BwdB(n,j) = BwdB(n,j+1) + 2*scoring.gapB(n,j+1);
	} else {
	    BwdB(n,j) = infty_score_t::neg_infty;
	}
	BwdA(n,j)=infty_score_t::neg_infty;
	Bwd(n,j) = BwdB(n,j);
    }
    
    Bwd(0,0) = infty_score_t::neg_infty;   

    // recurse
    for(size_type i=n; i>0;) { // for i=n-1 downto 0
	--i;
	
	size_t minj = min_col(i);
	size_t maxj = std::min(max_col(i),m-1);
		
	for(size_type j=maxj+1; j>minj;) { // for j=maxj downto minj
	    --j;
	    
	    //note: trace_allowed(i,j) holds (by j loop range)

	    Bwd(i,j) = infty_score_t::neg_infty;
	    BwdA(i,j) = infty_score_t::neg_infty;
	    BwdB(i,j) = infty_score_t::neg_infty;
	    	    
	    if (match_arrow_allowed(i+1,j+1)) {   // match i+1 ~ j+1 (!)
		Bwd(i,j) = Bwd(i+1,j+1) + UBM(i+1,j+1);
	    }
	    if (deletion_arrow_allowed(i+1,j)) { // delete i+1 after j
		BwdA(i,j) = std::max(BwdA(i,j),BwdA(i+1,j)+2*scoring.gapA(i+1,j));
		BwdA(i,j) = std::max(BwdA(i,j),Bwd(i+1,j) +2*scoring.gapA(i+1,j)+2*scoring.indel_opening());

		Bwd(i,j)=std::max(Bwd(i,j),BwdA(i,j));
	    }
	    if (insertion_arrow_allowed(i,j+1)) { // insert j+1 after i
		BwdB(i,j) = std::max(BwdB(i,j),BwdB(i,j+1)+2*scoring.gapB(i,j+1));
		BwdB(i,j) = std::max(BwdB(i,j),Bwd(i,j+1) +2*scoring.gapB(i,j+1)+2*scoring.indel_opening());

		Bwd(i,j)=std::max(Bwd(i,j),BwdB(i,j));
	    }
	}
    }
}

Gecode::ModEvent
AlignmentScore::prune(Gecode::Space& home, 
		      const InftyScoreRRMatrix &Fwd,
		      const InftyScoreRRMatrix &FwdA,
		      //const InftyScoreRRMatrix &FwdB,
		      const InftyScoreRRMatrix &Bwd,
		      const InftyScoreRRMatrix &BwdA,
		      //const InftyScoreRRMatrix &BwdB,
		      const ScoreMatrix &UBM
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
	size_t minj = min_col(i);
	size_t maxj = max_col(i);
		
	for(size_type j=minj; j<=maxj; j++) {
	    // prune MD variable for match i~j or deletion of i after j
	    if (match_or_deletion_allowed(i,j)) {
		infty_score_t ub = Fwd(i,j)+Bwd(i,j);
		
		//here, we are not interested in insertions!
		//if (j>0 && insertion_arrow_allowed(i,j) && j<m && insertion_arrow_allowed(i,j+1)) {
		//    ub = std::max(ub, FwdB(i,j)+BwdB(i,j)-2*scoring.indel_opening()); 
		//}
		
		if (i>0 && deletion_arrow_allowed(i,j) && i<n && deletion_arrow_allowed(i+1,j)) {
		    ub = std::max(ub, FwdA(i,j)+BwdA(i,j)-2*scoring.indel_opening()); 
		}
		
		if ( (!ub.is_finite()) || (ub < (infty_score_t) Score.min())) {
		    ret |= MD[i].nq(home,(int)j);
		    //std::cerr << "MD["<<i<<"].nq(home,(int)"<<j<<") "<< ret << " " << Fwd(i,j)+Bwd(i,j) << std::endl;
		}
	    }
	}
    }
    
    ////////////////////////////////////////
    // prune on M

    // prune impossible deletions
    for(size_type i=1; i<=n; i++) {
	bool del=false; // is deletion possible?
	size_t minj = min_col(i);
	size_t maxj = max_col(i);
	for(size_type j=minj; !del && j<=maxj; j++) {
	    
	    // upper bound for deletion of i after j
	    infty_score_t ubd = FwdA(i,j)+Bwd(i,j);
	    
	    if (i>0 && deletion_arrow_allowed(i,j) && i<n && deletion_arrow_allowed(i+1,j)) {
		ubd = std::max(ubd,FwdA(i,j)+BwdA(i,j)-2*scoring.indel_opening());
	    }
	    
	    if ( ubd.is_finite() && (ubd.finite_value() >= Score.min())) {
		del=true;
	    }
	}
	if (!del) {
	    ret |= M[i].nq(home,0); // deletion of i not possible
	    //std::cerr << "M["<<i<<"].nq(home,0) "<<ret<<std::endl;
	}
    }
    
    // prune impossible matches
    for(size_type i=1; i<=n; i++) {
	bool match=false; // is match possible?
	size_t minj = std::max((size_t)1,min_col(i));
	size_t maxj = max_col(i);
	for(size_type j=minj; !match && j<=maxj; j++) {
	    
	    if (match_arrow_allowed(i,j)) {
		// upper bound for match i~j (Note: we don't need the matrices
		// FwdA,FwdB,BwdA,BwdB here, since i~j and Bwd(i,j) is general)
		infty_score_t ubm = Fwd(i-1,j-1)+UBM(i,j)+Bwd(i,j);
		
		if ( ubm.is_finite() && (ubm.finite_value() >= Score.min())) {
		    match = true; // match is possible
		}
	    }
	}
	if (!match) {
	    ret |= M[i].nq(home,1); // match of i not possible
	    // std::cerr << "M["<<i<<"].nq(home,1) "<<ret<<std::endl;
	}
    }
    

    return ret;
}


void
AlignmentScore::choice(RNAalignment &s,
		       const InftyScoreRRMatrix &Fwd,
		       const InftyScoreRRMatrix &Bwd,
		       const SizeVec &traceA,
		       const SizeVec &traceB,
		       score_t trace_score,
		       const ScoreMatrix &UBM,
		       const ScoreMatrix &match_scores) const {
    
    // We choose a variable and domain change in four steps
    // STEP 0: decide whether the trace from the forward matrix is a new
    //         lower bound.
    //         Then, return; the brancher will enumerate the new solution.
    // STEP 1: decide whether to enumerate M or MD variable. 
    //         In case, select M var
    // STEP 2: select MD variable MD[pos]
    // STEP 3: select domain change of MD[pos]
    
    // NOTE: even more severe than in prune, UBM does not contain the
    // best upper bounds anymore to save time, we use UBM nevertheless
    // instead of re-computing the bounds
    
    const size_t n=seqA.length();
    const size_t m=seqB.length();
    
    // ------------------------------------------------------------
    // STEP 0: check whether the trace is a new solution
    //
    if (s.choice_data.new_lower_bound) {
	return; //choice done
    }
    
    // ------------------------------------------------------------
    // STEP 1: check whether to enumerate M variables
    //
    
    //strategy: whenever a MD[i] variable is assigned, but M[i] is not
    // enumerate M[i]
    for (size_t i=0; i<=n; i++) {
	if (MD[i].assigned() && !M[i].assigned()) {
	    s.choice_data.pos=i;
	    s.choice_data.enum_M=true;
	    return; // choice done
	}
    }
    
    s.choice_data.enum_M=false;
    
    // ------------------------------------------------------------
    // STEP 2: select pos / variable MD[pos] 
    //

    // TODO: revise the variable selection strategy
    
    // determine position with largest weight
    
    if (debug_out) std::cout <<"Determine choice"<<std::endl;
    //if (debug_out) print_vars(std::cout);
    
    vector<score_t> weights;
    score_t maxweight=numeric_limits<score_t>::min();
    
    size_t total_size=0; // count total domain size
    
    weights.resize(n+1);
    for (size_t i=0; i<=n; i++) {
	
	total_size += MD[i].size();
	
	#define MAXIMIZE
	// weight by maximal base pair bound
#ifdef MAXIMIZE 
	weights[i]=numeric_limits<score_t>::min();
#else
	weights[i]=0;
#endif
	
	if (!s.MD[i].assigned()) {
	    for (size_t j=s.MD[i].min(); j<=(size_t)s.MD[i].max(); j++) {
		if (s.MD[i].in((int)j)) {
#ifdef MAXIMIZE 
		    weights[i]=max(weights[i],
				   UBM(i,j)-match_scores(i,j) //ub_match(i,j,false) //CHECK: -2*scoring.base_match(i,j) or -match_scores(i,j)
				   +
				   (score_t)s.MD[i].size());
#else
		    weights[i] +=
			UBM(i,j)-match_scores(i,j) //ub_match(i,j,false) //CHECK: -2*scoring.base_match(i,j) or -match_scores(i,j)
			+
			1;
#endif
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
    
    // select mid position
    pos=best_left_end + (best_run_len/2); // split the longest run
    

    assert(0<=pos && pos<=n);
    
    
    // ------------------------------------------------------------
    // STEP 3: select domain change of MD[pos] 
    //
    
    // step 3a) determine value val for MD[pos] that is compatible with
    // trace. If pos is matched val is traceA[pos]. If pos is deleted by
    // trace, we need to determine val, such that pos is deleted between
    // val and val+1
    //
    size_t val = (pos==0)?0:traceA[pos];
    //
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
	
	if (debug_out) print_vars(std::cout);
    }
    
    
    if (debug_out) std::cout << "Choose pos:"<<pos<<" val:"<<val<<std::endl;
    
    // step 3b) determine values for splitting choice
    
    // determine minimal and maximal bound accross all j in MD[pos]
    //
    score_t minb=numeric_limits<score_t>::max();
    score_t maxb=numeric_limits<score_t>::min();

    if (debug_out) std::cout << "CHOICE AT "<< pos <<std::endl;
    for (size_t j=s.MD[pos].min(); j<=(size_t)s.MD[pos].max(); j++) {
	infty_score_t b;
	if (s.MD[pos].in(j) && (b=Fwd(pos,j)+Bwd(pos,j)).is_finite()) {
	    minb=min(minb,b.finite_value());
	    maxb=max(maxb,b.finite_value());
	}
    }
    
    // use strategy with fixed percentage !?
    // float fraction=total_size/(float)(n*m);
    int percentage=80; // max((float)0.5,1-sqrt(fraction));

    // determine the j in MD[pos], where the bound is greater or equal
    // than minb+percentage*(maxb-minb)
    
    s.choice_data.values.resize(0);
    s.choice_data.values.reserve(m/4);
    
    for (size_t j=s.MD[pos].min(); j<=(size_t)s.MD[pos].max(); j++) {
	infty_score_t b;
	if (s.MD[pos].in(j) ) {
	    if (debug_out) 
		std::cout << j << ":" << b<<" ";
	    if ((b=Fwd(pos,j)+Bwd(pos,j)).is_finite() && b.finite_value()>=minb+(percentage*(maxb-minb))/100) { 
		s.choice_data.values.push_back(j);
	    } 
	}
    }
    
    // make sure that the domain gets smaller due to the above choice. If it does not
    // split the domain.
    //
    if ( s.MD[pos].size()==s.choice_data.values.size() ) {
	
	size_t minval=s.MD[pos].min();
	size_t maxval=s.MD[pos].max();
    	
	size_t medval = minval + (maxval-minval)/2 ;
		
	// choose according to trace.  make sure that we always
	// include the j in the trace for the left branch. Recall: val
	// is the value of MD[pos] that is compatible to trace.
	if (val<=medval) {
	    maxval=medval;
	} else {
	    minval=medval+1;
	}
	
	// set values accordingly
	s.choice_data.values.resize(0);
	for (size_t j=minval; j<=maxval; j++) {
	    s.choice_data.values.push_back(j);
	}
    }
    
    if (debug_out) { 
	std::cout << s.MD[pos] << " IN:";
	for (size_t i=0; i< s.choice_data.values.size(); ++i) {
	    std::cout << " " << s.choice_data.values[i];
	}
	std::cout << std::endl;
    }
    
    // copy choice to the space
    s.choice_data.pos=pos;
    s.choice_data.val=val;    
}


Gecode::ModEvent
AlignmentScore::fix_vars_to_trace(Gecode::Space &home,
				  size_t start,
				  size_t end,
				  const SizeVec &traceA,
				  const SizeVec &traceB
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
			       const SizeVec &traceA,
			       const SizeVec &traceB,
			       const BoolMatrix &tight) {
    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
    
    size_t n=seqA.length();
    
    // determine tight runs
    
    int last_assigned=0;
    bool run_is_tight=true;
    
    for (size_t i=1; i<=n; i++) {
	if (MD[i].assigned() && M[i].assigned()) {// end of a run
	    
	    // if the whole run is tight
	    if (i>(size_t)(last_assigned+1) && run_is_tight) {
		
		// we can fix the whole run
		/*
		print_vars(std::cout);
		for (size_t k=0; k<=n; k++) {
		    if (traceA[k]!=0) {
			bool t=tight.get(k,traceA[k]);
			std::cout <<k<<"~"<<traceA[k]<<":"<<t<<" ";
		    }
		}
		std::cout << std::endl;
		*/
		//if (debug_out) std::cerr << "fix "<<(last_assigned+1)<<"-"<<(i-1)<<std::endl;
		ret |= fix_vars_to_trace(home,last_assigned+1,i-1,traceA,traceB);
	       
	    }
	    
	    last_assigned=i;
	    run_is_tight=true;
	    
	} else { // within a run
	    // the run is tight at i if only if all possible matches
	    // of i with some value val are tight
	    if (M[i].in(1)) { // if match is possible
		for (int val=MD[i].min(); run_is_tight && val<=MD[i].max(); val++) {
		    run_is_tight &= (!MD[i].in(val)) || tight.get(i,val);
		}
	    }
	}
    }
    
    
    // if the whole run from last_assigned to the end is tight
    if ((size_t)(last_assigned+1)<n+1 && run_is_tight) {
	// we can fix the whole run
	
	if (debug_out) std::cerr << "fix "<<(last_assigned+1)<<"-"<<n<<std::endl;
			
	ret |= fix_vars_to_trace(home,last_assigned+1,n,traceA,traceB);
	
    }
    //print_vars(std::cout);
    return ret;
}

void
AlignmentScore::prune_decided_arc_matches(ScoreMatrix &considered_ams, ScoreMatrix &match_scores) {
    const size_t n=seqA.length();
    const size_t m=seqB.length();

    const BasePairs &bpsA = arc_matches.get_base_pairsA();
    const BasePairs &bpsB = arc_matches.get_base_pairsB();    

    considered_ams.resize(bpsA.num_bps(), bpsB.num_bps());
    
    match_scores.resize(n+1,m+1);
    
    // initialize match_score matrix
    for(size_type i=1; i<=n; i++) {
	for(size_type j=std::max(1,MD[i].min()); j<=(size_t)MD[i].max(); j++) {
	    match_scores(i,j)=2*scoring.basematch(i,j);
	}
    }
    
    
    for(size_t i=1; i<=n; ++i) {
	for(size_t j=std::max(1,MD[i].min()); j<=(size_t)MD[i].max(); ++j) {	    
	    if (match_arrow_allowed(i,j)) {
		
		bool forced_ij=match_forced(i,j);
		
		// iterate through all arc matches to the right
		// if match_forced(i,j) transfer arcmatch score to the right end
		// otherwise if other right end is forced transfer arcmatch score to match i,j
		
		for(ArcMatchIdxVec::const_iterator it=arc_matches.common_left_end_list(i,j).begin();
		    arc_matches.common_left_end_list(i,j).end() != it; ++it ) {
		    
		    const ArcMatch &am = arc_matches.arcmatch(*it);
		    
		    const Arc &arcA = am.arcA();
		    const Arc &arcB = am.arcB();
		    
		    size_t k=arcA.right();
		    size_t l=arcB.right();
		    
		    if (!match_arrow_allowed(k,l)) continue;
		    
		    bool forced_kl=match_forced(k,l);
		    
		    if ( forced_ij && forced_kl ) {
			match_scores(i,j) += scoring.arcmatch(am);
			match_scores(k,l) += scoring.arcmatch(am);
			considered_ams.set(arcA.idx(),arcB.idx(),0); // don't consider arc match anymore
		    }
		    else if ( forced_ij && !forced_kl ) {
			match_scores(k,l) += 2 * scoring.arcmatch(am);
			considered_ams.set(arcA.idx(),arcB.idx(),0); // don't consider arc match anymore
		    }
		    else if ( !forced_ij && forced_kl ) {
			match_scores(i,j) += 2 * scoring.arcmatch(am);
			considered_ams.set(arcA.idx(),arcB.idx(),0); // don't consider arc match anymore
		    } else {
			considered_ams.set(arcA.idx(),arcB.idx(),scoring.arcmatch(am)); // consider arc match
		    }
		}
	    }
	    // else: for !allowed_match(i,j) nothing is to do
	}
    }
    // end of determining considered arcs and transferring scores
    ////////////////////////////////////////////////////////////

}

void
AlignmentScore::print_trace(std::ostream &out, const SizeVec &traceA,
			    const SizeVec &traceB, score_t trace_score) const {
    const size_t n=seqA.length();
    const size_t m=seqB.length();

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

Gecode::ExecStatus
AlignmentScore::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
    
    if (debug_out) {
	std::cerr << "AlignmentScore::propagate " <<std::endl;
    }

    if (debug_out) {
	std::cout << "Begin propagation"<<std::endl;
	print_vars(std::cout);
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
    
    // Matrix. contains arcmatch score of every arc match with allowed left and right match if we still consider
    // the arc match. Contains neg_infty otherwise.
    // NOTE: entries where one match in the arc match is not allowed are undefined!!!
    ScoreMatrix considered_ams;

    // matrix for base match scores
    ScoreMatrix match_scores;
    
    prune_decided_arc_matches(considered_ams,match_scores);
    
    
    ////////////////////////////////////////////////////////////
    // precompute upper bounds and determine tightness
    //
    Matrix<bool> tight;
    tight.resize(n+1,m+1);
    
    // precompute upper bounds for single matches i~j
    // for use in forward and backward algorithm
    ScoreMatrix UBM(n+1,m+1);
    for (size_t i=1; i<=n; i++) {
	for (size_t j=MD[i].min(); j<=(size_t)MD[i].max(); j++) {
	    
	    bool tight_ij=true;
	    UBM(i,j)=ub_match(i,j,considered_ams,match_scores,&tight_ij);
	    tight.set(i,j,tight_ij);
	    
	}
	//std::cout << i << ": " << tight_i << std::endl;
    }
    ////////////////////////////////////////////////////////////

    // construct range vector for DP matrices
    std::vector<RowRangeMatrix<infty_score_t>::size_pair_type> ranges;
    for (size_t i=0; i<(size_t)MD.size(); i++) {
	ranges.push_back(RowRangeMatrix<infty_score_t>::size_pair_type(min_col(i),max_col(i)));
    }
    
    // -------------------- RUN FORWARD ALGORITHM
    //InftyScoreRRMatrix Fwd(n+1,m+1); // ForWarD matrix
    //InftyScoreRRMatrix FwdA(n+1,m+1); // ForWarD matrix, A_i gapped
    //InftyScoreRRMatrix FwdB(n+1,m+1); // ForWarD matrix, B_j gapped
    InftyScoreRRMatrix Fwd(ranges); // ForWarD matrix
    InftyScoreRRMatrix FwdA(ranges); // ForWarD matrix, A_i gapped
    InftyScoreRRMatrix FwdB(ranges); // ForWarD matrix, B_j gapped
    
    forward_algorithm(home,Fwd,FwdA,FwdB,UBM);
    
    //if (debug_out) {std::cout << "Fwd"<<std::endl<<Fwd <<std::endl;}
    
    // the forward algorithm yields an upper bound
    // ==> constrain the score variable
    if (Fwd(n,m).is_finite()) {
	if (Gecode::me_failed(Score.lq(home,(int)Fwd(n,m).finite_value()))) {
	    if (debug_out) std::cout << "Bounding after Fwd failed: "<<(int)Fwd(n,m).finite_value()<<std::endl;
	    return Gecode::ES_FAILED;
	}
    } else {
	if (debug_out) std::cout << "Fwd has infinite bound: fail."<<std::endl;
	return Gecode::ES_FAILED;
    }
    
    // trace vectors. idea one entry per position gives position in
    // the other sequence to that the pos is matched or 0 (in case of
    // gap).
    SizeVec traceA;
    SizeVec traceB;
    
    // -------------------- BACKTRACE FWD
    backtrace_forward(home,Fwd,FwdA,FwdB,UBM,traceA,traceB);
    
    // test whether for all vars in a run the bound is tight
    // the variables in a "tight run" can be fixed to trace
    
    
    // -------------------- FIX TIGHT RUNS
    ret |= fix_tight_runs(home,traceA,traceB,tight);
    
    // -------------------- Lower bounding

    // evaluate alignment for obtaining lower bound
    score_t trace_score = evaluate_trace(traceA,traceB,considered_ams,match_scores);
    
    
    if (debug_out) {
	print_trace(std::cout,traceA,traceB,trace_score);
    }

    // constrain Score by trace_score, which is a lower bound
    
    // test whether the trace can improve the lower bound
    RNAalignment &rahome = static_cast<RNAalignment&>(home);
    if (Score.min()<trace_score) {
	rahome.choice_data.new_lower_bound=true;
	rahome.choice_data.best_traceA=traceA;
	rahome.choice_data.best_traceB=traceB;
	rahome.choice_data.best_trace_score=trace_score;	
    }

    if (Gecode::me_failed(Score.gq(home,(int)trace_score))) {
     	if (debug_out) std::cout << "Setting of new lower bound failed."<<std::endl;
      	return Gecode::ES_FAILED;
    }
    
    
    // -------------------- RUN BACKWARD ALGORITHM
    //InftyScoreMatrix Bwd(n+1,m+1); 
    //InftyScoreMatrix BwdA(n+1,m+1); 
    //InftyScoreMatrix BwdB(n+1,m+1);
    InftyScoreRRMatrix Bwd(ranges); 
    InftyScoreRRMatrix BwdA(ranges); 
    InftyScoreRRMatrix BwdB(ranges);

    backward_algorithm(home,Bwd,BwdA,BwdB,UBM);
    //if (debug_out) {std::cout << "Bwd"<<std::endl<<Bwd <<std::endl;}

    if(debug_out && Fwd(n,m).is_finite() && Bwd(0,0).is_finite() && Fwd(n,m).finite_value()!=Bwd(0,0).finite_value()) {
	std::cerr << "Clash"<<std::endl;
	exit(-1);
    }


    if (false && debug_out) {
	// ------------------------------------------------------------
	// DEBUGGING: check trace
	for (size_t i=1; i<=n; i++) {
	    std::cout << "TRACE " << i <<" " << traceA[i] <<" ";
	    if (traceA[i]!=0) {
		size_t j = traceA[i];
		score_t matchscore=evaluate_tracematch(traceA,traceB,considered_ams,match_scores,i,j);
		
		cout
		    << match_arrow_allowed(i,j)<<" "
		    <<Fwd(i-1,j-1)<<"+"<< matchscore<<"["<<ub_match(i,j,considered_ams,match_scores)<<"]"<<"+"<<Bwd(i,j)<<" = "
		    <<Fwd(i-1,j-1)+matchscore+Bwd(i,j);
	    }
	    std::cout << std::endl;
	}
    }
    
    if (debug_out && Gecode::me_failed(ret)) {
	std::cerr << "Fail before pruning."<<std::endl;
	return Gecode::ES_FAILED;
    }
    
    // -------------------- PRUNE VARIABLES
    ret |= prune(home,
		 Fwd,
		 FwdA,
		 //FwdB,
		 Bwd,
		 BwdA,
		 //BwdB,
		 UBM);
    
    if (debug_out) {
	std::cout << "After Pruning:"<<std::endl;
	print_vars(std::cout);
    }
    if (Gecode::me_failed(ret)) {
	if (debug_out) std::cout << "Fail when pruning."<<std::endl;
	return Gecode::ES_FAILED;
    }

    
    // test whether all vars fixed, then subsume (can we subsume earlier?)
    if ( all_vars_fixed() ) {
	return home.ES_SUBSUMED(*this);
    }

    // -------------------- select CHOICE for the space
    if (!Gecode::me_modified(ret)) {
	// don't call choice if propagate will be called again anyway (due to modifications)
	choice(static_cast<RNAalignment&>(home),Fwd,Bwd,traceA,traceB,trace_score,UBM,match_scores);
    }

    // -------------------- pass some debugging information to space
    

    // -------------------- return
    return Gecode::me_modified(ret)?Gecode::ES_NOFIX:Gecode::ES_FIX;
}
