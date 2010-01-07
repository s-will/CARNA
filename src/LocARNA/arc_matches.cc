#include "arc_matches.hh"

#include <fstream>
#include <sstream>
#include <algorithm>

class tuple5 {
public:
  typedef std::vector<int>::size_type size_type;

  size_type i;
  size_type j;
  size_type k;
  size_type l;
  score_t score;
  tuple5(size_type i_,size_type j_,size_type k_,size_type l_,score_t score_)
    :i(i_), j(j_), k(k_), l(l_), score(score_)
  {}
};


bool ArcMatches::is_valid_arcmatch(const Arc &arcA,const Arc &arcB) const {
  return
    (size_type)abs((int)arcA.left()-(int)arcB.left()) <= max_pos_diff
    &&
    (size_type)abs((int)arcA.right()-(int)arcB.right()) <= max_pos_diff
    &&
    (size_type)abs((int)(arcA.right()-arcA.left()) - (int)(arcB.right()-arcB.left())) <= max_length_diff	
    &&
    constraints.allowed_edge(arcA.left(),arcB.left())
    &&
    constraints.allowed_edge(arcA.right(),arcB.right());	
}


void 
ArcMatches::init_inner_arc_matchs() {
  inner_arcmatch_idxs.resize(arc_matches_vec.size());
    
  for (ArcMatch::idx_type i=0; i<arc_matches_vec.size(); ++i) {
	
    const ArcMatch &am=arc_matches_vec[i];
    const Arc &arcA=am.arcA();
    const Arc &arcB=am.arcB();
	
    // set invalid, in case we don't find an inner arc
    inner_arcmatch_idxs[i] = arc_matches_vec.size();
	
    // find index of inner arc match
    const ArcMatchIdxVec &list =common_left_end_lists(arcA.left()+1,arcB.left()+1);
    for ( ArcMatchIdxVec::const_iterator it=list.begin(); list.end()!=it; ++it) {
	if ((arcmatch(*it).arcA().right()==arcA.right()-1)
	    &&
	    ((arcmatch(*it).arcB().right()==arcB.right()-1)))
	    {
		inner_arcmatch_idxs[i] = *it;
		break;
	    }
    }
  }
}


void
ArcMatches::sort_right_adjacency_lists() {
    
  //sorted_common_right_end_lists.resize(lenA+1,lenB+1);
    
  for (size_type i=1; i<=lenA ; i++) {
    for (size_type j=1; j<=lenB ; j++) {
	    
      ArcMatchIdxVec &list = common_right_end_lists(i,j);
	    
      std::sort(list.begin(),
		list.end(),
		lex_greater_left_ends(*this));
	    
      // 	    // ------------------------------------------------------------
      // 	    // generate the "sorted list" special data structure
      // 	    //
	    
      // 	    std::vector<ArcMatchVec> &sorted_list = sorted_common_right_end_lists(i,j);
	    
      // 	    ArcMatchVec::const_iterator start_it=list.begin();
      // 	    ArcMatchVec::const_iterator it=list.begin();

      // 	    for(; list.end()!=it;) {
      // 		if (start_it->arcA().left()!=it->arcA().left()) {
      // 		    sorted_list.push_back(ArcMatchVec(start_it,it));
      // 		    start_it=it;
      // 		}
      // 		++it;
      // 	    }
      // 	    if (list.size()>0)
      // 		sorted_list.push_back(ArcMatchVec(start_it,it));
    }
  }
}


ArcMatches::ArcMatches(const Sequence &seqA_,
		       const Sequence &seqB_,
		       const std::string &arcmatch_scores_file,
		       int probability_scale,
		       size_type max_length_diff_, 
		       size_type max_pos_diff_,
		       const AnchorConstraints &constraints_
		       )
  : lenA(seqA_.length()),
    lenB(seqB_.length()),
    max_length_diff(max_length_diff_),
    max_pos_diff(max_pos_diff_),
    constraints(constraints_),
    maintain_explicit_scores(true)
{
  read_arcmatch_scores( arcmatch_scores_file, probability_scale );
}

ArcMatches::ArcMatches(const RnaData &rnadataA,
		       const RnaData &rnadataB,
		       double min_prob,
		       size_type max_length_diff_, 
		       size_type max_pos_diff_,
		       const AnchorConstraints &constraints_
		       )
  : lenA(rnadataA.get_sequence().length()),
    lenB(rnadataB.get_sequence().length()),
    bpsA(new BasePairs(&rnadataA,min_prob)),
    bpsB(new BasePairs(&rnadataB,min_prob)),
    max_length_diff(max_length_diff_),
    max_pos_diff(max_pos_diff_),
    constraints(constraints_),
    maintain_explicit_scores(false)
{

  // ----------------------------------------
  // initialize the vector for arc matchs and adjacency lists
    
  // iterate through all pairs of arc matches and 
  // generate entries for all relevant arc matches.
  //
  // Note that the probabilities of the arcs are already ok 
  // due to the filtering by BasePairs class.
  // Here, we will only check for difference heuristics

  common_left_end_lists.resize(lenA+1,lenB+1);
  common_right_end_lists.resize(lenA+1,lenB+1);
    
  for(size_type i=0; i<bpsA->num_bps(); i++) {
    const Arc *arcA = &bpsA->arc(i);
		    
    for(size_type j=0; j<bpsB->num_bps(); j++) {
			
      const Arc *arcB = &bpsB->arc(j);
	    
      // check whether arc match is valid
      if (!is_valid_arcmatch(*arcA,*arcB)) continue;    
	    
      size_type idx = arc_matches_vec.size();
	    
      // make entry in arc matches
      arc_matches_vec.push_back(ArcMatch(*arcA,*arcB,idx));
      
      // make entries in adjacency lists
      common_left_end_lists(arcA->left(),arcB->left()).push_back(idx);
      common_right_end_lists(arcA->right(),arcB->right()).push_back(idx);
    }
  }
    
  init_inner_arc_matchs();
    
  sort_right_adjacency_lists();
}

void ArcMatches::read_arcmatch_scores( const std::string &arcmatch_scores_file, int probability_scale ) {
    
  // try to open file
  std::ifstream in(arcmatch_scores_file.c_str());
    
  // catch error while opening
  if (!in.good()) {
    std::cerr << "Cannot open file "<<arcmatch_scores_file<<" for reading arcmatch-scores."<<std::endl;
    exit(-1);
  }
     
  size_type i;
  size_type j;
  size_type k;
  size_type l;
  score_t score;
    
  BasePairs::bpair_set_t arcsA;
  BasePairs::bpair_set_t arcsB;

  std::vector<tuple5> lines;
    
  // read all lines
  std::string line;
  size_type lineno=0;
  while( getline(in,line) ) {
    std::istringstream in(line);
    lineno++;
	
    in >> i >> j >> k >> l;
	
    if (probability_scale<0) {
      in >> score;
    } else {
      double prob;
      in >> prob;
      score = (score_t) (prob * (double)probability_scale);
    }
	
    if ( i<=0 || j<=0 || k<=0 || l<=0 
	 || i>j || j>lenA
	 || k>l || l>lenB)
      {
	std::cerr <<"Cannot read arc match scores. Invalid line "<<lineno <<": " <<line <<std::endl; 
	exit(-1);
      }
	
    lines.push_back(tuple5(i,j,k,l,score));
	
    arcsA.insert(BasePairs::bpair_t(i,j));
    arcsB.insert(BasePairs::bpair_t(k,l));
  }

  bpsA = new BasePairs(lenA, arcsA);
  bpsB = new BasePairs(lenB, arcsB);
    
  // ----------------------------------------
  // construct the vectors of arc matches and scores
    
  common_left_end_lists.resize(lenA+1,lenB+1);
  common_right_end_lists.resize(lenA+1,lenB+1);
    
  for (std::vector<tuple5>::iterator it=lines.begin(); lines.end()!=it; ++it) {
    const Arc &arcA=bpsA->arc(it->i,it->j);
    const Arc &arcB=bpsB->arc(it->k,it->l);

    // check whether arc match is valid
    if (!is_valid_arcmatch(arcA,arcB)) continue;    

    size_type idx = arc_matches_vec.size();

    arc_matches_vec.push_back(ArcMatch(arcA,arcB,idx));
    scores.push_back(it->score); // now the score has the same index as the corresponding arc match
	
    // make entries in adjacency lists
    common_left_end_lists(arcA.left(),arcB.left()).push_back(idx);
    common_right_end_lists(arcA.right(),arcB.right()).push_back(idx);
  }

  init_inner_arc_matchs();
    
  sort_right_adjacency_lists();
}


void ArcMatches::write_arcmatch_scores(const std::string &arcmatch_scores_file, const Scoring &scoring) const {
  // try to open file
  std::ofstream out(arcmatch_scores_file.c_str());
    
  // catch error while opening
  if (!out.good()) {
    std::cerr << "Cannot open file "<<arcmatch_scores_file<<" for writing arcmatch-scores."<<std::endl;
    exit(-1);
  }
    
  for (size_type i=0; i<num_arc_matches(); i++) {
	
    const Arc &arcA = arcmatch(i).arcA();
    const Arc &arcB = arcmatch(i).arcB();
    size_type al = arcA.left();
    size_type ar = arcA.right();
    size_type bl = arcB.left();
    size_type br = arcB.right();
	
    score_t score = scoring.arcmatch(arcmatch(i));
	
    out << al << " " << ar << " " << bl << " " << br << " " << score <<"\n";
  }
}


void ArcMatches::get_max_right_ends(size_type al,size_type bl,size_type *max_ar,size_type *max_br, bool no_lonely_pairs) const {     
  
  // for no lonely pairs, al,bl is the beginning of the inner base pair match
  // we need to find the largest base pair match with beginning al-1,bl-1
  // that has an inner bpm
  if (no_lonely_pairs) {
    al--;
    bl--;
  }

  for(ArcMatchIdxVec::const_iterator it=common_left_end_list(al,bl).begin();
      common_left_end_list(al,bl).end() != it; ++it ) {
	
    const ArcMatch &am = arcmatch(*it);

    // if lonely pairs are forbidden, consider only arc matchs that at least have an inner match
    if (no_lonely_pairs && ! exists_inner_arc_match(am)) continue;
	
	
    *max_ar=std::max(*max_ar,am.arcA().right());
    *max_br=std::max(*max_br,am.arcB().right());
  }

  if (no_lonely_pairs) {
    *max_ar--;
    *max_br--;
  }
}

void ArcMatches::get_min_right_ends(size_type al,size_type bl,size_type *min_ar,size_type *min_br) const { 
    
  for(ArcMatchIdxVec::const_iterator it=common_left_end_list(al,bl).begin();
      common_left_end_list(al,bl).end() != it; ++it ) {
    
    const ArcMatch &am = arcmatch(*it);
    
    *min_ar=std::min(*min_ar,am.arcA().right());
    *min_br=std::min(*min_br,am.arcB().right());
  }
}
