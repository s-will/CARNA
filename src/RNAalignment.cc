#include "RNAalignment.hh"

#include <LocARNA/rna_ensemble.hh>

using namespace LocARNA;


template<class T>
Gecode::Archive &
operator << (Gecode::Archive &e, const std::vector<T> &v) {
    e << (unsigned int)v.size();
    for (size_t i=0; i<v.size(); i++) e << v[i];
    return e;
}

template<class T>
Gecode::Archive &
operator >> (Gecode::Archive &e, std::vector<T> &v) {
    unsigned int size;
    e >> size;
    v.resize(size);
    
    for (size_t i=0; i<size; i++) e >> v[i];
    return e;
}

Gecode::Archive &
operator << (Gecode::Archive &e, const RNAalignment::ChoiceData &cd) {
    return
	e << cd.enum_M
	  << cd.pos
	  << cd.val
	  << cd.new_lower_bound
	  << cd.values
	  << cd.best_traceA
	  << cd.best_traceB
	  << (int)cd.best_trace_score;
}

Gecode::Archive &
operator >> (Gecode::Archive &e, RNAalignment::ChoiceData &cd) {
    int bts;

    e >> cd.enum_M
      >> cd.pos
      >> cd.val
      >> cd.new_lower_bound
      >> cd.values
      >> cd.best_traceA
      >> cd.best_traceB
      >> bts;
    cd.best_trace_score=bts;
    
    return e;
}
								      


RNAalignment::RNAalignment(const LocARNA::Sequence &seqA_, const LocARNA::Sequence &seqB_,
			   const LocARNA::ArcMatches &arcmatches_,
			   const LocARNA::AlignerParams &aligner_params, 
			   const LocARNA::Scoring &scoring,
			   int lower_score_bound,
			   int upper_score_bound,		 
			   bool opt_graphical_output)
    :
    seqA(seqA_),
    seqB(seqB_),
    arcmatches(arcmatches_),
    n(seqA.length()),
    m(seqB.length()),
    MD(*this,n+1,0,m), //we only need MD_1,...,MD_n ==> ignore MD_0
    M(*this,n+1,0,1),
    Score(*this,lower_score_bound,upper_score_bound),
    choice_data()
    //, discrepancy(0)
{
#  ifdef HAVE_GIST
    if (opt_graphical_output) 
	wind=new WinHandler(n+1,m+1,"Display variables status");
    else {
	wind=NULL;
    }
#  endif // HAVE_GIST


    //ignore MD_0 and M_0
    rel(*this,MD[0],Gecode::IRT_EQ,0);
    rel(*this,M[0],Gecode::IRT_EQ,1);
	
    // impose anchor constaints
    for(size_t i=1;i<=n; i++){
	int j=aligner_params.constraints.match_to_a(i);
	if (j>0) {
	    // there is an anchor i~j
	    rel(*this,MD[i],Gecode::IRT_EQ,j);
	    rel(*this,M[i],Gecode::IRT_EQ,1);
	}
    }
    
    /* restrict domains according to TraceController
       
       The trace controller gives information about the possible
       trace cells in each row of the alignment matrix.
    
       We show the relation of trace controller and domains by example:
       Assume in rows i-1 and i, we have the trace cells as marked by +
    
          |0123456789
       i-1|   ++++
       i  |     ++++
       
       I.e., tc.min_col(i-1)=3, tc.max_col(i-1)=6, tc.min_col(i)=5, tc.max_col(i)=8
       We conclude that
       dom(MD[i]) = [5..7] ! ( there is no deletion or match of i after or with j=8
       and dom(M[i])=[0..1]

       From
          |0123456789
       i-1|  +++
       i  |     ++++
       
       we conclude that
       dom(MD[i])={5} and dom(M[i])={1}
    */
    const LocARNA::TraceController &tc = aligner_params.trace_controller;

    for (size_t i=1; i<=n; i++) {
	// MD[i] :>= tc.min_col(i)
	rel(*this,MD[i],Gecode::IRT_GQ,tc.min_col(i));
	// MD[i] :<= tc.max_col(i-1)+1
	rel(*this,MD[i],Gecode::IRT_LQ,tc.max_col(i-1));
	
	if (tc.max_col(i-1)==tc.min_col(i)+1) {
	    rel(*this,M[i],Gecode::IRT_EQ,1); // M[i]:=1
	}
	if (tc.min_col(i-1)==tc.min_col(i) && tc.min_col(i)==tc.max_col(i)) {
	    rel(*this,M[i],Gecode::IRT_EQ,0); // M[i]:=0
	}
    }
    
    

    AlignmentScore::post(*this,seqA,seqB,arcmatches,aligner_params,scoring,
			 MD,M,Score);
	
    // simple relations
	
    for (size_t i=1; i<=n; i++) {
	// the MD variables are sorted
	rel(*this,MD[i],Gecode::IRT_GQ,MD[i-1]);
	    
	// M[i] implies M[i]>M[i-1]
	Gecode::BoolVar greater(*this,0,1);
	rel(*this,MD[i],Gecode::IRT_GR,MD[i-1],greater);
	rel(*this,M[i],Gecode::IRT_LQ,greater); // M[i] implies greater
    }
	
	
    // test the case where a good score is known
    //rel(*this,Score,Gecode::IRT_GR, 33592);
	
    if (custom_branching) {
	RNAalignBrancher::post(*this);
    } else {
	//  first enumerate MD, then the rest
	    
	// resort MD vector, such that we start enumerating in the middle
	// Ideally sort like balanced binary tree in array
	Gecode::IntVarArgs MD_resorted(MD.size());
	for (size_t i=0; i<(size_t)MD.size(); i++) { 
	    MD_resorted[i] = MD[(i+MD.size()/2)%MD.size()];
	}
	Gecode::branch(*this, MD_resorted, Gecode::INT_VAR_SIZE_MAX(), Gecode::INT_VAL_MED());
    }
	
    Gecode::branch(*this, M, Gecode::INT_VAR_SIZE_MAX(), Gecode::INT_VAL_MAX());
	
}



bool RNAalignment::all_assigned() const {
    bool all_assigned=true;

    for (size_t i=1; i<=seqA.length(); i++) {
	all_assigned &= MD[i].assigned();
	all_assigned &= M[i].assigned();
    }

    return all_assigned;
}

Alignment
RNAalignment::to_alignment() const {

    Alignment alignment(seqA,seqB);

    size_t j=1;
    for (size_t i=1; i<=seqA.length(); i++) {

	for (;j<(size_t)MD[i].val();j++) {
	    alignment.append(-1,j);
	}

	if (M[i].val()==1) {
	    alignment.append(i,j);
	    j++;
	} else {
	    alignment.append(i,-1);
	}
    }

    for (;j<=seqB.length();j++) {
	alignment.append(-1,j);
    }

    return alignment;
}


void
RNAalignment::print(std::ostream& out) const {

    size_t width=60;

    if (all_assigned()) {
	LocARNA::Alignment alignment=to_alignment();
	MultipleAlignment ma(alignment, false);
	ma.write(out,width);

    } else {
	out << "Matches/Deletions:    ";
	for (size_t i=0; i<=seqA.length(); i++) {
	    if (!(M[i].assigned() && MD[i].assigned())) {
		out <<i<<(M[i].assigned()?(M[i].val()==0?"g":"~"):"?")<<MD[i]<<", ";
	    }
	}
	out << std::endl;
    }

    out << "Score:      " << Score << std::endl;


    // ////////////////////////////////////////
    // list undecided base pairs

    for (ArcMatches::const_iterator it=arcmatches.begin(); arcmatches.end()!=it; ++it) {
	const Arc &arcA=it->arcA();
	const Arc &arcB=it->arcB();

	if (
	    // arc match possible
	    M[arcA.left()].in(1) && MD[arcA.left()].in(arcB.left())
	    &&
	    M[arcA.right()].in(1) && MD[arcA.right()].in(arcB.right())

	    && // not both ends are fixed
	    !((M[arcA.left()].assigned() && MD[arcA.left()].assigned())
	      || (M[arcA.right()].assigned() && MD[arcA.right()].assigned()))
	     ) {

	    out << arcA << "?" << arcB <<" ";
	}
    }
    out << std::endl;

#  ifdef HAVE_GIST
    // ////////////////////////////////////////
    // update display
    if (wind!=NULL) wind->update(MD,M);
#  endif // HAVE_GIST

}

Gecode::ModEvent
RNAalignment::RNAalignBrancher::
fix_vars_from_trace(Gecode::Space& home,
		    const std::vector<unsigned int> &traceA,
		    const std::vector<unsigned int> &traceB
		    ) const {
    RNAalignment& s = static_cast<RNAalignment&>(home);

    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
    
    int last=0;
    
    for (size_t i=1; i<traceA.size(); i++) {
	if (traceA[i]!=0) { // match i
	    ret |= Gecode::Int::IntView(s.MD[i]).eq(home,(int)traceA[i]);
	    ret |= Gecode::Int::BoolView(s.M[i]).eq(home,1);
	    last=traceA[i];
	} else { // deletion of i
	    ret |= Gecode::Int::IntView(s.MD[i]).eq(home,last);
	    ret |= Gecode::Int::BoolView(s.M[i]).eq(home,0);
	}
    }
    return ret;
}
					       

Gecode::ExecStatus 
RNAalignment::RNAalignBrancher::commit(Gecode::Space& home, 
			     const Gecode::Choice& _c,
			     unsigned int a) {
    RNAalignment& s = static_cast<RNAalignment&>(home);
    const Choice& c = static_cast<const Choice&>(_c);
    
    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
	    
    const RNAalignment::ChoiceData &cd=c.cd;
	    
    if (cd.new_lower_bound) {
	s.choice_data.new_lower_bound=false;
	
	if (a==0) {
	    //std::cout << "Commit to new better trace: "
	    //<< cd.best_trace_score << " " << s.Score << std::endl;
	
	    ret = Gecode::Int::IntView(s.Score).eq(home, (int) cd.best_trace_score);
	    
	    ret |= fix_vars_from_trace(home,cd.best_traceA,cd.best_traceB);
	    
	} else {
	    ret = Gecode::Int::IntView(s.Score).gr(home, (int) cd.best_trace_score);
	}
    }
    else if (cd.enum_M) {
	if (a==0) {
	    ret = Gecode::Int::BoolView(s.M[cd.pos]).eq(home, 1);
	} else {
	    ret = Gecode::Int::BoolView(s.M[cd.pos]).eq(home, 0);
	}
    } else {
		
	//Gecode::Iter::Ranges::Singleton r((int)cd.minval,(int)cd.maxval);
		
	// note the const cast is necessary due to an ill specified Gecode interface.  
	// &cd.values[0] points to the array encapsulated by the vector,
	// this is a HACK, since it depends on the stl implementation.
	// However, I would rather blame Gecode to require an C-style array here:)
	Gecode::Iter::Values::Array values_iter(const_cast<int *>(&cd.values[0]),cd.values.size());
	
	// std::cerr << "minval "<<cd.minval<<std::endl;
	// std::cerr << "maxval "<<cd.maxval<<std::endl;
	// for (size_t i=0; i<cd.values.size(); i++)
	//     std::cerr << cd.values[i]<<" ";
	// std::cerr << std::endl;
	// std::cerr << "before MD["<<cd.pos<<"] " << s.MD[cd.pos]<<std::endl;
	if (a==0) {
	    //ret = Gecode::Int::IntView(s.MD[cd.pos]).inter_r(home, r,false);
	    ret = Gecode::Int::IntView(s.MD[cd.pos]).inter_v(home, values_iter, false);
	} else {
	    //ret = Gecode::Int::IntView(s.MD[cd.pos]).minus_r(home, r, false);
	    ret = Gecode::Int::IntView(s.MD[cd.pos]).minus_v(home, values_iter, false);
		    
	    // EXPERIMENTAL limiting of discrepancy
	    // s.discrepancy++;
	    // if (s.discrepancy>discrepancy_limit) { return Gecode::ES_FAILED; } 
	}
	//std::cerr << "after MD["<<cd.pos<<"] " << s.MD[cd.pos]<<std::endl;
    }
	
    return Gecode::me_failed(ret)
	? Gecode::ES_FAILED
	: Gecode::ES_OK;
}




