#include "RNAalignment.hh"

using namespace LocARNA;

RNAalignment::RNAalignment(const LocARNA::Sequence &seqA_, const LocARNA::Sequence &seqB_,
	     const LocARNA::ArcMatches &arcmatches_,
	     const LocARNA::AlignerParams &aligner_params, 
	     const LocARNA::Scoring &scoring,
	     bool opt_graphical_output)
    :
    seqA(seqA_),
    seqB(seqB_),
    arcmatches(arcmatches_),
    n(seqA.length()),
    m(seqB.length()),
    MD(*this,n+1,0,m), //we only need MD_1,...,MD_n ==> ignore MD_0
    M(*this,n+1,0,1),
    Score(*this,Gecode::Int::Limits::min,Gecode::Int::Limits::max),
    choice_data()
    //, discrepancy(0)
{
    if (opt_graphical_output) 
	wind=new WinHandler(n+1,m+1,"Display variables status");
    else {
	wind=NULL;
    }


    //ignore MD_0 and M_0
    rel(*this,MD[0],Gecode::IRT_EQ,0);
    rel(*this,M[0],Gecode::IRT_EQ,1);
	
    /*
    //update domains of MD according to anchor constaints
    for(size_t i=1;i<=n; i++){
    AnchorConstraints::size_pair_t range = 
    aligner_params.constraints.get_range_seqA()[i];
    rel(*this,MD[i],Gecode::IRT_GQ,range.first);
    rel(*this,MD[i],Gecode::IRT_LQ,range.second);
    }
    */

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
	for (size_type i=0; i<(size_type)MD.size(); i++) { 
	    MD_resorted[i] = MD[(i+MD.size()/2)%MD.size()];
	}
	Gecode::branch(*this, MD_resorted, Gecode::INT_VAR_SIZE_MAX, Gecode::INT_VAL_MED);
    }
	
    Gecode::branch(*this, M, Gecode::INT_VAR_SIZE_MAX, Gecode::INT_VAL_MAX);
	
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

void RNAalignment::print_clustal_format(std::ostream& out_s) const{

    size_t width=60;

    if (all_assigned()) {
	to_alignment().write_clustal(out_s,
				     width,
				     (LocARNA::infty_score_t)Score.val(),
				     false,false,true,false);
    }
}

void
RNAalignment::print_pp_format(std::ostream& out_s,
			      const LocARNA::BasePairs& bpsA, const LocARNA::BasePairs& bpsB,
			      const LocARNA::Scoring& scoring,
			      const LocARNA::AnchorConstraints& seq_constraints) const {

    size_t width=60;

    if (all_assigned()) {

	to_alignment().write_pp(out_s,
				bpsA, bpsB,
				scoring,
				seq_constraints,
				width);
    }

}


void
RNAalignment::print(std::ostream& out) const {

    size_t width=60;

    if (all_assigned()) {

	to_alignment().write(out,
			     width,
			     (LocARNA::infty_score_t)Score.val()
			     );

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


    // ////////////////////////////////////////
    // update display
    if (wind!=NULL) wind->update(MD,M);

}

Gecode::ExecStatus 
RNAalignment::RNAalignBrancher::commit(Gecode::Space& home, 
			     const Gecode::Choice& _c,
			     unsigned int a) {
    RNAalignment& s = static_cast<RNAalignment&>(home);
    const Choice& c = static_cast<const Choice&>(_c);
    
    /*
    // split in eq and nq
    if (a==0) {
    return Gecode::me_failed(Gecode::Int::IntView(s.MD[c.pos]).eq(home, (int)c.val))
    ? Gecode::ES_FAILED
    : Gecode::ES_OK;
    } else {
    return Gecode::me_failed(Gecode::Int::IntView(s.MD[c.pos]).nq(home, (int)c.val))
    ? Gecode::ES_FAILED
    : Gecode::ES_OK;
    }
    */
	    
    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
	    
    const RNAalignment::ChoiceData &cd=c.cd;
	    
    if (cd.enum_M) {
	if (a==0) {
	    ret = Gecode::Int::BoolView(s.M[cd.pos]).eq(home, 1);
	} else {
	    ret = Gecode::Int::BoolView(s.M[cd.pos]).eq(home, 0);
	}
    } else {
		
	Gecode::Iter::Ranges::Singleton r((int)cd.minval,(int)cd.maxval);
		
	// note the const cast is necessary due to an ill specified Gecode interface.  
	// &cd.values[0] points to the array encapsulated by the vector,
	// this is a HACK, since it depends on the stl implementation.
	// However, I would rather blame Gecode to require an C-style array here:)
	Gecode::Iter::Values::Array values_iter(const_cast<int *>(&cd.values[0]),cd.values.size());
		
	if (a==0) {
	    ret = Gecode::Int::IntView(s.MD[cd.pos]).inter_r(home, r,false);
	    ret = Gecode::Int::IntView(s.MD[cd.pos]).inter_v(home, values_iter, false);
	} else {
	    ret = Gecode::Int::IntView(s.MD[cd.pos]).minus_r(home, r, false);
	    ret = Gecode::Int::IntView(s.MD[cd.pos]).minus_v(home, values_iter, false);
		    
	    // EXPERIMENTAL limiting of discrepancy
	    // s.discrepancy++;
	    // if (s.discrepancy>discrepancy_limit) { return Gecode::ES_FAILED; } 
	}
    }
	
    return Gecode::me_failed(ret)
	? Gecode::ES_FAILED
	: Gecode::ES_OK;
}




