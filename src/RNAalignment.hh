#ifndef RNA_ALIGNMENT
#define RNA_ALIGNMENT

#include <fstream>
#include <iostream>
#include "WinHandler.hh"
#include "alignment_score.hh"
#include "LocARNA/alignment.hh"


// EXPERIMENTAL: therefore only global var
//const size_t discrepancy_limit=7;


const bool custom_branching=true;

typedef std::vector<int>::size_type size_type;

class RNAalignBrancher;


/**
 * \brief %RNA Alignment
 *
 *
 *
 */
class RNAalignment : public Gecode::Space {
    friend class AlignmentScore;
      
protected:
    
    // only used for output
    const LocARNA::Sequence &seqA; 
    const LocARNA::Sequence &seqB;
    const LocARNA::ArcMatches &arcmatches;
    
    const size_t n;
    const size_t m;
    
public:
    
    //! choice class for the brancher RNAalignBranch. Confer commit
    //! method of the brancher for how this data is actually used.
    class ChoiceData {
    public:
	
	bool enum_M;   //!< whether to enumerate an M or MD variable
	size_t pos;    //!< position of variable in M or MD
	size_t val;    //!< selected value of variable (unused by current strategy)
	size_t minval; //!< minimal selected domain value
	size_t maxval; //!< maximal selected domain value
	std::vector<int> values; //!< selected domain values
	
	// An int vector is used for constructing the choice of domain values
	// by the AlignmentScore propagator and communicating this to the
	// space.  This is definitely not optimized for speed!!! One should
	// better use ranges instead of single values and use specific Gecode
	// mechanisms for doing this. Using the vector is a workaround, since
	// I (SW) could not make Gecode's BitSet working (or what else is more
	// appropriate?).
	
	
	ChoiceData()
	    : enum_M(false), pos(0),val(0), minval(0), maxval(0), values()
	{}
	
	ChoiceData(const ChoiceData &cd)
	    : enum_M(cd.enum_M), pos(cd.pos),val(cd.val), 
	      minval(cd.minval), maxval(cd.maxval), values(cd.values)
	{}
	
    };
    
protected:
   
    WinHandler* wind;
    
    //! MD[i] is position of match or deletion in row i. We model only
    //! match and deletion positions explicitely. Insertion positions are
    //! derived! Compare the handling of MD and M by the propagator AlignmentScore. 
    Gecode::IntVarArray MD;  
    Gecode::BoolVarArray M; //!< M[i] is true iff i~MD[i] is a match
    
    Gecode::IntVar Score;

    //size_t discrepancy; // experimental for simulating LDS
    
    //! the propagator AlignmentScore selects the choice for the brancher RNAalignBranch 
    //! and communcates its choice in this field of the space
    ChoiceData choice_data;

 
public:
    RNAalignment(const LocARNA::Sequence &seqA_, const LocARNA::Sequence &seqB_,
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

    /// Constructor for cloning \a s
    RNAalignment(bool share, RNAalignment& s) : 
	Gecode::Space(share,s),
	seqA(s.seqA),
	seqB(s.seqB),
	arcmatches(s.arcmatches),
	n(s.n),
	m(s.m),
	wind(s.wind),
	choice_data(s.choice_data)
	//, discrepancy(s.discrepancy)
    {
	MD.update(*this, share, s.MD);
	M.update(*this, share, s.M);
	Score.update(*this,share,s.Score);
    }
    
    /// Generate a new copy of the space
    virtual Gecode::Space*
    copy(bool share) {
	return new RNAalignment(share,*this);
    }
    
    //! test whether all variables M and MD are assigned
    bool all_assigned() const;
    
    //! convert variable valuation to Alignment object
    //! pre: all variables MD and M are assigned
    LocARNA::Alignment
    to_alignment() const;
    
    //! print solution in clustal format 
    void
    print_clustal_format(std::ostream& out_s) const;

    //! print solution in pp format
    void 
    print_pp_format(std::ostream& out_s,
		    const LocARNA::BasePairs& bpsA, const LocARNA::BasePairs& bpsB, 
		    const LocARNA::Scoring& scoring, 
		    const LocARNA::AnchorConstraints& seq_constraints) const;

    /// Print solution
    virtual 
    void
    print(std::ostream& out) const;
    
    virtual void constrain(const Space& _best) {
	const RNAalignment& best = static_cast<const RNAalignment&>(_best);
	rel(*this,Score,Gecode::IRT_GR,best.Score);
    }

    /**
     * \brief %Custom Brancher for RNA Alignment
     *
     *
     *
     */
    class RNAalignBrancher : public Gecode::Brancher {
    private:
	mutable size_t start;

	class Choice : public Gecode::Choice {
	    // Must be refined by inheritance such that the information
	    // stored inside a choice is sufficient to redo a commit
	    // performed by a particular brancher.
	    
	public:
	    const RNAalignment::ChoiceData cd;
	    
	    Choice(const Brancher& b, const RNAalignment::ChoiceData &cd_)
		: Gecode::Choice(b,2),
		  cd(cd_)
	    {}
	    
	    virtual size_t size(void) const {
		return sizeof(Choice);
	    }
	};
	
	
	RNAalignBrancher(Gecode::Space &home) : Gecode::Brancher(home),start(1) {}
	RNAalignBrancher(Gecode::Space &home,bool share,RNAalignBrancher &b) 
	    : Gecode::Brancher(home,share,b),start(b.start) {}
    
    public:

	// return true, if there are unassigned vars left for this brancher
	//        false, otherwise
	virtual bool status(const Gecode::Space& home) const {
	    const RNAalignment& s = static_cast<const RNAalignment&>(home);
	    
	    for (size_t i=start; i<(size_type)s.MD.size(); i++) {
		if (! s.MD[i].assigned()) {start=i;return true;}
	    }
	    return false;
	}
	
	// returns choice
	virtual Gecode::Choice* choice(Gecode::Space& home) {
	    const RNAalignment& s = static_cast<const RNAalignment&>(home);
	
	    // the decision about choice is moved to the propagator, which
	    // returned it's choice in the field choice of the space s
	    
	    //std::cout << "CHOICE "<<s.pos<<" "<<s.val<<" "<<s.minval<<" "<<s.maxval<<" "<<s.MD[s.pos]<<std::endl;
	    return new Choice(*this,s.choice_data);
	}
	
	// commits the choice c and alternative a 
	virtual Gecode::ExecStatus commit(Gecode::Space& home, 
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
    
	virtual Gecode::Actor* copy(Gecode::Space& home, bool share) {
	    return new (home) RNAalignBrancher(home, share, *this);
	}
    
	static void post(RNAalignment& home) {
	    (void) new (home) RNAalignBrancher(home);
	}
    
	virtual size_t dispose(Gecode::Space&) {
	    return sizeof(*this);
	}
    };
};
#endif // RNA_ALIGNMENT
