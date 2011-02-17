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
    
    bool enum_M;
    size_t pos;
    size_t val;
    size_t minval;
    size_t maxval;

    WinHandler* wind;
    
    //! MD[i] is position of match or deletion in row i. We model only
    //! match and deletion positions explicitely. Insertion positions are
    //! derived! Compare the handling of MD and M by the propagator AlignmentScore. 
    Gecode::IntVarArray MD;  
    Gecode::BoolVarArray M; //!< M[i] is true iff i~MD[i] is a match
    
    Gecode::IntVar Score;

    //size_t discrepancy; // experimental for simulating LDS


 
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
	Score(*this,Gecode::Int::Limits::min,Gecode::Int::Limits::max)
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
	enum_M(s.enum_M),
	pos(s.pos),
	val(s.val),
	minval(s.minval),
	maxval(s.maxval)
	//, discrepancy(s.discrepancy)
    {
	MD.update(*this, share, s.MD);
	M.update(*this, share, s.M);
	Score.update(*this,share,s.Score);
	wind=s.wind;
    }
    
    /// Perform copying during cloning
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
	private:
	public:
	    bool enum_M;
	    size_t pos;
	    size_t val;
	    size_t minval;
	    size_t maxval;
	

	    Choice(const Brancher& b,bool enum_M_,size_t pos_, size_t val_,size_t minval_,size_t maxval_)
		: Gecode::Choice(b,2),
		  enum_M(enum_M_),
		  pos(pos_),val(val_),
		  minval(minval_),
		  maxval(maxval_)
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
	
	    // move decision about choice to the propagator
	
	    //std::cout << "CHOICE "<<s.pos<<" "<<s.val<<" "<<s.minval<<" "<<s.maxval<<" "<<s.MD[s.pos]<<std::endl;
	    return new Choice(*this,s.enum_M,s.pos,s.val,s.minval,s.maxval);
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

	    if (c.enum_M) {
		if (a==0) {
		    ret = Gecode::Int::BoolView(s.M[c.pos]).eq(home, 1);
		} else {
		    ret = Gecode::Int::BoolView(s.M[c.pos]).eq(home, 0);
		}
	    } else {
	    
	    
		Gecode::Iter::Ranges::Singleton r((int)c.minval,(int)c.maxval);

		if (a==0) {
		    ret = Gecode::Int::IntView(s.MD[c.pos]).inter_r(home, r,false);
		} else {
		    ret = Gecode::Int::IntView(s.MD[c.pos]).minus_r(home, r,false);
		
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
