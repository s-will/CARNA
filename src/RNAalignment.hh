#ifndef RNA_ALIGNMENT
#define RNA_ALIGNMENT


#include "alignment_score.hh"




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

    const size_t n;
    const size_t m;
    
    const size_t undef;
    
    size_t pos;
    size_t val;
    size_t minval;
    size_t maxval;

    WinDisplay* wind;
    
    
protected:
    Gecode::IntVarArray M;
    Gecode::IntVarArray G;
    Gecode::SetVarArray H;
    
    Gecode::IntVar Score;

public:
    RNAalignment(const Sequence &seqA, const Sequence &seqB, const ArcMatches &arcmatches,const AlignerParams &aligner_params, const Scoring &scoring)
	:
	n(seqA.length()),
	m(seqB.length()),
	undef(m+1),
	M(*this,n+1,1,m+1), // assume that undef is m+1! // we only need M_1,...,M_n ==> ignore M_0
	G(*this,n+1,0,m+1), // ignore 0 as done for M
	H(*this,n+1,Gecode::IntSet::empty,1,m,0), // here, we want H_0(!)...H_n, init such that empty sets allowed, maximal sets are {1..m}
	Score(*this,Gecode::Int::Limits::min,Gecode::Int::Limits::max)
    {
      wind=new WinDisplay(n+1,m+1,"Display variables status");
      
      //ignore M_0
      rel(*this,M[0],Gecode::IRT_EQ,undef);
      rel(*this,G[0],Gecode::IRT_EQ,undef);

	AlignmentScore::post(*this,seqA,seqB,arcmatches,aligner_params,scoring,
			     M,G,H,Score);
	
	
	//  first enumerate M, then the rest
	
	// resort M vector, such that we start enumerating in the middle
	// Ideally sort like balanced binary tree in array
	Gecode::IntVarArgs M_resorted(M.size());
        for (size_type i=0; i<(size_type)M.size(); i++) { 
	    M_resorted[i] = M[(i+M.size()/2)%M.size()];
	}

	// test the case where a good score is known
	// rel(*this,Score,IRT_GR, 4000);
	
	if (custom_branching) {
	    RNAalignBrancher::post(*this);
	} else {
	    Gecode::branch(*this, M_resorted, Gecode::INT_VAR_SIZE_MAX, Gecode::INT_VAL_MED);
	}
	
	Gecode::branch(*this, G, Gecode::INT_VAR_SIZE_MAX, Gecode::INT_VAL_MED);
	
	// we should show: here, the H variables are already fixed
	// (this still depends on the propagator, otw. fix Hs to upper bound i sufficient)
	// Then: branching over H should be dropped/replaced
	Gecode::branch(*this, H, Gecode::SET_VAR_SIZE_MAX, Gecode::SET_VAL_MED_INC); // no set splitting possible
	
    }

    /// Constructor for cloning \a s
    RNAalignment(bool share, RNAalignment& s) : 
	Gecode::Space(share,s),
	n(s.n),
	m(s.m),
	undef(s.undef),
	pos(s.pos),
	val(s.val),
	minval(s.minval),
	maxval(s.maxval)
    {
	M.update(*this, share, s.M);
	G.update(*this, share, s.G);
	H.update(*this, share, s.H);
	Score.update(*this,share,s.Score);
	wind=s.wind;   
    }

    /// Perform copying during cloning
    virtual Gecode::Space*
    copy(bool share) {
	return new RNAalignment(share,*this);
    }
    
    /// Print solution
    virtual void
    print(std::ostream& out) const {
	std::cout << "SOLUTION" << std::endl;
	//std::cout << "Matches:    " << M << std::endl;
	//std::cout << "Deletions:  " << G << std::endl;
	//std::cout << "Insertions: " << H << std::endl;
	std::cout << "Score:      " << Score << std::endl;
	
	wind->update(M,G,H);
    }

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
	size_t pos;
	size_t val;
	size_t minval;
	size_t maxval;
	

	Choice(const Brancher& b,size_t pos_, size_t val_,size_t minval_,size_t maxval_)
	    : Gecode::Choice(b,2),
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
	
	for (size_t i=start; i<(size_type)s.M.size(); i++) {
	    if (! s.M[i].assigned()) {start=i;return true;}
	}
	return false;
    }
    
    // returns choice
    virtual Gecode::Choice* choice(Gecode::Space& home) {
	const RNAalignment& s = static_cast<const RNAalignment&>(home);
	
	// move decision about choice to the propagator
	
	std::cout << "CHOICE "<<s.pos<<" "<<s.val<<" "<<s.minval<<" "<<s.maxval<<" "<<s.M[s.pos]<<std::endl;
	return new Choice(*this,s.pos,s.val,s.minval,s.maxval);
    }
    
    // commits the choice c and alternative a 
    virtual Gecode::ExecStatus commit(Gecode::Space& home, 
				      const Gecode::Choice& _c,
				      unsigned int a) {
	const RNAalignment& s = static_cast<const RNAalignment&>(home);
	const Choice& c = static_cast<const Choice&>(_c);
	
	/*
	// split in eq and nq
	if (a==0) {
	    return Gecode::me_failed(Gecode::Int::IntView(s.M[c.pos]).eq(home, (int)c.val))
		? Gecode::ES_FAILED
		: Gecode::ES_OK;
	} else {
	    return Gecode::me_failed(Gecode::Int::IntView(s.M[c.pos]).nq(home, (int)c.val))
		? Gecode::ES_FAILED
		: Gecode::ES_OK;
		}
	*/
	if (a==0) {
	    // ACHTUNG: undef values!!!
	    return Gecode::me_failed(
				     Gecode::Int::IntView(s.M[c.pos]).gq(home, (int)c.minval)
				     |
				     Gecode::Int::IntView(s.M[c.pos]).lq(home, (int)c.maxval)
				     )
		? Gecode::ES_FAILED
		: Gecode::ES_OK;
	} else {
	    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
	    
	    for (size_t val=c.minval; val<=c.maxval; val++)
		ret|=Gecode::Int::IntView(s.M[c.pos]).nq(home, (int)val);
	    
	    return Gecode::me_failed(ret)
		? Gecode::ES_FAILED
		: Gecode::ES_OK;
	}
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
