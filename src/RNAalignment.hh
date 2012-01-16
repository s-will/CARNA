#ifndef RNA_ALIGNMENT
#define RNA_ALIGNMENT

#include "config.h"

#include <limits>

#include <fstream>
#include <iostream>

#ifdef HAVE_GIST
#include "WinHandler.hh"
#endif // HAVE_GIST

#include "alignment_score.hh"

#include "locarna.hh"

// EXPERIMENTAL: therefore only global var
//const size_t discrepancy_limit=7;

const bool custom_branching=true;


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
	
	//! \brief selected domain values
	//! @note An int vector is used for constructing the choice of domain values
	//! by the AlignmentScore propagator and communicating this to the
	//! space.  This is not optimized for speed! One should
	//! better use ranges instead of single values and use specific Gecode
	//! mechanisms for doing this. Using the vector is a workaround, since
	//! I (SW) could not make Gecode's BitSet working (or what else is more
	//! appropriate?).
	std::vector<int> values;
	
	//! signal whether the currenty best trace yields a new lower
	//! bound
	bool new_lower_bound;
	std::vector<size_t> best_traceA;
	std::vector<size_t> best_traceB;
	LocARNA::score_t best_trace_score;
	
	ChoiceData()
	    : enum_M(false), pos(0),val(0), values(),
	      new_lower_bound(false),
	      best_traceA(),
	      best_traceB(),
	      best_trace_score(std::numeric_limits<LocARNA::score_t>::min())
	{}
	
    };
    
protected:
   
#ifdef HAVE_GIST
    WinHandler* wind;
#endif

    
    //! MD[i] is position of match or deletion in row i. We model only
    //! match and deletion positions explicitely. Insertion positions are
    //! derived! Compare the handling of MD and M by the propagator AlignmentScore. 
    Gecode::IntVarArray MD;  
    Gecode::BoolVarArray M; //!< M[i] is true iff i~MD[i] is a match
    
    Gecode::IntVar Score;

    //size_t discrepancy; // experimental for simulating LDS
    
    //! the propagator AlignmentScore selects the choice for the
    //! brancher RNAalignBranch and communicates its choice in this
    //! field of the space
    ChoiceData choice_data;

 
public:
    /** 
     * Set up the constraint model
     * 
     * @param seqA_ sequence A 
     * @param seqB_ sequence B
     * @param arcmatches_ the arc matches
     * @param aligner_params parameters for the aligner
     * @param scoring scoring object
     * @param lower_score_bound lower bound on the score (use Gecode::Int::Limits::min for no bound)
     * @param upper_score_bound upper bound on the score (use Gecode::Int::Limits::max for no bound)
     * @param opt_graphical_output support of gist graphical output
     * 
     */
    RNAalignment(const LocARNA::Sequence &seqA_, const LocARNA::Sequence &seqB_,
		 const LocARNA::ArcMatches &arcmatches_,
		 const LocARNA::AlignerParams &aligner_params, 
		 const LocARNA::Scoring &scoring,
		 int lower_score_bound,
		 int upper_score_bound,		 
		 bool opt_graphical_output);

    //! \brief Constructor for cloning s
    RNAalignment(bool share, RNAalignment& s) : 
	Gecode::Space(share,s),
	seqA(s.seqA),
	seqB(s.seqB),
	arcmatches(s.arcmatches),
	n(s.n),
	m(s.m),
#ifdef HAVE_GIST
	wind(s.wind),
#endif
	choice_data(s.choice_data)
	//, discrepancy(s.discrepancy)
    {
	MD.update(*this, share, s.MD);
	M.update(*this, share, s.M);
	Score.update(*this,share,s.Score);
    }
    
    //! \brief Generate a new copy of the space
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
    print_clustal_format(std::ostream& out_s, size_t output_width) const;

    //! print solution in pp format
    void 
    print_pp_format(std::ostream& out_s,
		    const LocARNA::BasePairs& bpsA, const LocARNA::BasePairs& bpsB, 
		    const LocARNA::Scoring& scoring, 
		    const LocARNA::AnchorConstraints& seq_constraints,
		    const size_t output_width,
		    bool alifold_consensus_dp
		    ) const;

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
	
	/** 
	 * \brief Fix variables according to trace
	 * 
	 * Assigns the variables in the space such that the space is
	 * solved and its solution corresponds to the trace
	 *
	 * @param home the space
	 * @param traceA trace vector for A
	 * @param traceB trace vector for B
	 * 
	 * @return gecode exec status
	 */
	Gecode::ModEvent
	fix_vars_from_trace(Gecode::Space& home,
			    const std::vector<size_t> &traceA,
			    const std::vector<size_t> &traceB) const;

    public:

	// return true, if there are unassigned vars left for this brancher
	//        false, otherwise
	virtual bool status(const Gecode::Space& home) const {
	    const RNAalignment& s = static_cast<const RNAalignment&>(home);
	    
	    for (size_t i=start; i<(size_t)s.MD.size(); i++) {
		if (! s.MD[i].assigned()) {start=i;return true;}
	    }
	    return false;
	}
	
	// returns choice
	virtual Gecode::Choice* choice(Gecode::Space& home) {
	    const RNAalignment& s = static_cast<const RNAalignment&>(home);
	
	    // the decision about choice is moved to the propagator, which
	    // returned it's choice in the field choice of the space s
	    
	    //std::cout << "CHOICE "<<s.pos<<" "<<s.val<<" "<<s.MD[s.pos]<<std::endl;
	    return new Choice(*this,s.choice_data);
	}
	
	/// commits the choice c and alternative a 
	virtual Gecode::ExecStatus commit(Gecode::Space& home, 
			     const Gecode::Choice& _c,
			     unsigned int a);
    
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
