#include <gecode/minimodel.hh>
#include <gecode/search.hh>

// #include "alignmentScorePropagator/scoreAlignment.hh"

using namespace Gecode;

/**
 * \brief %RNA Alignment
 *
 *
 *
 */
class RNAalignment : public Space {
    const size_t n;
    const size_t m;

    const size_t undef;

protected:
    IntVarArray M;
    IntVarArray G;
    SetVarArray H;
    
    IntVar Score;

public:
    RNAalignment(const RnaData &rna_data_R,const RnaData &rna_data_S,AlignmentParams &alignment_params)
	:
	n(rna_data_R.size()),
	m(rna_data_R.size()),
	undef(m+1),
	M(*this,n+1,1,m+1), // assume that undef is m+1! // we only need M_1,...,M_n ==> ignore M_0
	G(*this,n+1,0,m+1), // ignore G_0 as for M
	H(*this,n+1,IntSet::empty,1,m,0) // here, we want H_0(!)...H_n
    {
	//ignore M_0, G_0
	M[0]=undef;
	G[0]=undef;
		
	scoreAlignment::post(*this,rna_data_R,rna_data_S,alignment_params,M,G,H,Score);
	
	// presumably reasonable: first enumerate M, then the rest
	branch(*this, M, INT_VAR_SIZE_MAX, INT_VAL_MED); // suggestion: split largest domain
	branch(*this, G, INT_VAR_SIZE_MAX, INT_VAL_MED);

	// we should show: here, the H variables are already fixed
	// (this still depends on the propagator, otw. fix Hs to upper bound i sufficient)
	// Then: branching over H should be dropped/replaced
	branch(*this, H, SET_VAR_SIZE_MAX, SET_VAL_MED_INC); // no set splitting possible
    }

    /// Constructor for cloning \a s
    RNAalignment(bool share, RNAalignment& s) : Space(share,s) {
	M.update(*this, share, s.M);
	G.update(*this, share, s.G);
	H.update(*this, share, s.H);
    }

    /// Perform copying during cloning
    virtual Space*
    copy(bool share) {
	return new RNAalignment(share,*this);
    }

    /// Print solution
    virtual void
    print(std::ostream& out) const {
	out << "RNAalignment::print: printing not supported yet."
    }
};

/** \brief Main-function
 *  \relates RNAalignment
 */
int
main(int argc, char* argv[]) {
    
    /*
      Get RnaData objects and AlignmentParams similar to LocARNA
    */
    
    RNAalignment* s = new RNAalignment(rna_data_R,rna_data_S,alignment_params);
    DFS<RNAalignment> e(s);
    
    RNAalignment* ex = e.next();
    if (ex != NULL) {
	ex->print(std::cout);
	delete ex;
    } else {
	std::cout << "No solution" << std::endl;
    }
    
    return 0;
}
