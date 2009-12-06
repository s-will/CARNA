#ifndef NEQ_HH_
#define NEQ_HH_

#include <gecode/kernel.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>

typedef Gecode::BinaryPropagator<Gecode::Int::IntView, Gecode::Int::PC_INT_DOM> GC_BinProp;

class NEq : public GC_BinProp {
private:

protected:
    // variables of superclass
    using GC_BinProp::x0;
    using GC_BinProp::x1;
    
    /// Constructor for cloning \a p
    NEq(Gecode::Space& home, bool share, NEq& p);

    /// Constructor for posting \a p
    NEq(Gecode::Space& home,
	Gecode::Int::IntView x0,
	Gecode::Int::IntView x1);

public:
    /// post a binary neighbor constraint
    static Gecode::ExecStatus post(Gecode::Space& home,
				   Gecode::Int::IntView x0,
				   Gecode::Int::IntView x1);
	    
    /// Copy propagator during cloning
    virtual Gecode::Actor* copy(Gecode::Space& home, bool share);
    
    /// Cost function
    virtual Gecode::PropCost cost(void) const;
    
    /// Perform propagation
    virtual Gecode::ExecStatus propagate(Gecode::Space& home, const Gecode::ModEventDelta&);
};


#endif /*NEQ_HH_*/
