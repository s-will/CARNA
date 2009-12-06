#include "NEq.hh"
//#include <iostream>


NEq::NEq(Gecode::Space& home,
	 Gecode::Int::IntView x0,
	 Gecode::Int::IntView x1)
    : GC_BinProp(home,x0,x1)
{}

NEq::NEq(Gecode::Space& home,
	 bool share,
	 NEq& p) :
    GC_BinProp(home,share,p)
{}

Gecode::ExecStatus
NEq::post(Gecode::Space& home,
	  Gecode::Int::IntView x0, Gecode::Int::IntView x1)
{
    if (same(x0,x1)) {
	return Gecode::ES_FAILED;
    } else {
	new (home) NEq(home,x0,x1);
    }
    return Gecode::ES_OK;
}

Gecode::Actor*
NEq::copy(Gecode::Space& home, bool share) {
    return new (home) NEq(home,share,*this);
}

Gecode::PropCost
NEq::cost(void) const {
    if (x0.assigned() || x1.assigned()) {
	return Gecode::PropCost::unary(Gecode::PropCost::LO);
    }
    return Gecode::PropCost::binary(Gecode::PropCost::LO);
}

Gecode::ExecStatus
NEq::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {

    Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
    if (x0.assigned()) {
    	ret = x1.nq(home,x0.val());
    } else if (x1.assigned()) {
    	ret = x0.nq(home,x1.val());
    }
    
    GECODE_ME_CHECK(ret);
    
    if (x0.assigned() || x1.assigned()) {	// no further propagation possible
	return Gecode::ES_SUBSUMED(*this, home);
    }
    
    // else both are unassigned
    return Gecode::ES_FIX;
}
