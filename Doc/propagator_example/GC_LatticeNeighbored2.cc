/*
 *  Main authors:
 *     Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 *
 *  Contributing authors:
 *     Sebastian Will http://www.bioinf.uni-freiburg.de/~will/
 *
 *  Copyright:
 *     Martin Mann, 2007
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#include "cpsp/GC_LatticeNeighbored2.hh"
#include "cpsp/GC_StlSetRangeIterator.hh"

#include <biu/assertbiu.hh>

#ifdef TMPOUT
	#include <iostream>
#endif




  GC_LatticeNeighbored2::GC_LatticeNeighbored2(	Gecode::Space* home,
  								Gecode::Int::IntView x0,
  								Gecode::Int::IntView x1,
  								const biu::LatticeFrame * lattice_)
    : GC_BinProp(*home,x0,x1),
    	lattice(lattice_)
  {


  }

  GC_LatticeNeighbored2::GC_LatticeNeighbored2(	Gecode::Space* home,
  								bool share,
  								GC_LatticeNeighbored2& p) :
    	GC_BinProp(*home,share,p),
    	lattice(p.lattice)
  {}

  Gecode::ExecStatus
  GC_LatticeNeighbored2::post(	Gecode::Space* home,
  						Gecode::Int::IntView x0, Gecode::Int::IntView x1,
  						const biu::LatticeFrame * lattice)
  {
	if (same(x0,x1)) {
		return Gecode::ES_FAILED;
    } else {
		new (*home) GC_LatticeNeighbored2( home,x0,x1,lattice);
    }
    return Gecode::ES_OK;
  }


  Gecode::Actor*
  GC_LatticeNeighbored2::copy(Gecode::Space& home, bool share) {
    return new (home) GC_LatticeNeighbored2(&home,share,*this);
  }

  Gecode::PropCost
  GC_LatticeNeighbored2::cost(void) const {
    if (x0.assigned() || x1.assigned()) {
		return Gecode::PropCost::unary(Gecode::PropCost::HI);
    }
    return Gecode::PropCost::binary(Gecode::PropCost::HI);
  }

  Gecode::ExecStatus
  GC_LatticeNeighbored2::propagate(Gecode::Space& home, const Gecode::ModEventDelta&) {
  	//std::cerr <<"vorher : "<<x0 <<" / " <<x1 <<std::endl;
		// pruefe ob minimale propagationsschranke erreicht
//	if (x0.size() < minPropSize && x1.size() < minPropSize)
//		return Gecode::ES_SUBSUMED;


  	Gecode::ModEvent ret = Gecode::ME_GEN_NONE;
	if (x0.assigned()) {
    	ret = removeNonNeighbors(&home,x1,x0.val());
    }else if (x1.assigned()) {
    	ret = removeNonNeighbors(&home,x0,x1.val());
    } else {
    	ret = removeNonNeighbors(&home);
    }
  	//std::cerr <<"danach : "<<x0 <<" / " <<x1 <<std::endl;
    	// pruefen auf ME_GEN_FAILED
    GECODE_ME_CHECK(ret);


	if (x0.assigned() || x1.assigned()) {	// no further propagation possible
		return Gecode::ES_SUBSUMED(*this, home);
	}
	// else both are unassigned
	return Gecode::ES_FIX;
//    return (ret ==	Gecode::ME_GEN_ASSIGNED
//    			?	Gecode::ES_SUBSUMED(*this, sizeof(*this))
//    			:	Gecode::ES_FIX);
  }

	Gecode::ModEvent
	GC_LatticeNeighbored2::removeNonNeighbors(	Gecode::Space* home,
										Gecode::Int::IntView& x,
										int center)
	{
		biu::LatticeFrame::index_set data;
		data = getIndexedNeighbors(center,data);
		GC_StlSetRangeIterator dataIt(&data);
		return x.inter_r(*home, dataIt);
	}

/* DDS SPECIFIC BUT CURRENTLY NOT SUPPORTED IN GECODE 1.3.1
	Gecode::VarIter*
	GC_LatticeNeighbored2::vars(void) const {
		if (cancelProp || x0.assigned() || x1.assigned())
			return new Gecode::EmptyVarIter();
		else
	      return new Gecode::BinaryVarIter<Gecode::Int::IntView,Gecode::Int::IntView>(x0,x1);

	}
*/	

	Gecode::ModEvent
	GC_LatticeNeighbored2::removeNonNeighbors(Gecode::Space* home) {

		const biu::LatticeFrame::index_set& neighborhood(lattice->getIndexedNeighborhood());

			// a wird referenz auf kleinere domain
		Gecode::Int::IntView* a = &x0, *b = &x1;
		if(a->size() > b->size()) {
			a = &x1, b = &x0;
		}


			// alte groessen der variablen
		unsigned int aSize = a->size(), bSize = b->size();

			// aenderungsstatus der variablen
		Gecode::ModEvent	aStat = Gecode::ME_GEN_NONE,
							bStat = Gecode::ME_GEN_NONE;

			// a domain einschraenken und neue b domain (newDom) aufbauen
		Gecode::Int::ViewValues<Gecode::Int::IntView> aIt(*a);
		biu::LatticeFrame::index_set newDom, toRem;
		biu::LatticeFrame::index_set::iterator whereToPut = toRem.begin();
		biu::LatticeFrame::index_set::const_iterator neigh;
		bool noNeighInB = true;
		while(aIt()) {	// laufe durch alle elemente von a
				// init
			noNeighInB = true;
				// nachbarn testen
			for(neigh = neighborhood.begin(); neigh != neighborhood.end(); neigh++) {
				if (b->in( *neigh +aIt.val())) {	// dann nachbar in b vorhanden
					noNeighInB = false;	// merken dass nachbar gefunden
					newDom.insert(*neigh+aIt.val());	// nachbar in neue domain aufnehmen
				}
			}
				// altes element loeschen wenn noetig
			if (noNeighInB) {
					// zu entfernendes aVal aus a merken, wenn kein nachbar in b
				whereToPut = toRem.insert(whereToPut, aIt.val());
			}
				// naechstes element via iterator
			++aIt;
		}
			// alle elemente in a loeschen, die keine nachbarn in b haben
		GC_StlSetRangeIterator remDomIt(&toRem);
		aStat = a->minus_r(*home,remDomIt);

			// b domain angleichen
		if (aStat != Gecode::ME_GEN_FAILED) {	// nur wenn a nicht leer
			GC_StlSetRangeIterator newDomIt(&newDom);
			bStat = b->inter_r(*home,newDomIt);
		} else {
			return Gecode::ME_GEN_FAILED;	// a domain leer
		}

			// rueckgabewert bestimmen
		if ( bStat == Gecode::ME_GEN_FAILED ) {
			return Gecode::ME_GEN_FAILED;	// dann b domain leer geworden
		}
		if ( aSize == a->size() && bSize == b->size())
			return Gecode::ME_GEN_NONE;		// nixs passiert
		if ( aStat == Gecode::ME_GEN_ASSIGNED
				|| bStat == Gecode::ME_GEN_ASSIGNED)
			return Gecode::ME_GEN_ASSIGNED;
		return Gecode::Int::ME_INT_DOM;	// werte haben sich geaendert
	}

	biu::LatticeFrame::index_set&
	GC_LatticeNeighbored2::getIndexedNeighbors(
			const biu::LatticeFrame::index_type center, biu::LatticeFrame::index_set& toFill) const
	{
		const biu::LatticeFrame::index_set& neighborhood(lattice->getIndexedNeighborhood());

		toFill.clear();
			// alle nachbarn generieren
		for(biu::LatticeFrame::index_set::const_iterator it = neighborhood.begin();
				it != neighborhood.end(); it++) {
			toFill.insert(center+(*it));
		}
		return toFill;
	}


