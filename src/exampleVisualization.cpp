#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include <gecode/gist.hh>
#include <gecode/set.hh>
#include "WinDisplay.cpp"


using namespace Gecode;
using namespace std;

/**
 * \brief %RNA Alignment
 *
 *
 *
 */
class RNAalignment : public Space {

  size_t n;
  size_t m;
  
  size_t undef;

  WinDisplay* wind;

protected:  
  IntVarArray M;
  IntVarArray G;
  SetVarArray H;

  // used as temp variables, maybe they can be removed with better constraint posting

  BoolVarArray b;
  IntVarArray card;
  

public:
  RNAalignment(size_t n1, size_t m1)
    : n(n1),
      m(m1),
      undef(m1+1),
      M(*this,n+1,1,m+1), // assume that undef is m+1! // we only need M_1,...,M_n ==> ignore M_0
      G(*this,n+1,0,m+1), // ignore G_0 as for M
      H(*this,n+1,IntSet::empty,1,m,0), // here, we want H_0(!)...H_n
      b(*this,n+1,0,1),
      card(*this,n+1,0,m+1)
  {
    
    wind=new WinDisplay(m+1,n+1,"Display variables status");
    
    M[0].init(*this,undef,undef);
    G[0].init(*this,undef,undef);
    
    for (int i=0;i<n+1;i++){
      // exacty one M[i], G[i] is undef
      post (*this, tt( imp(M[i]!=undef, G[i]==undef)  ));
      post (*this, tt( imp(G[i]!=undef, M[i]==undef)  ));
      //      post (*this, ff( G[i]==undef &&  M[i]==undef));
      
      // at least one variable has a defined value
      cardinality(*this,H[i],card[i]);
      rel(*this, card[i],IRT_EQ,0,b[i]);      
      post (*this, ff(M[i]==undef && 
		      G[i]==undef && 
		      b[i]));
    }
      

    //sorted M
    for (int i=0;i<n+1;i++){
      for (int j=i+1;j<n+1;j++){
	
	post (*this, tt( imp(M[i]!=undef && 
			     M[j]!=undef,
			     M[i]<M[j])));
      }
    }

    //sorted G
    for (int i=0;i<n+1;i++){
      for (int j=i+1;j<n+1;j++){
	
	post (*this, tt( imp(G[i]!=undef && 
			     G[j]!=undef,
			     G[i]<G[j])));
      }
    }
    
    // a relation between M and G (missing H[i].in(M[i]+1) to be correct)
    for (int i=0;i<n;i++){
      post(*this, tt(
		     imp(M[i]!=undef,
			 M[i+1]==M[i]+1 ||
			 G[i+1]==M[i]
			 )));
    }



    branch(*this, M, INT_VAR_SIZE_MAX, INT_VAL_MED); // suggestion: split largest domain
    branch(*this, G, INT_VAR_SIZE_MAX, INT_VAL_MED);
    branch(*this, H, SET_VAR_SIZE_MAX, SET_VAL_MED_INC); // no set splitting possible
    }

  /// Constructor for cloning \a s
  RNAalignment(bool share, RNAalignment& s) : Space(share,s){
    M.update(*this, share, s.M);
    G.update(*this, share, s.G);
    H.update(*this, share, s.H);
    b.update(*this, share, s.b);
    card.update(*this, share, s.card);
    wind=s.wind;   
    n=s.n;
    m=s.m;
    undef=s.undef;
  }

  /// Perform copying during cloning
  virtual Space*
  copy(bool share) {
    return new RNAalignment(share,*this);
  }

  /// Print solution
  virtual void
  print(std::ostream& os) const {
    os << "Data drawn in windows:\n";
    wind->update(M,G,H);    
  }



};




/** \brief Main-function
 *  \relates RNAalignment
 */
int
main(int argc, char* argv[]) {


  RNAalignment* s = new RNAalignment(10,10);
  Gist::Print<RNAalignment> p("Node explorer");
  Gist::Options o;
  o.inspect.click(&p);

  Gist::dfs(s,o);

  cout << "end\n";
  return 0;
}


