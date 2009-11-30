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
public:
  IntVarArray q;
  SetVarArray H;

  WinDisplay* wind;
  WinDisplay* wind2;


public:
  RNAalignment(size_t n, WinDisplay* wd, WinDisplay* wd2)
    : q(*this,n,0,n-1),
      H(*this, n,IntSet::empty,1,n) // here, we want H_0(!)...H_n
    {

    wind=wd;
    wind2=wd2;
    for (size_t i = 0; i<n; i++)
      for (size_t j = i+1; j<n; j++) {
        post(*this, q[i] != q[j]);
        post(*this, q[i]+i != q[j]+j);
        post(*this, q[i]-i != q[j]-j);
	dom(*this, H[i], SRT_DISJ, IntSet(0,i-2));
	dom(*this, H[i], SRT_SUP, i+2, n);
      }
    branch(*this, q, INT_VAR_SIZE_MIN, INT_VAL_MIN);
  }

  /// Constructor for cloning \a s
  RNAalignment(bool share, RNAalignment& s) : Space(share,s) {
    q.update(*this, share, s.q);
    H.update(*this,share,s.H);    
    wind=s.wind;    
    wind2=s.wind2;
  }

  /// Perform copying during cloning
  virtual Space*
  copy(bool share) {
    return new RNAalignment(share,*this);
  }

  /// Print solution
  virtual void
  print(std::ostream& os) const {
    os << "queens:\t";
    for (int i = 0; i < q.size(); i++) {
      os << q[i] << ", ";
      if ((i+1) % 10 == 0)
        os << std::endl << "\t";
    }
    os << std::endl;

    wind->update(q);    
    wind2->update(H);    

    //    while (!(main_disp.is_closed)){
    //  main_disp.wait(20);
    //}

  }



};




/** \brief Main-function
 *  \relates RNAalignment
 */
int
main(int argc, char* argv[]) {


  WinDisplay* d=new WinDisplay(10,10,"Display int var array");
  WinDisplay* d2=new WinDisplay(10,10,"Display set var array");

  //  d->display();

  RNAalignment* s = new RNAalignment(10,d,d2);
  Gist::Print<RNAalignment> p("Node explorer");
  Gist::Options o;
  o.inspect.click=&p;

  //  MyInspector i = new MyInspector();

  Gist::dfs(s,o);

  d->close();
  d2->close();

  cout << "end\n";
  return 0;
/*  
  DFS<RNAalignment> e(s);
  
  RNAalignment* ex = e.next();
  if (ex != NULL) {
    ex->print(std::cout);
    delete ex;
  } else {
    std::cout << "No solution" << std::endl;
  }


  return 0;
*/
}


