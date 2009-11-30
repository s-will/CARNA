#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include <gecode/gist.hh>
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

  WinDisplay* wind;


public:
  RNAalignment(size_t n, WinDisplay* wd)
    : q(*this,n,0,n-1) {

    wind=wd;
    for (size_t i = 0; i<n; i++)
      for (size_t j = i+1; j<n; j++) {
        post(*this, q[i] != q[j]);
        post(*this, q[i]+i != q[j]+j);
        post(*this, q[i]-i != q[j]-j);
      }
    branch(*this, q, INT_VAR_SIZE_MIN, INT_VAL_MIN);
  }

  /// Constructor for cloning \a s
  RNAalignment(bool share, RNAalignment& s) : Space(share,s) {
    q.update(*this, share, s.q);    
    wind=s.wind;    
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


  WinDisplay* d=new WinDisplay(10,10);

  d->display();

  RNAalignment* s = new RNAalignment(10,d);
  Gist::Print<RNAalignment> p("Node explorer");
  Gist::Options o;
  o.inspect.click=&p;

  //  MyInspector i = new MyInspector();

  Gist::dfs(s,o);

  d->close();
  delete d;

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


