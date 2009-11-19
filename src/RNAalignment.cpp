#include <gecode/minimodel.hh>
#include <gecode/search.hh>

using namespace Gecode;

/**
 * \brief %RNA Alignment
 *
 *
 *
 */
class RNAalignment : public Space {
protected:
  IntVarArray q;
public:
  RNAalignment(size_t n)
    : q(*this,n,0,n-1) {
    
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
  }

  /// Perform copying during cloning
  virtual Space*
  copy(bool share) {
    return new RNAalignment(share,*this);
  }

  /// Print solution
  virtual void
  print(std::ostream& os) const {
    os << "queens\t";
    for (int i = 0; i < q.size(); i++) {
      os << q[i] << ", ";
      if ((i+1) % 10 == 0)
        os << std::endl << "\t";
    }
    os << std::endl;
  }
};

/** \brief Main-function
 *  \relates RNAalignment
 */
int
main(int argc, char* argv[]) {
  RNAalignment* s = new RNAalignment(10);
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


