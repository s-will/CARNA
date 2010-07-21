#ifndef STRAL_SCORE_HH
#define STRAL_SCORE_HH


#include <math.h>


#include "matrices.hh"
#include "sequence.hh"
#include "alphabet.hh"
#include "rna_data.hh"

//
//! implements the stral-like scoring function
//!
class StralScore {
    
    typedef std::vector<double> p_vec_t;
    typedef Sequence::size_type size_type;
    
    Sequence seqA;
    Sequence seqB;
    
    p_vec_t p_upA;   // probability paired upstream seq A
    p_vec_t p_downA; // probability paired downstream seq A
    p_vec_t p_unA;   // probability unpaired seq A
    
    p_vec_t p_upB;   // probability paired upstream seq B
    p_vec_t p_downB; // probability paired downstream seq B
    p_vec_t p_unB;   // probability unpaired seq B
    
    const Matrix<double> &sim_mat;
    const Alphabet<char> &alphabet;
    double pf_struct_weight;
    double gap_opening;
    double gap_extension;
    
private:
    void init_prob_vecs(const RnaData &rna,
			p_vec_t &p_up,
			p_vec_t &p_down,
			p_vec_t &p_un);
public:
    
    StralScore(const RnaData &rnaA,
	       const RnaData &rnaB, 
	       const Matrix<double> &sim_mat_,
	       const Alphabet<char> &alphabet_,
	       double pf_struct_weight_,
	       double gap_opening_,
	       double gap_extension_
	       )
	: seqA(rnaA.get_sequence()),
	  seqB(rnaB.get_sequence()),
	  sim_mat(sim_mat_),
	  alphabet(alphabet_),
	  pf_struct_weight(pf_struct_weight_),
	  gap_opening(gap_opening_),
	  gap_extension(gap_extension_)
    {
	// initialize the vectors
	init_prob_vecs(rnaA,p_upA,p_downA,p_unA);
	init_prob_vecs(rnaB,p_upB,p_downB,p_unB);
    }
    
    double sigma(size_type i, size_type j) const {
	// compute the average similarity over all pairs of alignment rows
	//
	// ATTENTION: the treatment of gaps and unknown nucleotide symbols
	// in the aligned alignments is still very AD HOC 
	//
	int pairs=0;
	double seq_score=0;
	for (size_type k=0; k<seqA.get_rows(); k++) {
	    for (size_type l=0; l<seqB.get_rows(); l++) {
		if (alphabet.in(seqA[i][k]) && alphabet.in(seqB[j][l])) {
		    seq_score += sim_mat(alphabet.idx(seqA[i][k]),alphabet.idx(seqB[j][l]));
		    pairs++;
		}
	    }
	}
	if (pairs!=0) seq_score /= pairs;
	
	double res=
	    pf_struct_weight * ( sqrt( p_downA[i]*p_downB[j] ) + sqrt( p_upA[i]*p_upB[j] ) )
	    //    + sqrt( std::max(0.0,p_unA[i]*p_unB[j]) ) * seq_score;
	    + seq_score;
	/* ATTENTION: in the StrAl paper it is claimed that not weighting the sequence score is beneficial,
	   i.e. effectively p_unA[i]==p_unB[j]==1 in above return statement.
	*/
	
	//std::cout << "sigma(" << i << "," << j << ")=" << res << " " << seq_score << std::endl;
	
	return res;
    }
    
    double alpha() const {
	return gap_opening;
    }
    double beta() const {
	return gap_extension;
    }
    
    //! in memoriam :), but it reverses the scoring (i.e. it then scores reverted sequences)
    void 
    rewerse();

};


#endif //STRAL_SCORE_HH
