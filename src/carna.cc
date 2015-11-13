/*
  CARNA --- Constraint-based Alignment of RNA
  
  A constraint-programming based approach for alignment of RNA with
  structures of unlimited complexity.
  
  Authors: Alessandor Dal Palu, Mathias Moehl, Sebastian Will


  The Carna-algorithm has been published at the conference CP 2010.
  
  Alessandro Dal Palu, Mathias Möhl, and Sebastian Will.  A propagator
  for maximum weight string alignment with arbitrary pairwise
  dependencies. In Proceedings of the 16th International Conference on
  Principles and Practice of Constraint Programming (CP-2010), page 8,
  2010.

*/


// include config.h
#include "../config.h"

#include <limits>

#include "LocARNA/sequence.hh"
#include "LocARNA/scoring.hh"
#include "LocARNA/basepairs.hh"
#include "LocARNA/alignment.hh"
#include "LocARNA/aligner.hh"
#include "LocARNA/rna_data.hh"
#include "LocARNA/arc_matches.hh"
#include "LocARNA/match_probs.hh"
#include "LocARNA/ribosum.hh"
#include "LocARNA/ribofit.hh"
#include "LocARNA/anchor_constraints.hh"
//#include "LocARNA/sequence_annotation.hh"
#include "LocARNA/trace_controller.hh"
#include "LocARNA/ribosum85_60.icc"
#include "LocARNA/global_stopwatch.hh"
#include "LocARNA/pfold_params.hh"
#include "LocARNA/rna_ensemble.hh"

namespace LocARNA {
#include <LocARNA/ribosum85_60.icc>
}



#include <gecode/search.hh>
#include <gecode/minimodel.hh>

#ifdef HAVE_GIST
#include <gecode/gist.hh>
#endif

#include "RnaAlignment.hh"

using namespace std;
//using namespace Gecode;
//using namespace LocARNA;

const std::string
VERSION_STRING = (std::string)PACKAGE_STRING;

// ------------------------------------------------------------
// BEGIN declaration of command line options
// Parameter

double min_prob; // only pairs with a probability of at least min_prob are taken into account

// maximal ratio of number of base pairs divided by sequence
// length. This serves as a second filter on the "significant"
// base pairs.
double max_bps_length_ratio;

int match_score;
int mismatch_score;
int indel_score;
int indel_opening_score;
//int temperature;
int struct_weight;

int tau_factor; // contribution of sequence similarity in an arc match (in percent)

bool opt_no_lonely_pairs; // no lonely pairs option

bool struct_local; // allow exclusions for maximizing alignment of connected substructures
bool sequ_local; // maximize alignment of subsequences

std::string free_endgaps; //!< specification of free end gaps,
// order left end sequence 1, right 1, left 2, right 2
// e.g. "+---" allows free end gaps at the left end of the first alignment string
// ; "----" forbids free end gaps


int max_diff; // maximal difference for positions of alignment edges
// (only used for ends of arcs)
int max_diff_am; //maximal difference between two arc ends, -1 is off

// maximal difference for alignment traces, at arc match
// positions
int max_diff_at_am;

// only consider arc matchs where
//   1. for global (bl-al)>max_diff || (br-ar)<=max_diff    (if max_diff>=0)
//   2. for local (ar-al)-(br-bl)<=max_diff_am              (if max_diff_am>=0)
// except when there is no additional computation of M matrices necessary,
// this occurs if arcs are left-incident with larger arcs where 1 and 2 hold


int exclusion_score; // Score contribution per exclusion
// set to zero for unrestricted structure locality

// double prob_exp_f(int seqlen) {return 1.0/(2*seqlen);} // expected probability of a base pair (null-model)

double exp_prob;
bool opt_exp_prob;

int output_width;

// ------------------------------------------------------------
// File arguments
std::string fileA;
std::string fileB;

std::string clustal_out;
bool opt_clustal_out;

std::string pp_out;
bool opt_pp_out;

bool opt_alifold_consensus_dp; //!< whether to compute consensus dp by alifold

// ------------------------------------------------------------
//
// Options
//
#include <LocARNA/options.hh>


bool opt_help;
bool opt_version;
bool opt_verbose;

bool opt_local_file_output=false;
//bool opt_local_output;
//bool opt_pos_output;

bool opt_write_structure;

// stacking is not implemented in Carna
bool opt_stacking=false;

bool opt_gist;

std::string ribosum_file;
bool use_ribosum;

// use estimated ribosum-like base and arc match scores adapted
// to sequence indentity; overrides use_ribosum and ribosum_file
bool opt_ribofit; 
    
// bool opt_probcons_file;
// std::string probcons_file;

// bool opt_mea_alignment;

// bool opt_write_matchprobs;
// bool opt_read_matchprobs;
// std::string matchprobs_file;

// bool opt_write_arcmatch_scores;

// bool opt_read_arcmatch_scores;
// bool opt_read_arcmatch_probs;
// std::string arcmatch_scores_file;

// int match_prob_method;



//double min_am_prob; // only matched arc-pair with a probability of at least min_am_prob are taken into account
//double min_bm_prob; // only matched base-pair with a probability of at least min_bm_prob are taken into account

// int kbest_k;

std::string seq_constraints_A;
std::string seq_constraints_B;

bool opt_ignore_constraints;

int pf_struct_weight;

int c_d;

int time_limit;

bool opt_lower_bound;
int lower_score_bound;
bool opt_upper_bound;
int upper_score_bound;

// bool opt_mea_gapcost;
// int mea_alpha;
// int mea_beta;
// int mea_gamma;
// int probability_scale;

// bool opt_eval;

LocARNA::option_def my_options[] = {
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Scoring parameters"},

    {"match",'m',0,O_ARG_INT,&match_score,"50","score","Match score"},
    {"mismatch",'M',0,O_ARG_INT,&mismatch_score,"0","score","Mismatch score"},
    {"ribosum-file",0,0,O_ARG_STRING,&ribosum_file,"RIBOSUM85_60","f","Ribosum file"},
    {"use-ribosum",0,0,O_ARG_BOOL,&use_ribosum,"true","bool","Use ribosum scores"},
    {"indel",'i',0,O_ARG_INT,&indel_score,"-350","score","Indel score"},
    {"indel-opening",0,0,O_ARG_INT,&indel_opening_score,"-500","score","Indel opening score"},
    {"struct-weight",'s',0,O_ARG_INT,&struct_weight,"200","score","Maximal weight of 1/2 arc match"},
    {"exp-prob",'e',&opt_exp_prob,O_ARG_DOUBLE,&exp_prob,O_NODEFAULT,"prob","Expected probability"},
    {"tau",'t',0,O_ARG_INT,&tau_factor,"0","factor","Tau factor in percent"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling output"},
    
    {"gist",0,&opt_gist,O_NO_ARG,0,O_NODEFAULT,"",
#ifdef HAVE_GIST
     "Use gist for interactive/graphical search."
#else
     "Use gist for graphical search (feature disabled, recompile to activate)."
#endif
    },

    {"width",'w',0,O_ARG_INT,&output_width,"120","columns","Output width"},
    {"clustal",0,&opt_clustal_out,O_ARG_STRING,&clustal_out,O_NODEFAULT,"file","Clustal output"},
    {"pp",0,&opt_pp_out,O_ARG_STRING,&pp_out,O_NODEFAULT,"file","PP output"},
    // {"local-output",'L',&opt_local_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment"},
    // {"pos-output",'P',&opt_pos_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment positions"},
    
#ifdef HAVE_LIBRNA
    {"alifold-consensus-dp",0,&opt_alifold_consensus_dp,O_NO_ARG,0,O_NODEFAULT,"","Compute consensus dot plot by alifold"},
#endif
    
    {"write-structure",0,&opt_write_structure,O_NO_ARG,0,O_NODEFAULT,"","Write guidance structure in output"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Heuristics for speed accuracy trade off"},

    {"min-prob",'p',0,O_ARG_DOUBLE,&min_prob,"0.01","prob","Minimal probability"},
    {"max-bps-length-ratio",0,0,O_ARG_DOUBLE,&max_bps_length_ratio,"0.0","factor","Maximal ratio of #base pairs divided by sequence length (default: no effect)"},
    {"max-diff-am",'D',0,O_ARG_INT,&max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&max_diff,"-1","diff","Maximal difference for alignment cuts"},
    {"max-diff-at-am",0,0,O_ARG_INT,&max_diff_at_am,"-1","diff","Maximal difference for alignment traces, only at arc match positions"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Constraints"},

    {"noLP",0,&opt_no_lonely_pairs,O_NO_ARG,0,O_NODEFAULT,"","No lonely pairs (only used when predicing ensemble porobabilities and for compatibility with locarna; otherwise no effect)"},
    {"anchorA",0,0,O_ARG_STRING,&seq_constraints_A,"","string","Anchor constraints sequence A."},
    {"anchorB",0,0,O_ARG_STRING,&seq_constraints_B,"","string","Anchor constraints sequence B."},
    {"ignore-constraints",0,&opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},
    
    {"lb",0,&opt_lower_bound,O_ARG_INT,&lower_score_bound,O_NODEFAULT,"score","Lower score bound"},
    {"ub",0,&opt_upper_bound,O_ARG_INT,&upper_score_bound,O_NODEFAULT,"score","Upper score bound"},
    

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling Gecode"},
    {"c_d",0,0,O_ARG_INT,&c_d,"1","distance","Recomputation distance"},
    {"time-limit",0,0,O_ARG_INT,&time_limit,"300000","time","Time limit in ms (always search first solution; turn off by 0)."},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Standard options"},

    {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,"","This help"},
    {"version",'V',&opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","RNA sequences and pair probabilities"},


    {"",0,0,O_SECTION_HIDE,0,O_NODEFAULT,"","Hidden Options"},
    {"ribofit",0,0,O_ARG_BOOL,&opt_ribofit,"false","bool","Use Ribofit base and arc match scores (overrides ribosum)"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Input_files RNA sequences and pair probabilities"},

    {"",0,0,O_ARG_STRING,&fileA,O_NODEFAULT,"file 1","Input file 1  ;; input files can be fasta, dotplot, pp, or aln files"},
    {"",0,0,O_ARG_STRING,&fileB,O_NODEFAULT,"file 2","Input file 2  ;; and can speficy structure and anchor constraints."},
    {"",0,0,0,0,O_NODEFAULT,"",""}
};

// end option declaration
// ------------------------------------------------------------

/** \brief Main-function
 *  \relates RnaAlignment
 */
int
main(int argc, char* argv[]) {
    
    // ----------------------------------------
    // BEGIN process options
    //

    // ------------------------------------------------------------
    // Process options

    bool process_success=process_options(argc,argv,my_options);

    if (opt_help) {
	cout << VERSION_STRING<<endl;

	cout << "A tool for pairwise Alignment of RNA."<<endl<<endl;

	print_help(argv[0],my_options);

	// cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."<<endl<<endl;
	exit(0);
    }

    if (opt_version || opt_verbose) {
	cout << VERSION_STRING<<endl;
	if (opt_version) exit(0); else cout <<endl;
    }


    if (!process_success) {
	std::cerr << "ERROR --- "
		  <<LocARNA::O_error_msg<<std::endl;
	printf("USAGE: ");
	print_usage(argv[0],my_options);
	printf("\n");
	exit(-1);
    }

    if (opt_verbose) {
	print_options(my_options);
    }

    //
    // END process options
    // ----------------------------------------


#ifndef HAVE_GIST
    if (opt_gist) {
	std::cerr << "Graphical output (option --gist) was disabled when compiling." << std::endl
		  << "For enabling this feature, please reconfigure and recompile." << std::endl;
	exit(-1);
    }
#endif


    // options consistency
    if (struct_weight<0) {
	std::cerr << "Structure weight must be greater equal 0."<<std::endl;
	exit(-1);
    }


    // ----------------------------------------
    // BEGIN construct parameter classes from command line options
    //

    // ----------------------------------------  
    // Ribosum matrix
    //
    LocARNA::RibosumFreq *ribosum=NULL;
    LocARNA::Ribofit *ribofit=NULL;
    
    if (opt_ribofit) {
	ribofit = new LocARNA::Ribofit_will2014;
    }
    
    if (use_ribosum) {
	if (ribosum_file == "RIBOSUM85_60") {
	    if (opt_verbose) {
                std::cout <<"Use built-in ribosum."<<std::endl;
	    }
	    ribosum = new LocARNA::Ribosum85_60;
	} else {
	    ribosum = new LocARNA::RibosumFreq(ribosum_file);
	}
    }
    
    // ------------------------------------------------------------
    // Get input data and generate data objects
    //
    LocARNA::PFoldParams pfparams(opt_no_lonely_pairs,false); //no_lonely_pairs,opt_stacking
    
    LocARNA::RnaData *rna_dataA=0;
    try {
	rna_dataA = new LocARNA::RnaData(fileA,min_prob,max_bps_length_ratio,pfparams); 
    } catch (LocARNA::failure &f) {
	std::cerr << "ERROR:\tfailed to read from file "<<fileA <<std::endl
		  << "\t"<< f.what() <<std::endl;
	return -1;
    }
    
    LocARNA::RnaData *rna_dataB=0;
    try {
	rna_dataB = new LocARNA::RnaData(fileB,min_prob,max_bps_length_ratio,pfparams);
    } catch (LocARNA::failure &f) {
	std::cerr << "ERROR: failed to read from file "<<fileB <<std::endl
		  << "       "<< f.what() <<std::endl;
	if (rna_dataA) delete rna_dataA;
	return -1;
    }
    
    const LocARNA::Sequence &seqA=rna_dataA->sequence();
    const LocARNA::Sequence &seqB=rna_dataB->sequence();

    LocARNA::size_type lenA=seqA.length();
    LocARNA::size_type lenB=seqB.length();

    // --------------------
    // handle max_diff restriction

    LocARNA::TraceController trace_controller(seqA,seqB,NULL,max_diff,false);

    // ------------------------------------------------------------
    // Handle constraints (optionally)
    
    LocARNA::AnchorConstraints 
	seq_constraints(lenA,
			seqA.annotation(LocARNA::MultipleAlignment::AnnoType::anchors).single_string(),
			lenB,
			seqB.annotation(LocARNA::MultipleAlignment::AnnoType::anchors).single_string());
    
    if (opt_verbose) {
	if (! seq_constraints.empty()) {
	    std::cout << "Found sequence constraints."<<std::endl;
	}
    }
    
    // ----------------------------------------
    // construct set of relevant arc matches
    //
    LocARNA::ArcMatches *arc_matches=NULL;
    // always initialize from RnaData (reading in arc-matches could be supported later)
    arc_matches = new LocARNA::ArcMatches(*rna_dataA,
					  *rna_dataB,
					  min_prob,
					  (max_diff_am!=-1)
					  ?(LocARNA::size_type)max_diff_am
					  :std::max(lenA,lenB),
					  max_diff_at_am!=-1
					  ? (LocARNA::size_type)max_diff_at_am
					  : std::max(lenA,lenB),
					  trace_controller,
					  seq_constraints
					  );

    LocARNA::BasePairs bpsA = arc_matches->get_base_pairsA();
    LocARNA::BasePairs bpsB = arc_matches->get_base_pairsB();

    // ----------------------------------------
    // report on input in verbose mode
    if (opt_verbose) {
	std::cout << "Sequence A: "<<std::endl;
	seqA.write(cout);
	std::cout<<" (Length:"<< seqA.length()<<", Basepairs:"<<bpsA.num_bps() << ")" <<std::endl;

	std::cout << "Sequence B: "<<std::endl;
	seqB.write(cout);
	std::cout<<" (Length:"<< seqB.length()<<", Basepairs:"<<bpsB.num_bps() << ")" <<std::endl;

	cout <<std::endl
	     <<"Base Pair Matches: "<<arc_matches->num_arc_matches() << "." <<std::endl;
	// cout << "Base Identity: "<<(seq_identity(seqA,seqB)*100)<<endl;
    }

    // ----------------------------------------
    // construct scoring

    double my_exp_probA = opt_exp_prob?exp_prob:LocARNA::prob_exp_f(lenA);
    double my_exp_probB = opt_exp_prob?exp_prob:LocARNA::prob_exp_f(lenB);


    // ----------------------------------------
    // Scoring Parameter
    //
    LocARNA::ScoringParams 
	scoring_params(match_score,
		       mismatch_score,
		       indel_score,
		       indel_score, // = indel loop score (not used by CARNA)
		       indel_opening_score,
		       indel_opening_score, // = indel loop opening score (not used by CARNA)
                       ribosum,
                       ribofit,
		       0, // unpaired penalty (not configurable in CARNA)
		       struct_weight,
		       tau_factor,
		       exclusion_score,
		       my_exp_probA,
		       my_exp_probB,
		       0, // temperature
		       false, // opt_stacking not implemented in Carna
		       false, // opt_new_stacking not implemented in Carna
		       false, // mea
		       0,0,0,0 // alpha, beta,gamma,probabilit_scale (unused in Carna)
		       );

    LocARNA::Scoring scoring(seqA,
		    seqB,
		    *rna_dataA,
		    *rna_dataB,
		    *arc_matches,
		    0L,
		    scoring_params,
		    false // no Boltzmann weights
		    );    

    // ------------------------------------------------------------
    // construct the constraint model / root space
    //
    
    const RnaAlignmentParams &rna_alignment_params =
	RnaAlignment::create
	(
	 LocARNA::Aligner::create()
	 . seqA(seqA)
	 . seqB(seqB)
	 . arc_matches(*arc_matches)
	 . scoring(scoring)
	 . no_lonely_pairs(false) // no_lonely_pairs does not make sense in carna
	 . struct_local(struct_local)
	 . sequ_local(sequ_local)
	 . free_endgaps(free_endgaps)
	 . max_diff_am(max_diff_am)
	 . max_diff_at_am(max_diff_at_am)
	 . trace_controller(trace_controller)
	 . min_am_prob(0)
	 . min_bm_prob(0)
	 . stacking(false)
	 . constraints(seq_constraints)
	 )
	. lower_score_bound(opt_lower_bound?lower_score_bound:Gecode::Int::Limits::min)
	. upper_score_bound(opt_upper_bound?upper_score_bound:Gecode::Int::Limits::max)
	. gist(opt_gist)
	. output_width(output_width)
	;
        
    RnaAlignment *rna_alignment_space = 
	new RnaAlignment(rna_alignment_params);

    // ------------------------------------------------------------
    // run the search engine
    //
    if (opt_gist) {
#ifdef HAVE_GIST
	Gecode::Gist::Print<RnaAlignment> p("Node explorer");
	Gecode::Gist::Options o;
	o.inspect.click(&p);

	o.c_d =  c_d;

	//Gecode::Gist::dfs(s,o);
	Gecode::Gist::bab(rna_alignment_space,o);
#endif
    } else {
	Gecode::Search::Options o;

	o.c_d =  c_d;
	o.stop=0L;
	
	if (time_limit>0) {
	    Gecode::Search::Stop *timestop=0L;
	    // activate next line for "if time limit is exceeded, no solution is produced"
	    // timestop = new Gecode::Search::TimeStop(time_limit);
	    // activate instead for "the first solution is always produced"
	    timestop = new Gecode::Search::TimeStop(std::numeric_limits<unsigned long int>::max());
	    o.stop=timestop;
	}
	
	// construct engine
	Gecode::BAB<RnaAlignment> e(rna_alignment_space,o);
	
	bool first_solution=true;

	// ----------------------------------------
	// enumerate solutions
	RnaAlignment* ex=NULL;
	while ((ex = e.next()) && (ex != NULL)) {
	    // write each solution to stdout and optionally to files
	    // in clustal and pp format
	    
	    if (first_solution && time_limit>0) {
		((Gecode::Search::TimeStop *)o.stop)->limit(time_limit);
		first_solution=false;
	    }
	    
	    ex->print(std::cout);

	    if (!clustal_out.empty()){
		ofstream outfile;
		outfile.open(clustal_out.c_str(),ios::out | ios::trunc);
		if (outfile.good()) {
		    
		    assert(ex->all_assigned());
		    LocARNA::Alignment alignment=ex->to_alignment();
			
		    LocARNA::MultipleAlignment ma(alignment,opt_local_file_output);
		    
		    outfile << "CLUSTAL W --- "<<PACKAGE_STRING;
		    
		    // for legacy, clustal files of pairwise alignments contain the score 
		    if (seqA.num_of_rows()==1 && seqB.num_of_rows()==1)
			outfile  <<" --- Score: " << ex->score();
		    outfile <<std::endl<<std::endl;
		    
		    ma.write(outfile,output_width);
		    outfile.close();
		} else {
		    std::cerr << "Cannot write solution to file "<<clustal_out<<std::endl;
		}
	    }

	    if (!pp_out.empty()){
		ofstream outfile;
		outfile.open(pp_out.c_str(),ios::out | ios::trunc);
		
		if (outfile.good()) {
		    assert(ex->all_assigned());
		    
		    LocARNA::Alignment alignment=ex->to_alignment();

		    if (opt_alifold_consensus_dp) {
			LocARNA::PFoldParams pfparams(opt_no_lonely_pairs,false); // opt_stacking
			
			LocARNA::MultipleAlignment ma(alignment,opt_local_file_output);
			LocARNA::RnaEnsemble ens(ma,pfparams,false,true); // alifold the alignment
			// construct consensus rna data from ensemble
			LocARNA::RnaData consensus(ens,min_prob,max_bps_length_ratio,pfparams);
			consensus.write_pp(outfile); // write alifold dot plot
		    } else {
			// compute averaged consensus base pair probabilities
			LocARNA::RnaData consensus(*rna_dataA,
						   *rna_dataB,
						   alignment,
						   my_exp_probA,
						   my_exp_probB,
						   opt_local_file_output);
			consensus.write_pp(outfile); // write averaged dot plot
		    }
		    
		    outfile.close();
		} else {
		    std::cerr << "Cannot write solution to file "<<pp_out<<std::endl;
		}
	    }
	    if (ex) delete ex;
	}
	

	Gecode::Search::Statistics stats = e.statistics();

	cout << "Nodes:    "<<stats.node << std::endl
	     << "Fail:     "<<stats.fail << std::endl
	     << "Depth:    "<<stats.depth << std::endl
	    //<< "Restarts: "<<stats.restart << std::endl
	    //<< "Nogoods:  "<<stats.nogood
	     << std::endl;

	if (e.stopped()) {
	    cout << "Time limit exceeded."<<std::endl;
	}

	if (o.stop) delete o.stop;

    }
    
    if (rna_alignment_space) delete rna_alignment_space;
    if (arc_matches) delete arc_matches;
    if (ribosum) delete ribosum;
    return 0;
}
