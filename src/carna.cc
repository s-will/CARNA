/*
  CARNA --- Constraint-based Alignment of RNA
  
  A constraint-programming based approach for alignment of RNA with
  structures of unlimited complexity.
  
  Authors: Alessandor Dal Palu, Mathias Moehl, Sebastian Will


  The Carna-algorithm was published at the conference CP 2010.
  
  Alessandro Dal Palu, Mathias Möhl, and Sebastian Will.  A propagator
  for maximum weight string alignment with arbitrary pairwise
  dependencies. In Proceedings of the 16th International Conference on
  Principles and Practice of Constraint Programming (CP-2010), page 8,
  2010.

*/


// include config.h
#include "../config.h"

#include <limits>

// get all locarna headers
#include <locarna.hh>

#include <LocARNA/aux.hh>

namespace LocARNA {
#include <LocARNA/ribosum85_60.icc>
}



#include <gecode/search.hh>
#include <gecode/minimodel.hh>

#ifdef HAVE_GIST
#include <gecode/gist.hh>
#endif

#include "RNAalignment.hh"

using namespace std;
//using namespace Gecode;
//using namespace LocARNA;

const std::string
VERSION_STRING = (std::string)PACKAGE_STRING;

// ------------------------------------------------------------
// BEGIN declaration of command line options
// Parameter

double min_prob; // only pairs with a probability of at least min_prob are taken into account

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
std::string file1;
std::string file2;

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

//bool opt_local_output;
//bool opt_pos_output;

bool opt_write_structure;

// stacking is not implemented in Carna
// bool opt_stacking;

bool opt_gist;

std::string ribosum_file;
bool use_ribosum;

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

bool opt_time_limit;
int time_limit;

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
    // {"exclusion",'E',0,O_ARG_INT,&exclusion_score,"0","score","Exclusion weight"},
    //{"stacking",0,&opt_stacking,O_NO_ARG,0,O_NODEFAULT,"","Use stacking terms (needs stack-probs by RNAfold -p2)"},

    // {"",0,0,O_SECTION,0,O_NODEFAULT,"","Type of locality"},

    // {"struct-local",0,0,O_ARG_BOOL,&struct_local,"false","bool","Structure local"},
    // {"sequ-local",0,0,O_ARG_BOOL,&sequ_local,"false","bool","Sequence local"},
    // {"free-endgaps",0,0,O_ARG_STRING,&free_endgaps,"----","spec","Whether and which end gaps are free. order: L1,R1,L2,R2"},

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

    {"min-prob",'p',0,O_ARG_DOUBLE,&min_prob,"0.0005","prob","Minimal probability"},
    {"max-diff-am",'D',0,O_ARG_INT,&max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&max_diff,"-1","diff","Maximal difference for alignment cuts"},
    //{"min-am-prob",'a',0,O_ARG_DOUBLE,&min_am_prob,"0.0005","amprob","Minimal Arc-match probability"},
    //{"min-bm-prob",'b',0,O_ARG_DOUBLE,&min_bm_prob,"0.0005","bmprob","Minimal Base-match probability"},

    // {"",0,0,O_SECTION,0,O_NODEFAULT,"","Special sauce options"},
    // {"kbest",'k',0,O_ARG_INT,&kbest_k,"-1","k","Find k-best alignments"},

    // {"",0,0,O_SECTION,0,O_NODEFAULT,"","Options for controlling MEA score (experimental, under construction)."},

    // {"mea-alignment",0,&opt_mea_alignment,O_NO_ARG,0,O_NODEFAULT,"","Do MEA alignment."},
    // {"probcons-file",0,&opt_probcons_file,O_ARG_STRING,&probcons_file,O_NODEFAULT,"file","Probcons parameter file."},

    // {"match-prob-method",0,0,O_ARG_INT,&match_prob_method,"0","int","Method for computation of match probs."},
    // {"temperature",0,0,O_ARG_INT,&temperature,"150","int","Temperature for PF-computation."},
    // {"pf-struct-weight",0,0,O_ARG_INT,&pf_struct_weight,"200","weight","Structure weight in PF-computation."},

    // {"mea-gapcost",0,&opt_mea_gapcost,O_NO_ARG,0,O_NODEFAULT,"","Use gap cost in mea alignment."},
    // {"mea-alpha",0,0,O_ARG_INT,&mea_alpha,"0","weight","Weight alpha for MEA."},
    // {"mea-beta",0,0,O_ARG_INT,&mea_beta,"200","weight","Weight beta for MEA."},
    // {"mea-gamma",0,0,O_ARG_INT,&mea_gamma,"100","weight","Weight gamma for MEA."},
    // {"probability-scale",0,0,O_ARG_INT,&probability_scale,"10000","scale","Scale for probabilities/resolution of mea score."},

    // {"write-match-probs",0,&opt_write_matchprobs,O_ARG_STRING,&matchprobs_file,O_NODEFAULT,"file","Write match probs to file (don't align!)."},
    // {"read-match-probs",0,&opt_read_matchprobs,O_ARG_STRING,&matchprobs_file,O_NODEFAULT,"file","Read match probabilities from file."},

    // {"write-arcmatch-scores",0,&opt_write_arcmatch_scores,O_ARG_STRING,&arcmatch_scores_file,O_NODEFAULT,"file","Write arcmatch scores (don't align!)."},
    // {"read-arcmatch-scores",0,&opt_read_arcmatch_scores,O_ARG_STRING,&arcmatch_scores_file,O_NODEFAULT,"file","Read arcmatch scores."},
    // {"read-arcmatch-probs",0,&opt_read_arcmatch_probs,O_ARG_STRING,&arcmatch_scores_file,O_NODEFAULT,"file","Read arcmatch probabilities (weight by mea_beta/100)."},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Constraints"},

    {"noLP",0,&opt_no_lonely_pairs,O_NO_ARG,0,O_NODEFAULT,"","No lonely pairs (only compatibility with locarna; no effect)"},
    {"anchorA",0,0,O_ARG_STRING,&seq_constraints_A,"","string","Anchor constraints sequence A."},
    {"anchorB",0,0,O_ARG_STRING,&seq_constraints_B,"","string","Anchor constraints sequence B."},
    {"ignore-constraints",0,&opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},

    //{"",0,0,O_SECTION_HIDE,0,O_NODEFAULT,"","Mode of operation"},
    //{"eval",0,&opt_eval,O_NO_ARG,0,O_NODEFAULT,"","Turn on evaluation mode."},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling Gecode"},
    {"c_d",0,0,O_ARG_INT,&c_d,"1","distance","Recomputation distance"},
    {"time-limit",0,&opt_time_limit,O_ARG_INT,&time_limit,O_NODEFAULT,"time","Search time limit (always searches for first solution)."},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Standard options"},

    {"help",'h',&opt_help,O_NO_ARG,0,O_NODEFAULT,"","This help"},
    {"version",'V',&opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","RNA sequences and pair probabilities"},

    {"",0,0,O_ARG_STRING,&file1,O_NODEFAULT,"file 1","Basepairs input file 1"},
    {"",0,0,O_ARG_STRING,&file2,O_NODEFAULT,"file 2","Basepairs input file 2"},
    {"",0,0,0,0,O_NODEFAULT,"",""}
};

// end option declaration
// ------------------------------------------------------------

/** \brief Main-function
 *  \relates RNAalignment
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

	cout << "Copyright Alessandro Dal Palu, Mathias Moehl, Sebastian Will, 2005-2009"<<endl<<endl;

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

    if (use_ribosum) {
	if (ribosum_file == "RIBOSUM85_60") {
	    if (opt_verbose) {
		std::cout <<"Use built-in ribosum."<<std::endl;
	    }
	    ribosum = new LocARNA::Ribosum85_60();
	} else {
	    ribosum = new LocARNA::RibosumFreq(ribosum_file);
	}
    }


    // ----------------------------------------
    // Scoring Parameter
    //
    LocARNA::ScoringParams 
	scoring_params(match_score,
		       mismatch_score,
		       indel_score,
		       indel_opening_score,
		       ribosum,
		       struct_weight,
		       tau_factor,
		       exclusion_score,
		       opt_exp_prob?exp_prob:-1,
		       0, // temperature
		       false, // opt_stacking not implemented in Carna
		       false,
		       0,0,0,0
		       );

    // ------------------------------------------------------------
    // Get input data and generate data objects
    //

    LocARNA::RnaData rnadataA(file1,false,false); // false->opt_stacking not implemented in Carna
    LocARNA::RnaData rnadataB(file2,false,false); // false->opt_stacking not implemented in Carna

    LocARNA::Sequence seqA=rnadataA.get_sequence();
    LocARNA::Sequence seqB=rnadataB.get_sequence();

    LocARNA::size_type lenA=seqA.length();
    LocARNA::size_type lenB=seqB.length();

    // --------------------
    // handle max_diff restriction

    LocARNA::TraceController trace_controller(seqA,seqB,NULL,max_diff,false);

    // ------------------------------------------------------------
    // Handle constraints (optionally)

    std::string seqCA = seq_constraints_A;
    std::string seqCB = seq_constraints_B;

    if (!opt_ignore_constraints) {
	if ( seqCA=="" ) seqCA = rnadataA.get_seq_constraints();
	if ( seqCB=="" ) seqCB = rnadataB.get_seq_constraints();
    }

    LocARNA::AnchorConstraints 
	seq_constraints(seqA.length(),seqCA,
			seqB.length(),seqCB);

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
    arc_matches = new LocARNA::ArcMatches(rnadataA,
					  rnadataB,
					  min_prob,
					  (max_diff_am!=-1)?(LocARNA::size_type)max_diff_am:std::max(lenA,lenB),
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

    LocARNA::Scoring scoring(seqA,seqB,arc_matches,0L,&scoring_params);

    // ------------------------------------------------------------
    // parameter for the alignment
    //
    LocARNA::AlignerParams aligner_params(0, //no_lonely_pairs does not make sense in carna
					  struct_local,
					  sequ_local,
					  free_endgaps,
					  trace_controller,
					  max_diff_am,
					  0, // min_am_prob and
					  0, // min_bm_prob are not used in Carna
					  false, // opt_stacking not implemented in Carna
					  seq_constraints
					  );


    // ------------------------------------------------------------
    // construct the constraint model / root space
    //
    RNAalignment* s = new RNAalignment(seqA,seqB,
				       *arc_matches,
				       aligner_params,scoring,
				       opt_gist	
				       );
    

    // ------------------------------------------------------------
    // run the search engine
    //
    if (opt_gist) {
#ifdef HAVE_GIST
	Gecode::Gist::Print<RNAalignment> p("Node explorer");
	Gecode::Gist::Options o;
	o.inspect.click(&p);

	o.c_d =  c_d;

	//Gecode::Gist::dfs(s,o);
	Gecode::Gist::bab(s,o);
#endif
    } else {
	Gecode::Search::Options o;

	o.c_d =  c_d;
	o.stop=0L;
	
	if (opt_time_limit) {
	    Gecode::Search::Stop *timestop=0L;
	    // activate next line for "if time limit is exceeded, no solution is produced"
	    // timestop = new Gecode::Search::TimeStop(time_limit);
	    // activate instead for "the first solution is always produced"
	    timestop = new Gecode::Search::TimeStop(std::numeric_limits<unsigned long int>::max());
	    o.stop=timestop;
	}
	
	// construct engine
	Gecode::BAB<RNAalignment> e(s,o);
	
	bool first_solution=true;

	// ----------------------------------------
	// enumerate solutions
	RNAalignment* ex=NULL;
	while ((ex = e.next()) && (ex != NULL)) {
	    // write each solution to stdout and optionally to files
	    // in clustal and pp format
	    
	    if (first_solution && opt_time_limit) {
		((Gecode::Search::TimeStop *)o.stop)->limit(time_limit);
		first_solution=false;
	    }
	    
	    ex->print(std::cout);

	    if (!clustal_out.empty()){
		ofstream outfile;
		outfile.open(clustal_out.c_str(),ios::out | ios::trunc);
		if (outfile.good()) {
		    ex->print_clustal_format(outfile,output_width);
		    outfile.close();
		} else {
		    std::cerr << "Cannot write solution to file "<<clustal_out<<std::endl;
		}
	    }

	    if (!pp_out.empty()){
		ofstream outfile;
		outfile.open(pp_out.c_str(),ios::out | ios::trunc);
		if (outfile.good()) {
		    ex->print_pp_format(outfile,
					bpsA,
					bpsB,
					scoring, 
					seq_constraints, 
					output_width,
					opt_alifold_consensus_dp
					);
		    outfile.close();
		} else {
		   std::cerr << "Cannot write solution to file "<<pp_out<<std::endl;
		}
	    }

	    if (ex) delete ex;
	}
	

	Gecode::Search::Statistics stats = e.statistics();

	cout << "Nodes:  "<<stats.node << std::endl
	     << "Fail:   "<<stats.fail << std::endl
	     << "Depth:  "<<stats.depth << std::endl
	     << "Memory: "<<stats.memory << std::endl;

	if (e.stopped()) {
	    cout << "Time limit exceeded."<<std::endl;
	}

	if (o.stop) delete o.stop;

    }

    if (ribosum) delete ribosum;
    if (arc_matches) delete arc_matches;
    if (s) delete s;
    return 0;
}
