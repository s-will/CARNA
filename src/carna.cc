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
//
// Options
//
#include <LocARNA/options.hh>

//! \brief Structure for command line parameters of locarna
//!
//! Encapsulating all command line parameters in a common structure
//! avoids name conflicts and makes downstream code more informative.
//!
struct command_line_parameters {

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

    // int kbest_k;

    bool opt_ignore_constraints;

    int pf_struct_weight;

    int c_d;

    int time_limit;

    bool opt_lower_bound;
    int lower_score_bound;
    bool opt_upper_bound;
    int upper_score_bound;

    bool relaxed_anchors; //!< strict or relaxed anchor constraints

    // bool opt_mea_gapcost;
    // int mea_alpha;
    // int mea_beta;
    // int mea_gamma;
    // int probability_scale;

    // bool opt_eval;
};

//! \brief holds command line parameters  
command_line_parameters clp;


LocARNA::option_def my_options[] = {
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Scoring parameters"},

    {"match",'m',0,O_ARG_INT,&clp.match_score,"50","score","Match score"},
    {"mismatch",'M',0,O_ARG_INT,&clp.mismatch_score,"0","score","Mismatch score"},
    {"ribosum-file",0,0,O_ARG_STRING,&clp.ribosum_file,"RIBOSUM85_60","f","Ribosum file"},
    {"use-ribosum",0,0,O_ARG_BOOL,&clp.use_ribosum,"true","bool","Use ribosum scores"},
    {"indel",'i',0,O_ARG_INT,&clp.indel_score,"-350","score","Indel score"},
    {"indel-opening",0,0,O_ARG_INT,&clp.indel_opening_score,"-500","score","Indel opening score"},
    {"struct-weight",'s',0,O_ARG_INT,&clp.struct_weight,"200","score","Maximal weight of 1/2 arc match"},
    {"exp-prob",'e',&clp.opt_exp_prob,O_ARG_DOUBLE,&clp.exp_prob,O_NODEFAULT,"prob","Expected probability"},
    {"tau",'t',0,O_ARG_INT,&clp.tau_factor,"0","factor","Tau factor in percent"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling output"},
    
    {"gist",0,&clp.opt_gist,O_NO_ARG,0,O_NODEFAULT,"",
#ifdef HAVE_GIST
     "Use gist for interactive/graphical search."
#else
     "Use gist for graphical search (feature disabled, recompile to activate)."
#endif
    },

    {"width",'w',0,O_ARG_INT,&clp.output_width,"120","columns","Output width"},
    {"clustal",0,&clp.opt_clustal_out,O_ARG_STRING,&clp.clustal_out,O_NODEFAULT,"file","Clustal output"},
    {"pp",0,&clp.opt_pp_out,O_ARG_STRING,&clp.pp_out,O_NODEFAULT,"file","PP output"},
    // {"local-output",'L',&clp.opt_local_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment"},
    // {"pos-output",'P',&clp.opt_pos_output,O_NO_ARG,0,O_NODEFAULT,"","Output only local sub-alignment positions"},
    
#ifdef HAVE_LIBRNA
    {"alifold-consensus-dp",0,&clp.opt_alifold_consensus_dp,O_NO_ARG,0,O_NODEFAULT,"","Compute consensus dot plot by alifold"},
#endif
    
    {"write-structure",0,&clp.opt_write_structure,O_NO_ARG,0,O_NODEFAULT,"","Write guidance structure in output"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Heuristics for speed accuracy trade off"},

    {"min-prob",'p',0,O_ARG_DOUBLE,&clp.min_prob,"0.01","prob","Minimal probability"},
    {"max-bps-length-ratio",0,0,O_ARG_DOUBLE,&clp.max_bps_length_ratio,"0.0","factor","Maximal ratio of #base pairs divided by sequence length (default: no effect)"},
    {"max-diff-am",'D',0,O_ARG_INT,&clp.max_diff_am,"-1","diff","Maximal difference for sizes of matched arcs"},
    {"max-diff",'d',0,O_ARG_INT,&clp.max_diff,"-1","diff","Maximal difference for alignment cuts"},
    {"max-diff-at-am",0,0,O_ARG_INT,&clp.max_diff_at_am,"-1","diff","Maximal difference for alignment traces, only at arc match positions"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Constraints"},

    {"noLP",0,&clp.opt_no_lonely_pairs,O_NO_ARG,0,O_NODEFAULT,"","No lonely pairs (only used when predicing ensemble porobabilities and for compatibility with locarna; otherwise no effect)"},
    {"maxBPspan",0,0,O_ARG_INT,&max_bp_span,"-1","span","Limit maximum base pair span (default=off)"},
    {"relaxed-anchors",0,&clp.relaxed_anchors,O_NO_ARG,0,O_NODEFAULT,"","Relax anchor constraints (default=off)"},
    {"ignore-constraints",0,&clp.opt_ignore_constraints,O_NO_ARG,0,O_NODEFAULT,"","Ignore constraints in pp-file"},
    
    {"lb",0,&clp.opt_lower_bound,O_ARG_INT,&clp.lower_score_bound,O_NODEFAULT,"score","Lower score bound"},
    {"ub",0,&clp.opt_upper_bound,O_ARG_INT,&clp.upper_score_bound,O_NODEFAULT,"score","Upper score bound"},
    

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Controlling Gecode"},
    {"c_d",0,0,O_ARG_INT,&clp.c_d,"1","distance","Recomputation distance"},
    {"time-limit",0,0,O_ARG_INT,&clp.time_limit,"300000","time","Time limit in ms (always search first solution; turn off by 0)."},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Standard options"},

    {"help",'h',&clp.opt_help,O_NO_ARG,0,O_NODEFAULT,"","This help"},
    {"version",'V',&clp.opt_version,O_NO_ARG,0,O_NODEFAULT,"","Version info"},
    {"verbose",'v',&clp.opt_verbose,O_NO_ARG,0,O_NODEFAULT,"","Verbose"},

    {"",0,0,O_SECTION,0,O_NODEFAULT,"","RNA sequences and pair probabilities"},


    {"",0,0,O_SECTION_HIDE,0,O_NODEFAULT,"","Hidden Options"},
    {"ribofit",0,0,O_ARG_BOOL,&clp.opt_ribofit,"false","bool","Use Ribofit base and arc match scores (overrides ribosum)"},
    
    {"",0,0,O_SECTION,0,O_NODEFAULT,"","Input_files RNA sequences and pair probabilities"},

    {"",0,0,O_ARG_STRING,&clp.fileA,O_NODEFAULT,"file 1","Input file 1  ;; input files can be fasta, dotplot, pp, or aln files"},
    {"",0,0,O_ARG_STRING,&clp.fileB,O_NODEFAULT,"file 2","Input file 2  ;; and can speficy structure and anchor constraints."},
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

    if (clp.opt_help) {
	cout << VERSION_STRING<<endl;

	cout << "A tool for pairwise Alignment of RNA."<<endl<<endl;

	print_help(argv[0],my_options);

	// cout << "Report bugs to <will (at) informatik.uni-freiburg.de>."<<endl<<endl;
	exit(0);
    }

    if (clp.opt_version || clp.opt_verbose) {
	cout << VERSION_STRING<<endl;
	if (clp.opt_version) exit(0); else cout <<endl;
    }


    if (!process_success) {
	std::cerr << "ERROR --- "
		  <<LocARNA::O_error_msg<<std::endl;
	printf("USAGE: ");
	print_usage(argv[0],my_options);
	printf("\n");
	exit(-1);
    }

    if (clp.opt_verbose) {
	print_options(my_options);
    }

    //
    // END process options
    // ----------------------------------------


#ifndef HAVE_GIST
    if (clp.opt_gist) {
	std::cerr << "Graphical output (option --gist) was disabled when compiling." << std::endl
		  << "For enabling this feature, please reconfigure and recompile." << std::endl;
	exit(-1);
    }
#endif


    // options consistency
    if (clp.struct_weight<0) {
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
    
    if (clp.opt_ribofit) {
	ribofit = new LocARNA::Ribofit_will2014;
    }
    
    if (clp.use_ribosum) {
	if (clp.ribosum_file == "RIBOSUM85_60") {
	    if (clp.opt_verbose) {
                std::cout <<"Use built-in ribosum."<<std::endl;
	    }
	    ribosum = new LocARNA::Ribosum85_60;
	} else {
	    ribosum = new LocARNA::RibosumFreq(clp.ribosum_file);
	}
    }
    
    // ------------------------------------------------------------
    // Get input data and generate data objects
    //
    LocARNA::PFoldParams pfparams(clp.opt_no_lonely_pairs, false, max_bp_span, 2); //no_lonely_pairs,clp.opt_stacking, maxBPspan, dangling
    
    LocARNA::RnaData *rna_dataA=0;
    try {
	rna_dataA = new LocARNA::RnaData(clp.fileA,clp.min_prob,clp.max_bps_length_ratio,pfparams); 
    } catch (LocARNA::failure &f) {
	std::cerr << "ERROR:\tfailed to read from file "<<clp.fileA <<std::endl
		  << "\t"<< f.what() <<std::endl;
	return -1;
    }
    
    LocARNA::RnaData *rna_dataB=0;
    try {
	rna_dataB = new LocARNA::RnaData(clp.fileB,clp.min_prob,clp.max_bps_length_ratio,pfparams);
    } catch (LocARNA::failure &f) {
	std::cerr << "ERROR: failed to read from file "<<clp.fileB <<std::endl
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

    LocARNA::TraceController trace_controller(seqA,seqB,NULL,clp.max_diff,false);

    // ------------------------------------------------------------
    // Handle constraints (optionally)
    
    LocARNA::AnchorConstraints 
	seq_constraints(lenA,
			seqA.annotation(LocARNA::MultipleAlignment::AnnoType::anchors).single_string(),
			lenB,
			seqB.annotation(LocARNA::MultipleAlignment::AnnoType::anchors).single_string(),
                        !clp.relaxed_anchors
                        );
    
    if (clp.opt_verbose) {
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
					  clp.min_prob,
					  (clp.max_diff_am!=-1)
					  ?(LocARNA::size_type)clp.max_diff_am
					  :std::max(lenA,lenB),
					  clp.max_diff_at_am!=-1
					  ? (LocARNA::size_type)clp.max_diff_at_am
					  : std::max(lenA,lenB),
					  trace_controller,
					  seq_constraints
					  );

    LocARNA::BasePairs bpsA = arc_matches->get_base_pairsA();
    LocARNA::BasePairs bpsB = arc_matches->get_base_pairsB();

    // ----------------------------------------
    // report on input in verbose mode
    if (clp.opt_verbose) {
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

    double my_exp_probA = clp.opt_exp_prob?clp.exp_prob:LocARNA::prob_exp_f(lenA);
    double my_exp_probB = clp.opt_exp_prob?clp.exp_prob:LocARNA::prob_exp_f(lenB);


    // ----------------------------------------
    // Scoring Parameter
    //
    LocARNA::ScoringParams 
	scoring_params(clp.match_score,
		       clp.mismatch_score,
		       clp.indel_score,
		       clp.indel_score, // = indel loop score (not used by CARNA)
		       clp.indel_opening_score,
		       clp.indel_opening_score, // = indel loop opening score (not used by CARNA)
                       ribosum,
                       ribofit,
		       0, // unpaired penalty (not configurable in CARNA)
		       clp.struct_weight,
		       clp.tau_factor,
		       clp.exclusion_score,
		       my_exp_probA,
		       my_exp_probB,
		       0, // temperature
		       false, // clp.opt_stacking not implemented in Carna
		       false, // clp.opt_new_stacking not implemented in Carna
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
	 . struct_local(clp.struct_local)
	 . sequ_local(clp.sequ_local)
	 . free_endgaps(clp.free_endgaps)
	 . max_diff_am(clp.max_diff_am)
	 . max_diff_at_am(clp.max_diff_at_am)
	 . trace_controller(trace_controller)
	 . stacking(false)
	 . constraints(seq_constraints)
	 )
	. lower_score_bound(clp.opt_lower_bound?clp.lower_score_bound:Gecode::Int::Limits::min)
	. upper_score_bound(clp.opt_upper_bound?clp.upper_score_bound:Gecode::Int::Limits::max)
	. gist(clp.opt_gist)
	. output_width(clp.output_width)
	;
        
    RnaAlignment *rna_alignment_space = 
	new RnaAlignment(rna_alignment_params);

    // ------------------------------------------------------------
    // run the search engine
    //
    if (clp.opt_gist) {
#ifdef HAVE_GIST
	Gecode::Gist::Print<RnaAlignment> p("Node explorer");
	Gecode::Gist::Options o;
	o.inspect.click(&p);

	o.c_d =  clp.c_d;

	//Gecode::Gist::dfs(s,o);
	Gecode::Gist::bab(rna_alignment_space,o);
#endif
    } else {
	Gecode::Search::Options o;

	o.c_d =  clp.c_d;
	o.stop=0L;
	
	if (clp.time_limit>0) {
	    Gecode::Search::Stop *timestop=0L;
	    // activate next line for "if time limit is exceeded, no solution is produced"
	    // timestop = new Gecode::Search::TimeStop(clp.time_limit);
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
	    
	    if (first_solution && clp.time_limit>0) {
		((Gecode::Search::TimeStop *)o.stop)->limit(clp.time_limit);
		first_solution=false;
	    }
	    
	    ex->print(std::cout);

	    if (!clp.clustal_out.empty()){
		ofstream outfile;
		outfile.open(clp.clustal_out.c_str(),ios::out | ios::trunc);
		if (outfile.good()) {
		    
		    assert(ex->all_assigned());
		    LocARNA::Alignment alignment=ex->to_alignment();
			
		    LocARNA::MultipleAlignment ma(alignment,clp.opt_local_file_output);
		    
		    outfile << "CLUSTAL W --- "<<PACKAGE_STRING;
		    
		    // for legacy, clustal files of pairwise alignments contain the score 
		    if (seqA.num_of_rows()==1 && seqB.num_of_rows()==1)
			outfile  <<" --- Score: " << ex->score();
		    outfile <<std::endl<<std::endl;
		    
		    ma.write(outfile,clp.output_width);
		    outfile.close();
		} else {
		    std::cerr << "Cannot write solution to file "<<clp.clustal_out<<std::endl;
		}
	    }

	    if (!clp.pp_out.empty()){
		ofstream outfile;
		outfile.open(clp.pp_out.c_str(),ios::out | ios::trunc);
		
		if (outfile.good()) {
		    assert(ex->all_assigned());
		    
		    LocARNA::Alignment alignment=ex->to_alignment();

		    if (clp.opt_alifold_consensus_dp) {
			LocARNA::PFoldParams pfparams(clp.opt_no_lonely_pairs,false,max_bp_span,2);
			
			LocARNA::MultipleAlignment ma(alignment,clp.opt_local_file_output);
			LocARNA::RnaEnsemble ens(ma,pfparams,false,true); // alifold the alignment
			// construct consensus rna data from ensemble
			LocARNA::RnaData consensus(ens,clp.min_prob,clp.max_bps_length_ratio,pfparams);
			consensus.write_pp(outfile); // write alifold dot plot
		    } else {
			// compute averaged consensus base pair probabilities
			LocARNA::RnaData consensus(*rna_dataA,
						   *rna_dataB,
						   alignment,
						   my_exp_probA,
						   my_exp_probB,
						   clp.opt_local_file_output);
			consensus.write_pp(outfile); // write averaged dot plot
		    }
		    
		    outfile.close();
		} else {
		    std::cerr << "Cannot write solution to file "<<clp.pp_out<<std::endl;
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
