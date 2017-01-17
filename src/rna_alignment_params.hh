#ifndef RNA_ALIGNMENT_PARAMS
#define RNA_ALIGNMENT_PARAMS

#include <gecode/int.hh>
#include <LocARNA/params.hh>

/**
 * @brief parameters for RnaAlignment
 */
class RnaAlignmentParams: public LocARNA::AlignerParams {  
    friend class RnaAlignment;
    
protected:
    // specific parameters
    int lower_score_bound_;
    int upper_score_bound_;
    bool gist_;
    size_t output_width_;
    
public:
    /** 
     * Construct with default parameters
     */
    RnaAlignmentParams(const AlignerParams &params)
	: AlignerParams(params),
	  // specific parameters defaults:
	  lower_score_bound_(Gecode::Int::Limits::min),
	  upper_score_bound_(Gecode::Int::Limits::max),
	  gist_(false)
    {}
    
    
    // specific parameters set methods:
    
    /**
     * @brief set lower bound
     * @param lower_bound the lower bound
     */
    RnaAlignmentParams &
    lower_score_bound(int lower_score_bound) {lower_score_bound_=lower_score_bound; return *this;}

    /**
     * @brief set upper bound
     * @param lower_bound the upper bound
     */
    RnaAlignmentParams &
    upper_score_bound(int upper_score_bound) {upper_score_bound_=upper_score_bound; return *this;}

    /**
     * @brief set gist flag
     * @param gist whether to show search tree via gist
     */
    RnaAlignmentParams &
    gist(bool gist) {gist_=gist; return *this;}

    /**
     * @brief set output width
     * @param width the output width
     */
    RnaAlignmentParams &
    output_width(size_t output_width) {output_width_=output_width; return *this;}
    
};


#endif // RNA_ALIGNMENT_PARAMS

