#ifndef ANBAMTOOLS_HPP
#define ANBAMTOOLS_HPP

#include <htslib/sam.h>
#include <vector>
#include <memory>
#include <string>
#include "abg.hpp"

const char init_header[] = "@HD\tVN:1.4\tSO:unknown\n";
const uint8_t base_qual = (uint8_t)0;

/**
 * Self-contained struct for (independently) loading a BAM file (i.e. multi-threading). It contains an
 * 'init' and 'destroy' methods for cleaning.
 * 
 * @param samFile     SAM/BAM-file file pointer
 * @param bam_hdr_t   SAM/BAM header
 * @param bool        Whether SAM/BAM is indexed
 * @param hts_idx_t   SAM/BAM header pointer
 * @param bam1_t      SAM/BAM read-alignment pointer
 */
struct ansbam {
	samFile *fp;
	bam_hdr_t *header;
	bool index;
	hts_idx_t *idx;
	bam1_t* read;

	ansbam();
	void init(const std::string, bool);
	void destroy();
};


/**
 * Struct that accompanies the 'abg_generate' function (see below), for determining both whether a 
 * an AB-graph generation was successful, and the spanning-type.
 * 
 * @param bool    Whether an AB-graph was successfully generated
 * @param bool    Whether the corresponding read spans left-breakpoint
 * @param bool    Whether the corresponding read spans right-breakpoint
 */
struct abg_generate_msg {
	bool successful;
	bool spanning_l;
	bool spanning_r;
	abg_generate_msg():successful(true),spanning_l(true), spanning_r(true){};
	bool spanning() const{return spanning_l && spanning_r;};
};

/**
 * Function to project reference-coordinates on a given read. The projection is stored in the given 
 * vector such that each nucleotide in the read has a corresponding reference coordinate. A value of
 * '-1' indicates that the nucleotide does not exist in the reference. Additionally, the given pair 
 * of bool variables are changed to reflect if the read-alignment is clipped.
 * 
 * @param bam1_t       Read-alignment
 * @param bool         Read is left-clipped
 * @param bool         Read is right-clipped
 * @param vector<int>  Projection of reference-coordinatees on the read  
 */
void project_positions(
	bam1_t* alignment, 
	bool& clipped_l,
	bool& clipped_r,
	std::vector<int>& refcoords
	);

/**
 * Function to find corresponding start/end coordinates (breakpoints) of a given read given clipping
 * information and the projected reference coordinates. Breakpoints as stored in the given pointer 
 * of int-pairs. Failure to find the breakpoints as well as the spanning-type are reflected in the
 * given 'abg_generate_msg' variable.
 * 
 * @param int                         Start coordinate of breakpoint
 * @param int                         End coordinate of breakpoint
 * @param bool                        Read-alignment is left-clipped
 * @param bool                        Read-alignment is right-clipped
 * @param abg_generate_msg            Status message
 * @param vector<int>                 Projected positions
 * @param unique_ptr<pair<int,int>>   Varianble to store corresponding breakpoints
 * 
 */
void get_breakpoints(
	const int& start, 
	const int& end,
	const bool& clipped_l,
	const bool& clipped_r,
	abg_generate_msg& msg,
	std::vector<int>& projectedcoords,
	std::unique_ptr<std::pair<int, int>>& ptr
	);

/**
 * Function to generage an AB-graph given the start/end breakpoint coordinates and read-alignment.
 * The AB-graph is inferred from the status message and corresponding sequence of the edge.
 * 
 * @param bam1_t            Read-alignment
 * @param int               Start of breakpoint coordinate
 * @param int               End of breakpoint coordinate
 * @param abg_generate_msg  Status message
 * @param seq               Variable to store corresponding sequence of the edge
 */
void abg_generate(
	bam1_t* read,
	const int rstart,
	const int rend,
	abg_generate_msg& msg,
	std::string& seq
	);

/**
 * 
 */
void write2bam(
    samFile* outfile, 
    const std::string& region_bed, 
    const std::vector<uint32_t>& indeces, 
    const std::vector<abg>& abg_reads, 
    const std::vector<std::string>& edge_seqs
    );


class abg_block_iter {
    public:
        ansbam bam_inst;
        std::vector<abg_block> loaded_blocks;

        abg_block_iter(const std::string&, const uint32_t&, const uint32_t&, const bool&);
        bool next();
        void destroy();

    private:
        abg_block current_block;
        uint32_t block_size;
        uint32_t maxcov;
        bool hp;
};

#endif