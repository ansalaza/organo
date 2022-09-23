#include <string>
#include "organo_opts.hpp"


/**
 * Function to generate AB-graphs for all reads in a given BAM-file that span/contained target 
 * regions in the BED-file. Optionally, provide a two-column TSV mapping a read-name with its 
 * corresponding haplotype (int) ID.
 * 
 * The resulting AB-graphs are compressed in individual BAM-alignment entries and redirected to
 * stdout.
 * 
 * @params string
 * @params string
 * @params string
 * @params organo_opts
 * 
 * @return void
 */
void generate(
	const std::string bam, 
	const std::string bed, 
	const std::string haplotypes, 
	const organo_opts& params
	);