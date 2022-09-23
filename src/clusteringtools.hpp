#include <vector>
#include "organo_opts.hpp"
#include "andistmatrix.hpp"

void cluster_task(
	const organo_opts& params, 
	andistmatrix& distmatrix, 
	std::vector<uint32_t>& spanning, 
	std::vector<std::vector<uint32_t>>& ccs
	);