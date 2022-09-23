#include <string>
#include <unordered_map>

void load_ha_tagged_reads(
	const std::string& file, 
	std::unordered_map<std::string, int>& map
	);