#include <string>
#include <vector>
#include <algorithm>
#include "utils.hpp"


void homopolymer_compressor(const char* ptr, std::string& str)
{
    char last = *ptr;
    uint32_t count = 1;
    ++ptr;
    while(*ptr != '\0'){
        if(*ptr == last) {
            ++count;
            ++ptr;
        }
        else {
            str.push_back(last);
            count = 1;
            last = *ptr;
        }
    }
    if(last != '\0') str.push_back(last);
}

void group_count_task(std::vector<int>& elements, std::vector<std::pair<int,int>>& counts)
{
	if(!elements.empty()){
		sort(elements.begin(), elements.end());
		auto it = elements.begin();
		auto value = *it;
		int count = 1;
		++it;
		while(it != elements.end()){
			if(*it == value) ++count;
			else {
				counts.emplace_back(std::make_pair(value, count));
				value = *it;
				count = 1;
			}
			++it;
		}
		counts.emplace_back(std::make_pair(value, count));
		sort(counts.begin(), counts.end(), [](const std::pair<int, int>& i, const std::pair<int, int>& j) -> bool{ return i.second > j.second; });
	}
}