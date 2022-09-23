#include <vector>
#include <exception>
#include <iostream>
#include <cstddef>
#include "andistmatrix.hpp"

andistmatrix::andistmatrix(uint32_t _n): n(_n), d_n((n * (n - 1)) / 2)
{
	distmat.resize(d_n, 1.0);
}

void andistmatrix::set_dist(uint32_t i, uint32_t j, double d)
{
	if(i == j) throw std::exception();
	int a = i < j ? i : j;
	int b = i > j ? i : j;
	distmat[(static_cast<std::ptrdiff_t>(2*n-3-(a))*(a)>>1)+(b)-1] = d;
}

double andistmatrix::get_dist(uint32_t i, uint32_t j) const 
{
	if(i == j) throw std::exception();
	int a = i < j ? i : j;
	int b = i > j ? i : j;
	return distmat[(static_cast<std::ptrdiff_t>(2*n-3-(a))*(a)>>1)+(b)-1];
}

uint32_t andistmatrix::representative(std::vector<uint32_t>& indeces) const
{
	if(indeces.size() < 2) return indeces[0];
	else{
		//same order as indeces
		std::vector<std::pair<double, uint32_t>> sum_dists(indeces.size(), std::make_pair(0.0,0));
		//iterate through each index
		for(uint32_t i = 0; i < indeces.size(); ++i){
			for(uint32_t j = 0; j < indeces.size(); ++j){
				if(i != j) {
					sum_dists[i].first += get_dist(indeces[i], indeces[j]);
					++sum_dists[i].second;
				}
			}
		}

		
		uint32_t smallest_index = 0;
		double smallest_dist = sum_dists[0].first / sum_dists[0].second;
		for(uint32_t i = 1; i < indeces.size(); ++i){
			double dist = sum_dists[i].first / sum_dists[i].second;
			if(dist < smallest_dist) {
				smallest_index = i;
				smallest_dist = dist;
			}
		}
		return indeces[smallest_index];
	}
}

uint32_t andistmatrix::size()const{return distmat.size();}