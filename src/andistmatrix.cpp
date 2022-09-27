#include <vector>
#include <exception>
#include <iostream>
#include <cstddef>
#include <cmath>
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

double andistmatrix::binned_kde(const int& max_bandwidth, const int& min_support, const double& error, std::vector<uint32_t>& indeces) const
{

	std::vector<std::vector<double>> bin_dists;
	for(uint32_t i = 0; i < 100; ++i) bin_dists.emplace_back(std::vector<double>());

	for(uint32_t i = 0; i < indeces.size(); ++i){
		for(uint32_t j = i + 1; j < indeces.size(); ++j){
			auto d = get_dist(indeces[i], indeces[j]);
			if(d == 1.0) d = 0.99;
			bin_dists[std::floor(d*100)].emplace_back(d);
		}
	}

	int error_i = 0;
	uint32_t start_i = 0;
	if(bin_dists[0].size() >= min_support) start_i = error_i;
	else {
		while(start_i < 100) if(bin_dists[start_i].size() < min_support) ++start_i; else break;
	}
	//std::cerr << "start at " << start_i << '\n';
	int current_max = 0;
	int bandwidth_iter = 0;
	while(bandwidth_iter < max_bandwidth && start_i < bin_dists.size()){
		int acc_max = 0;
		for(uint32_t i = 0; i < max_bandwidth; ++i) acc_max += bin_dists[start_i + i].size();
		//std::cerr << "index: " << start_i << ',' << acc_max << '\n';
		if(acc_max == 0 || acc_max < current_max) ++bandwidth_iter;
		else {
			bandwidth_iter = 0;
			current_max = acc_max;
		}
		++start_i;
	}

	//std::cerr << "max at: " << (start_i - max_bandwidth - 1) << '\n';
	double mean_sum = 0.0;
	uint32_t mean_n = 0;

	uint32_t i;
	if(start_i < max_bandwidth - 1) i = 0; else i = start_i - max_bandwidth - 1;

	for(; i < start_i; ++i){
		mean_n += bin_dists[i].size();
		for(const auto& d : bin_dists[i]) mean_sum += d;
	}

	if(mean_n == 0) return 0.0; else return (mean_sum / mean_n) + 0.005;

}

uint32_t andistmatrix::size()const{return distmat.size();}