#include <vector>
#include <cstdint>
#include "organo_opts.hpp"
#include "andistmatrix.hpp"
#include "fastcluster.h"
#include "utils.hpp"

struct cluster{
	int label;
	double dist;
	int n;
	cluster(){};
	cluster(int l, double d, int _n):label(l), dist(d), n(_n){};
};

void merge_singleton_clusters(const organo_opts& params, std::vector<uint32_t>& spanning, double*& distmat, int*& labels, std::vector<std::pair<int,int>>& cluster_labels_counts)
{
	std::vector<cluster> cluster_dists(cluster_labels_counts.size());
	uint32_t i, j;

	for(i = 0; i < cluster_labels_counts.size(); ++i){
		int subj_label = cluster_labels_counts[i].first;
		double subj_distsum = 0.0;
		int subj_n = 0;
		if(labels[i] == subj_label){
			for(j = i + 1; j < spanning.size(); ++j){
				if(labels[j] == subj_label){
					subj_distsum += distmat[(static_cast<std::ptrdiff_t>(2*spanning.size()-3-(i))*(i)>>1)+(j)-1];
					++subj_n;
				}
			}
		}
		cluster_dists[i] = cluster(subj_label, subj_distsum, subj_n);
	}


	uint32_t singleton_start_i = 0;
	while(singleton_start_i < cluster_labels_counts.size() && cluster_labels_counts[singleton_start_i].second >= params.mincov) ++singleton_start_i;

	uint32_t x, y;
	for(x = singleton_start_i; x < cluster_dists.size(); ++x){
		int min_target_label = -1;
		double min_target_distsum;
		int min_target_n;

		for(y = 0; y < singleton_start_i; ++y){
			cluster acc = cluster_dists[x];
			int target_label = cluster_dists[y].label;

			for(i = 0; i < spanning.size(); ++i){
				if(labels[i] == acc.label){
					for(j = 0; j < spanning.size(); ++j){
						if(labels[j] == target_label){
							int a = i < j ? i : j;
							int b = i > j ? i : j;
							acc.dist += distmat[(static_cast<std::ptrdiff_t>(2*spanning.size()-3-(a))*(a)>>1)+(b)-1];
							++acc.n;
						}
					}
				}
			}

			if(min_target_label < 0 || (acc.dist/acc.n) < (min_target_distsum/min_target_n)) {
				min_target_label = target_label;
				min_target_distsum = acc.dist;
				min_target_n = acc.n;
			}
		}
		for(i = 0; i < spanning.size(); ++i) if(labels[i] == cluster_dists[x].label) labels[i] = min_target_label;	
	}
}

void cluster_task(const organo_opts& params, andistmatrix& distmatrix, std::vector<uint32_t>& spanning, std::vector<std::vector<uint32_t>>& ccs)
{
	if(spanning.empty()) return;
	else{
		if(params.maxalleles == 1) ccs.emplace_back(spanning);
		else if(spanning.size() == 1){
			std::vector<uint32_t> cc = {spanning.front()};
			ccs.emplace_back(cc);
		}
		else if(spanning.size() == 2){
			std::vector<uint32_t> cc;
			if(distmatrix.get_dist(spanning[0], spanning[1]) <= params.maxerror){
				cc.emplace_back(spanning.front());
				cc.emplace_back(spanning.back());
				ccs.emplace_back(cc);
			}
			else{
				cc.emplace_back(spanning.front());
				ccs.emplace_back(cc);
				cc[0] = spanning.back();
				ccs.emplace_back(cc);
			}
		}
		else {
			uint32_t total_d = (spanning.size() * (spanning.size() - 1)) / 2;
			double* distmat = new double[total_d];
			uint32_t k,i,j;
			for(i=k=0; i < spanning.size(); ++i) {
				for(j=i+1; j < spanning.size(); ++j) {
			    	distmat[k] = distmatrix.get_dist(spanning[i], spanning[j]);
			    	++k;
			  	}
			}
			int* labels = new int[spanning.size()];
			int* merge = new int[2*(spanning.size()-1)];
			double* height = new double[spanning.size()-1];
			hclust_fast(spanning.size(), distmat, HCLUST_METHOD_AVERAGE, merge, height);
			for(i = spanning.size() - 2; i > 0; --i) if(height[i] <= params.maxerror) break;
			int clusters = spanning.size() - i - 1;
			//output all alleles
			if(params.maxalleles <= 0 || clusters < params.maxalleles) cutree_cdist(spanning.size(), merge, height, (double)params.maxerror, labels);
			//guide cluster merging to meet specified max allele number
			else{
				std::vector<int> cluster_labels;
				cutree_cdist(spanning.size(), merge, height, (double)params.maxerror, labels);
				for(i = 0; i < spanning.size(); ++i) cluster_labels.emplace_back(labels[i]);
				std::vector<std::pair<int,int>> cluster_labels_counts;
				group_count_task(cluster_labels, cluster_labels_counts);
				int clusters_min_n = 0;
				for(i = 0; i < cluster_labels_counts.size(); ++i) if(cluster_labels_counts[i].second >= params.mincov) ++clusters_min_n;
				if(clusters_min_n > 0 && clusters_min_n <= params.maxalleles) merge_singleton_clusters(params, spanning, distmat, labels, cluster_labels_counts);
				else{
					clusters = params.maxalleles + 1;
					cluster_labels.clear();
					cutree_k(spanning.size(), merge, params.maxalleles + 1, labels);
					for(i = 0; i < spanning.size(); ++i) cluster_labels.emplace_back(labels[i]);
					cluster_labels_counts.clear();
					group_count_task(cluster_labels, cluster_labels_counts);
					if(cluster_labels_counts.back().second >= params.mincov) cutree_k(spanning.size(), merge, params.maxalleles, labels);
					else merge_singleton_clusters(params, spanning, distmat,labels,cluster_labels_counts);
				}

				/**
				std::vector<int> cluster_labels;
				if(clusters == params.maxalleles){
					cutree_k(spanning.size(), merge, params.maxalleles, labels);
					for(i = 0; i < spanning.size(); ++i) cluster_labels.emplace_back(labels[i]);
					std::vector<std::pair<int,int>> cluster_labels_counts;
					group_count_task(cluster_labels, cluster_labels_counts);
					//check if second allele is not supported enough
					if(cluster_labels_counts.back().second < params.mincov){
						for(i = 0; i < spanning.size(); ++i) if(labels[i] == cluster_labels_counts.back().first) labels[i] = cluster_labels_counts.front().first;
					}
				}
				else{
					clusters = params.maxalleles + 1;
					cutree_k(spanning.size(), merge, params.maxalleles + 1, labels);
					
					for(i = 0; i < spanning.size(); ++i) cluster_labels.emplace_back(labels[i]);
					std::vector<std::pair<int,int>> cluster_labels_counts;
					group_count_task(cluster_labels, cluster_labels_counts);
					
					if(cluster_labels_counts.back().second >= params.mincov) cutree_k(spanning.size(), merge, params.maxalleles, labels);
					else merge_singleton_clusters(params, spanning, distmat,labels,cluster_labels_counts);
				}
				*/

			}

			std::vector<int> final_labels(spanning.size());
			for(i = 0; i < spanning.size(); ++i) final_labels[i] = labels[i];
			std::vector<std::pair<int,int>> final_label_counts;
			group_count_task(final_labels, final_label_counts);

			ccs.resize(final_label_counts.size());
			for(i = 0; i < final_label_counts.size(); ++i) {
				int label = (int)final_label_counts[i].first;
				for(j = 0; j < spanning.size(); ++j){
					if(labels[j] == label) ccs[i].emplace_back(spanning[j]);
				}
			}

			delete[] distmat;
			delete[] merge;
			delete[] height;
			delete[] labels;
		}
	}
}