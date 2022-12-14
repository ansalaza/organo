#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "organo_opts.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "pairwise_alignments.hpp"
#include "abg.hpp"
#include "andistmatrix.hpp"
#include "utils.hpp"

void spanning_alignment(wfa::WFAligner& aligner, std::string& query, abg& q_tag, std::string& ref)
{
	//query fully spans
	if(q_tag.spanning()) aligner.alignEnd2End(query, ref);
	//query spans left
	else if(q_tag.spanning_l) aligner.alignEndsFree(query, 0, 0, ref, 0, ref.size());
	//query spans right
	else if(q_tag.spanning_r) aligner.alignEndsFree(query, 0, 0, ref, ref.size(), 0);
	//query is contained
	else aligner.alignEndsFree(query, 0, 0, ref, ref.size(), ref.size());
}

void trimming_pairwise_alignment(wfa::WFAligner& aligner, const organo_opts& params, std::vector<abg>& abg_reads, std::vector<std::string>& seqs, andistmatrix& distmatrix)
{
    for(uint32_t i = 0; i < seqs.size(); ++i){
    	abg& subj_tag = abg_reads[i];
    	if(!subj_tag.spanning() && (subj_tag.spanning_l || subj_tag.spanning_r)){
    		std::string& subj = seqs[i];
    		for(uint32_t j = 0; j < seqs.size(); ++j){
				abg& target_tag = abg_reads[j];
				if(target_tag.spanning()){
					std::string& target = seqs[j];
		    		if(subj.size() > target.size()){
		    			uint32_t max_size = subj.size() > target.size() ? subj.size() : target.size();
			            double norm_dist;
			        	//no need for alignment
			        	if(subj == target) norm_dist = 0.0;
			        	else if(subj.size() == 0 || target.size() == 0) norm_dist = 1.0;
			        	else {
			        		//swap target, subj given situation (flank(s) is there, just not aligned)
	                		spanning_alignment(aligner, target, subj_tag, subj);
			        		norm_dist = aligner.getAlignmentScore() / (double)max_size;
			        	}
			        	distmatrix.set_dist(i, j, norm_dist);
		    		}
				}
	        }
    	}
    }
}

void flanking_kmers(const organo_opts& params, const std::string& cigar, std::vector<int>& kmers)
{
	int pos = 0;
	int k_acc = 0;

	for(uint32_t i = 0; i < cigar.size(); ++i){
		if(k_acc >= params.k) kmers.emplace_back(pos);
		auto& op = cigar[i];
		if(op == 'M') ++k_acc;
		else k_acc = 0;
		if(op != 'D') ++pos;
	}
	if(k_acc >= params.k) kmers.emplace_back(pos);

}

void find_bps_cigar(wfa::WFAligner& aligner, const organo_opts& params, std::vector<abg>& abg_reads, std::vector<std::string>& seqs, int& subj_i, std::vector<int>& target_indeces, std::vector<int>& bps_l, std::vector<int>& bps_r)
{
	abg& subj_tag = abg_reads[subj_i];
	std::string& subj_seq = seqs[subj_i];
	for(const auto& j : target_indeces){
		auto& target_seq = seqs[j];
		if(subj_seq.size() > target_seq.size()){
			std::string subj_seq_sub;
			if(subj_tag.spanning_r) subj_seq_sub = subj_seq.substr(subj_seq.size() - target_seq.size());
			else if(subj_tag.spanning_l) subj_seq_sub = subj_seq.substr(0, target_seq.size());
			if(subj_seq_sub == target_seq) {
				if(subj_tag.spanning_r) {
					bps_l.emplace_back(subj_seq.size() - subj_seq_sub.size() + params.k);
					bps_r.emplace_back(subj_seq.size());
				}
				else if(subj_tag.spanning_l) {
					bps_l.emplace_back(params.k);
					bps_r.emplace_back(subj_seq_sub.size());
				}
			}
			else{
				
				aligner.alignEnd2End(target_seq, subj_seq);
				std::string cigar = aligner.getAlignmentCigar();
				int all = 0;
				int match = 0;
				for(char& op : cigar) {
					if(op != 'I') ++all;
					if(op == 'M') ++match;
				}
				double sim = (double)match / all;
				if(sim > params.minsim){
					std::vector<int> kmers;
					flanking_kmers(params, cigar, kmers);
					if(!kmers.empty()){
						bps_l.emplace_back(kmers.front());
						bps_r.emplace_back(kmers.back());
					}
				}
			}
		}
	}
}

int maximum_weighted_window(std::vector<uint32_t>& weights, uint32_t& k){
	uint32_t acc_weight = 0;
	for(uint32_t i = 0; i < k; ++i) acc_weight += weights[i];
	uint32_t max_weight = acc_weight;
	uint32_t max_i = 0;
	for(uint32_t i = k; i < weights.size(); ++i){
		acc_weight = acc_weight - weights[i - k] + weights[i];
		if(acc_weight > max_weight) {
			max_weight = acc_weight;
			max_i = i - k + 1;
		}
	}
	if(max_weight < k) return -1; else return (int)max_i;
}

void flank_alignment(wfa::WFAligner& aligner_gap, const organo_opts& params, abg& tag, std::string& subseq, std::string& seq, std::vector<uint32_t>& base_weight)
{
	if(tag.spanning_l) aligner_gap.alignEndsFree(subseq, 0, 0, seq, seq.size(), seq.size());
	else aligner_gap.alignEndsFree(subseq, 0, 0, seq, seq.size(), seq.size());
	
	uint32_t index = 0;
	std::string cigar = aligner_gap.getAlignmentCigar();
	uint32_t matches = 0;
	std::vector<uint32_t> weight_indeces;
	for(uint32_t i = 0; i < cigar.size(); ++i){
		const char c = cigar[i];
		if(c == 'M' || c == 'X') {
			if(c == 'M') ++matches;
			weight_indeces.emplace_back(index);
		}
		if(cigar[i] != 'D') ++index;
	}
	if(((double)matches/subseq.size()) >= params.minsim) {
		for(const auto& w : weight_indeces) ++base_weight[w];
	}
}

void realignment2(wfa::WFAligner& aligner_gap, const organo_opts& params, std::vector<abg>& abg_reads, std::vector<std::string>& seqs)
{

	uint32_t k_flank = params.k;
	uint32_t max_spanning_size = 0;
	std::vector<uint32_t> spanning_indeces;
	for(uint32_t i = 0; i < seqs.size(); ++i){
		abg& subj_tag = abg_reads[i];
		std::string& subj_seq = seqs[i];
		if(subj_tag.spanning() && subj_seq.size() >= k_flank) {
			spanning_indeces.emplace_back(i);
			if(subj_seq.size() > max_spanning_size) max_spanning_size = subj_seq.size();
		}
	}

	//std::cerr << "spanning: " << spanning_indeces.size() << '\n';

	std::vector<uint32_t> nonspanning_indeces;

	if(!spanning_indeces.empty()){
		for(uint32_t i = 0; i < seqs.size(); ++i){
			abg& subj_tag = abg_reads[i];
			std::string& subj_seq = seqs[i];
			if(!subj_tag.spanning() && subj_seq.size() > max_spanning_size && (subj_tag.spanning_l || subj_tag.spanning_r)) nonspanning_indeces.emplace_back(i);
		}
	}

	//std::cerr << "nonspanning: " << nonspanning_indeces.size() << '\n';

	if(!nonspanning_indeces.empty() && !spanning_indeces.empty()){
		std::vector<std::string> subseqs_l;
		std::vector<std::string> subseqs_r;
		for(uint32_t i = 0; i < abg_reads.size(); ++i){
			if(seqs[i].size() >= k_flank){
				const auto& tag = abg_reads[i];
				if(tag.spanning_l) subseqs_l.emplace_back(seqs[i].substr(0, k_flank));
				if(tag.spanning_r) subseqs_r.emplace_back(seqs[i].substr(seqs[i].size() - k_flank));
			}
		}

		//std::cerr << "left: " << subseqs_l.size() << '\n';
		//std::cerr << "right: " << subseqs_r.size() << '\n';

		if(!subseqs_l.empty() && !subseqs_r.empty()){
			for(const auto& i : nonspanning_indeces){
				auto& tag = abg_reads[i];
				//std::cerr << tag.name << " size: " << seqs[i].size() << '\n';
				std::vector<uint32_t> base_weight(seqs[i].size(), 0);
				if(tag.spanning_l){
					for(auto& s : subseqs_r) flank_alignment(aligner_gap, params, tag, s, seqs[i], base_weight);
				}
				else {
					for(auto& s : subseqs_l) flank_alignment(aligner_gap, params, tag, s, seqs[i], base_weight);	
				}

				int start_i = maximum_weighted_window(base_weight, k_flank);
				if(start_i >= 0){
					//std::cerr << "index: " << start_i << '\n';
					if(tag.spanning_l) {
						if((uint32_t)start_i + k_flank < seqs[i].size()) tag.realigned.reset(new std::pair<int,int>(0, start_i + (int)k_flank));
						tag.spanning_r = true;
					}
					else {
						if(start_i > 0) tag.realigned.reset(new std::pair<int,int>(start_i, (int)seqs[i].size() - start_i));
						tag.spanning_l = true;
					}
					
				}
			}
		}

	}

}


void realignment(wfa::WFAligner& aligner_edit, wfa::WFAligner& aligner_gap, const organo_opts& params, std::vector<abg>& abg_reads, std::vector<std::string>& seqs)
{
	uint32_t max_spanning_size = 0;
	std::vector<int> spanning_indeces;
	for(uint32_t i = 0; i < seqs.size(); ++i){
		abg& subj_tag = abg_reads[i];
		std::string& subj_seq = seqs[i];
		if(subj_tag.spanning() && subj_seq.size() > 0) {
			spanning_indeces.emplace_back((int)i);
			if(subj_seq.size() > max_spanning_size) max_spanning_size = subj_seq.size();
		}
	}

	std::vector<int> nonspanning_indeces;

	for(uint32_t i = 0; i < seqs.size(); ++i){
		abg& subj_tag = abg_reads[i];
		std::string& subj_seq = seqs[i];
		if(!subj_tag.spanning() && subj_seq.size() > max_spanning_size && (subj_tag.spanning_l || subj_tag.spanning_r)) nonspanning_indeces.emplace_back((int)i);
	}

	if(!nonspanning_indeces.empty() && !spanning_indeces.empty()){
		for(auto& i : nonspanning_indeces){
			abg& subj_tag = abg_reads[i];
			std::vector<int> bps_l;
			std::vector<int> bps_r;
			find_bps_cigar(aligner_gap, params, abg_reads, seqs, i, spanning_indeces, bps_l, bps_r);
			if(!bps_l.empty()){
				std::vector<std::pair<int,int>> bps_l_counts;
				group_count_task(bps_l, bps_l_counts);
				std::vector<std::pair<int,int>> bps_r_counts;
				group_count_task(bps_r, bps_r_counts);
				int left = bps_l_counts.front().first - params.k;
				int right = bps_r_counts.front().first;
				if(left != 0 || right != (int)seqs[i].size()){
					subj_tag.realigned.reset(new std::pair<int,int>(left, right - left));
					if(subj_tag.spanning_l) subj_tag.spanning_r = true; else subj_tag.spanning_l = true;
				}
			}
		}
	}

}

void spanning_aware_pairwise_alignment(wfa::WFAligner& aligner, const organo_opts& params, std::vector<abg>& abg_reads, std::vector<std::string>& seqs, std::vector<uint32_t>& seqs_l, andistmatrix& distmatrix)
{
    for(uint32_t i = 0; i < abg_reads.size(); ++i){
    	abg& subj_tag = abg_reads[i];
        std::string& subj = seqs[i];
        const uint32_t& subj_l = seqs_l[i];
        for(uint32_t j = i + 1; j < abg_reads.size(); ++j){
			abg& target_tag = abg_reads[j];
            std::string& target = seqs[j];
            const uint32_t& target_l = seqs_l[j];
            //at least one spanning
            if(!subj_tag.spanning() && !target_tag.spanning()) distmatrix.set_dist(i, j, 1.0);
            else {
            	//original size difference
            	uint32_t min_size = subj_l < target_l ? subj_l : target_l;
				uint32_t max_size = subj_l > target_l ? subj_l : target_l;
        		uint32_t size_dist = max_size - min_size;
        		double norm_dist_o = size_dist / (double)max_size;
            	//no need for alignment
	        	if(size_dist == 0 && subj == target) distmatrix.set_dist(i, j, 0.0);
	        	//empty string
	        	else if(max_size == 0 || min_size == 0) distmatrix.set_dist(i, j, 1.0);
	        	else {
					max_size = subj.size() > target.size() ? subj.size() : target.size();
	        		double norm_dist;
	        		//both are spanning
	        		if(subj_tag.spanning() && target_tag.spanning()){
	        			spanning_alignment(aligner, target, target_tag, subj);
	        			double norm_dist_e = std::abs(aligner.getAlignmentScore()) / (double)max_size;
	                	norm_dist = norm_dist_e > norm_dist_o ? norm_dist_e : norm_dist_o;
	        		}
	                //subj is spanning
	                else {
	                	if(subj_tag.spanning() && !target_tag.spanning()) spanning_alignment(aligner, target, target_tag, subj);
	                	//target is spanning
	                	else spanning_alignment(aligner, subj, subj_tag, target);
	                	norm_dist = std::abs(aligner.getAlignmentScore()) / (double)max_size;
	            	}

	                distmatrix.set_dist(i, j, norm_dist);
	       		}
            }
        }
    }
}

void nonspanning_assignment_task(const organo_opts& params, andistmatrix& distmatrix, const std::vector<abg>& abg_reads, const std::vector<std::vector<uint32_t>>& ccs, std::vector<std::vector<uint32_t>>& ccs_expanded, std::vector<uint32_t>& unassigned)
{
	for(uint32_t i = 0; i < abg_reads.size(); ++i){
		if(!abg_reads[i].spanning()){
			std::vector<uint32_t> valid_ccs;
			//iterate through each cc and check if seq is assignable (error < threshold)
			for(uint32_t cc_i = 0; cc_i < ccs.size(); ++cc_i){
				double min_weight = -1.0;
				int min_weight_i = -1;
				auto& local_cc = ccs[cc_i];
				for(const auto& c : local_cc){
					double dist = distmatrix.get_dist(i, c);
					if(dist <= params.maxerror && (min_weight_i == -1 || dist < min_weight)){
						min_weight = dist;
						min_weight_i = (int)c;
					}
				}
				//is assignable to current cc
				if(min_weight >= 0.0) valid_ccs.emplace_back(cc_i);
			}
			//unique assignment
			if(valid_ccs.empty()) unassigned.emplace_back(i);
			else if(valid_ccs.size() == 1) ccs_expanded[valid_ccs.front()].emplace_back(i);
		}
	}
}

