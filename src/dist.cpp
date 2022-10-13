#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include <iostream>
#include <vector>
#include <htslib/sam.h>
#include <unordered_map>
#include <utility>
#include <string>
#include "antimestamp.hpp"
#include "kseq.h"
#include "BS_thread_pool.hpp"
#include "organo_opts.hpp"
#include "anbed.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "utils.hpp"

KSEQ_INIT(gzFile, gzread)

struct target_seqs {
	std::vector<std::string> seqs;
	std::vector<int> indeces;
};

void get_allele_info(const organo_opts& params, const char* comment_ptr, int& cov, double& af)
{
	if(comment_ptr != nullptr) {
		std::string value;
    	std::istringstream stream(comment_ptr);
    	while(std::getline(stream, value, ' ')) {
    		if(value.substr(0,5) == "DP:i:") cov = std::stoi(value.substr(5));
    		else if(value.substr(0,5) == "AF:f:") af = std::stod(value.substr(5));
    	}
	}
}

void dist( const std::vector<std::string>& fastas, const organo_opts& params)
{
	uint32_t total_inputs = fastas.size();
	//setup thread pool
 	BS::thread_pool pool(params.threads);
 	std::mutex out_mutex;
 	
 	std::cerr << "(" << antimestamp() << "): Initialised " << params.threads << " additional threads " << std::endl;

 	std::unordered_map<std::string, target_seqs> loaded_seqs;

 	KS_FULL_COMMENT = true;
 	for(uint32_t fasta_i = 0; fasta_i < total_inputs; ++fasta_i){
 		const std::string& fasta = fastas[fasta_i];
 		//LOAD FASTA
	    gzFile fp = gzopen(fasta.c_str(), "r");
	    //sequence pointer
	    kseq_t *seq = kseq_init(fp);
	    int seq_l;
	    while ((seq_l = kseq_read(seq)) >= 0 ) {
	    	int cov = params.mincov;
	    	double af = params.minaf;
	    	get_allele_info(params, seq->comment.s, cov, af);
	    	if(fasta_i == 0 || (cov >= params.mincov && af >= params.minaf)){
	    		std::string region = seq->name.s;
		    	std::string local_seq = "";

		    	if(seq_l > 0) local_seq = seq->seq.s;
		    	//compress if user provided
		    	if(params.hp) {
		    		std::string compressed;
		    		homopolymer_compressor(local_seq.c_str(), compressed);
		    		local_seq = compressed;
		    	}

		    	auto it = loaded_seqs.find(region);
		    	if(fasta_i == 0 && it == loaded_seqs.end()){
		    		target_seqs new_target;
		    		new_target.seqs.emplace_back(local_seq);
		    		new_target.indeces.emplace_back(fasta_i);
		    		loaded_seqs.insert({region, new_target});
		    	}
		    	else if(it != loaded_seqs.end()){
		    		it->second.seqs.emplace_back(local_seq);
		    		it->second.indeces.emplace_back(fasta_i);
		    	}

	    	}
	    }
	  	kseq_destroy(seq);
	  	gzclose(fp);
 	}
 	//move to a vector STL
 	std::vector<std::pair<std::string, target_seqs>> loaded_seqs_vect;
 	auto it = loaded_seqs.begin();
 	while(it != loaded_seqs.end()){
 		loaded_seqs_vect.emplace_back(std::make_pair(it->first, it->second));
 		it = loaded_seqs.erase(it);
 	}
 	std::cerr <<  '(' << antimestamp() << "): Loaded " << loaded_seqs_vect.size() << " total target sequences\n";
 	
 	pool.parallelize_loop(0, loaded_seqs_vect.size(),
	[&pool, &out_mutex, &total_inputs, &loaded_seqs_vect, &params](const int a, const int b){
		//init alignment objects for current thread
		wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
		for(int i = a; i < b; ++i){
			auto& local_block = loaded_seqs_vect[i];
			auto bed_region = BED(local_block.first);
			std::vector<std::pair<int,int>> group_counts;
			group_count_task(local_block.second.indeces, group_counts);
			std::vector<int> indeces_counts(total_inputs, 0);
			bool found = false;
			uint32_t max_size = 0;
			for(uint32_t j = 0; j < local_block.second.seqs.size(); ++j){
				if(local_block.second.indeces[j] != 0){
					found = true;
					auto& subj = local_block.second.seqs[j];
					double min_dist = 1.0;
					if(subj.size() > max_size) max_size = subj.size();
					for(uint32_t k = 0; k < local_block.second.seqs.size(); ++k){
						if(j != k && local_block.second.indeces[k] == 0) {
							double dist;
							auto& ref = local_block.second.seqs[k];
							if(subj == ref) dist = 0.0;
							else{
								aligner.alignEnd2End(subj, ref);
								uint32_t max_size = subj.size() > ref.size() ? subj.size() : ref.size();
								dist = aligner.getAlignmentScore() / (double)max_size;
							}
							if(dist < min_dist) min_dist = dist;
							if(ref.size() > max_size) max_size = ref.size();
						}
					}
					out_mutex.lock();
					std::cout << bed_region.chr << ':' << bed_region.start << '-' << bed_region.end << '\t' << bed_region.size() << '\t' << max_size << '\t' << local_block.second.indeces[j] << '\t' << min_dist << '\n';
					out_mutex.unlock();

				}
			}

			if(!found){
				for(uint32_t j = 0; j < local_block.second.seqs.size(); ++j) if(local_block.second.indeces[j] == 0) if(local_block.second.seqs[j].size() > max_size) max_size = local_block.second.seqs[j].size();
				for(uint32_t j = 0; j < local_block.second.seqs.size(); ++j){
					if(local_block.second.indeces[j] == 0){
						out_mutex.lock();
						std::cout << bed_region.chr << ':' << bed_region.start << '-' << bed_region.end << '\t' << bed_region.size() << '\t' << max_size << "\t0\t1.0\n";
						out_mutex.unlock();
					}
				}
			}
		}
	}).wait();
}