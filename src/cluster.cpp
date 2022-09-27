#include <iostream>
#include <vector>
#include <htslib/sam.h>
#include <thread>
#include <chrono>
#include <fstream>
#include "antimestamp.hpp"
#include "BS_thread_pool.hpp"
#include "abg.hpp"
#include "organo_opts.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "cluster.hpp"
#include "anbamtools.hpp"
#include "utils.hpp"
#include "andistmatrix.hpp"
#include "pairwise_alignments.hpp"
#include "clusteringtools.hpp"
#include "ansparc.hpp"

void output_seq_process(const std::string& region, std::vector<abg>& abg_reads, const uint32_t& total_reads, const uint32_t& cc_id, const std::vector<uint32_t>& cc_complete, const std::vector<uint32_t>& ccs_expanded, std::string& seq)
{
	uint32_t dp = cc_complete.size() + ccs_expanded.size();
	double freq = (double)dp / total_reads;
	std::cout << '>' << region << " CC:i:" << cc_id << " DP:i:" << dp << " AF:f:" << freq << " RN:Z:";
	for(uint32_t i = 0; i < cc_complete.size(); ++i){
		if(i != 0) std::cout << ';';
		std::cout << abg_reads[cc_complete[i]].name;
	}
	for(const auto& i : ccs_expanded) std::cout << ';' << abg_reads[i].name;
    std::cout << '\n' << seq << '\n';
}

void consensus_process(wfa::WFAlignerGapAffine2Pieces& aligner, const uint32_t& representative_i, std::vector<std::string>& sequences, const std::vector<uint32_t>& cc, std::string& consensus)
{
    if(cc.size() < 3) consensus = sequences[representative_i];
    else{
        uint32_t k = 3;
        uint32_t g = 2;
        uint32_t l = 2;
        float c = 2.0f;
        if(sequences.size() < 4) c = 1.0f;
        float t = 0.2;
        ansparc_graph graph;
        graph.init(sequences[representative_i], k, g, l);
        for(const auto& i : cc){
            if(i != representative_i){
                aligner.alignEnd2End(sequences[representative_i], sequences[i]);
                std::string cigar = aligner.getAlignmentCigar();
                graph.insert_alignment(cigar, sequences[i]);
            }
        }
        //adjust weights
        graph.adjust_weights(c, t);
        graph.consensus(consensus);
    }
}

void genotype_process(
	BS::thread_pool& pool, 
	std::vector<abg_block>& loaded_blocks,
	std::mutex& seq_block_mutex,
	std::ofstream& dist_output,
	const organo_opts& params
	)
{
	pool.parallelize_loop(0, loaded_blocks.size(),
	[&pool, &loaded_blocks, &seq_block_mutex, &params, &dist_output](const int a, const int b){
		wfa::WFAlignerEdit aligner_edit(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
		wfa::WFAlignerGapAffine2Pieces aligner_gap(2, 4, 2, 24, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryMed);

		for(int i = a; i < b; ++i){
			organo_opts local_params = params;
			//init alignment objects for current thread
			auto& local_block = loaded_blocks[i];
			std::vector<uint32_t> seqs_l(local_block.size());
			for(uint32_t j = 0; j < local_block.size(); ++j) seqs_l[j] = local_block.seqs[j].size();
			std::vector<std::string> hp_compressed_seq;
	    	if(!params.hp && params.hpd) {
	    		hp_compressed_seq.resize(local_block.size(), "");
	        	for(uint32_t j = 0; j < local_block.seqs.size(); ++j) if(local_block.seqs[j].size() > 0) homopolymer_compressor(local_block.seqs[j].c_str(), hp_compressed_seq[j]);
	      	}

			andistmatrix distmat(local_block.size());
			if(!params.hp && params.hpd) spanning_aware_pairwise_alignment(aligner_edit, params, local_block.reads, hp_compressed_seq, seqs_l, distmat);
			else spanning_aware_pairwise_alignment(aligner_edit, params, local_block.reads, local_block.seqs, seqs_l, distmat);
			std::vector<uint32_t> spanning;
			std::vector<uint32_t> spanning_nonempty;
			for(uint32_t s = 0; s < local_block.size(); ++s) {
				if(local_block.reads[s].spanning()) {
					spanning.emplace_back(s);
					if(local_block.seqs[s].size() > 0) spanning_nonempty.emplace_back(s);
				}
			}
			if(params.maxalleles > 0){
				double error_est = distmat.binned_kde(3, params.mincov > 1 ? params.mincov - 1 : 1, params.maxerror, spanning_nonempty);
				local_params.maxerror = error_est > params.maxerror ? error_est : params.maxerror;
			}
			std::vector<std::vector<uint32_t>> ccs;
			cluster_task(local_params, distmat, spanning, ccs);
			std::vector<std::vector<uint32_t>> ccs_expanded(ccs.size());
			std::vector<uint32_t> unassigned;
			nonspanning_assignment_task(local_params, distmat, local_block.reads, ccs, ccs_expanded, unassigned);
			uint32_t total_reads = 0.0;
			for(uint32_t j = 0; j < ccs.size(); ++j) total_reads += ccs[j].size() + ccs_expanded[j].size();
			for(uint32_t j = 0; j < ccs.size(); ++j) {
				auto& cc = ccs[j];
				auto& cc_e = ccs_expanded[j];
				uint32_t rep_i = distmat.representative(cc);
				std::string consensus;
				consensus_process(aligner_gap, rep_i, local_block.seqs, cc, consensus);
				seq_block_mutex.lock();
				output_seq_process(local_block.name, local_block.reads, total_reads, j, cc, cc_e, consensus);
				seq_block_mutex.unlock();
			}
		}
	}).wait();
}


void cluster(const std::string& bam, const organo_opts& params)
{
	//setup thread pool
 	BS::thread_pool pool(params.threads);
 	
 	std::cerr << "(" << antimestamp() << "): Initialised " << params.threads << " additional threads " << std::endl;

	std::ofstream dist_output;
  	if(params.tsv.size() > 0) dist_output.open(params.tsv);
  	else if(params.matrix.size() > 0) dist_output.open(params.matrix);

	uint32_t total_blocks_processed = 0;
	std::mutex seq_block_mutex;

	abg_block_iter iter(bam, params.block_size, params.maxcov, params.hp);

	while(iter.next()){
		std::cerr << "(" << antimestamp() << "): Loaded " << params.block_size << " target regions\n";
		genotype_process(pool, iter.loaded_blocks, seq_block_mutex, dist_output, params);
		total_blocks_processed += iter.loaded_blocks.size();
		std::cerr << "(" << antimestamp() << "): Processed " << total_blocks_processed << " target regions\n";
	}

  	iter.destroy();


	if(params.tsv.size() > 0 || params.matrix.size() > 0) dist_output.close();
}