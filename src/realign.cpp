#include <vector>
#include <htslib/sam.h>
#include <thread>
#include <chrono>
#include <memory>
#include <utility>
#include "realign.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "antimestamp.hpp"
#include "BS_thread_pool.hpp"
#include "abg.hpp"
#include "organo_opts.hpp"
#include "anbamtools.hpp"
#include "utils.hpp"
#include "andistmatrix.hpp"
#include "pairwise_alignments.hpp"


void realign_process(BS::thread_pool& pool, samFile* outfile, std::vector<abg_block>& loaded_blocks, std::mutex &seq_block_mutex, const organo_opts& params)
{
	pool.parallelize_loop(0, loaded_blocks.size(),
	[&pool, &outfile, &loaded_blocks, &seq_block_mutex, &params](const int a, const int b){
		//wfa::WFAlignerGapAffine2Pieces aligner(2, 4, 2, 24, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryMed);
		wfa::WFAlignerGapAffine aligner(4,6,2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryMed);
		wfa::WFAlignerEdit aligner_edit(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);

		for(int i = a; i < b; ++i){
			//init alignment objects for current thread
			auto& local_block = loaded_blocks[i];
			//realignment(aligner_edit, aligner, params, local_block.reads, local_block.seqs);
			realignment2(aligner, params, local_block.reads, local_block.seqs);
			std::vector<uint32_t> final_indeces(local_block.size());
			for(uint32_t j = 0; j < local_block.size(); ++j) final_indeces[j] = j;
			seq_block_mutex.lock();
			write2bam(outfile, local_block.name, final_indeces, local_block.reads, local_block.seqs);
			seq_block_mutex.unlock();
		}
	}).wait();
}

void realign(const std::string& bam, const organo_opts& params)
{
	//setup thread pool
 	BS::thread_pool pool(params.threads);
 	
 	std::cerr << "(" << antimestamp() << "): Initialised " << params.threads << " additional threads " << std::endl;

 	uint32_t total_blocks_processed = 0;
 	std::mutex seq_block_mutex;

 	//ansbam bam_inst;
 	//bam_inst.init(bam, false);

 	samFile* file = sam_open("-", "wb" );
  bam_hdr_t *hdr = bam_hdr_init();
  hdr->l_text = strlen(init_header);
  hdr->text = strdup(init_header);
  hdr->n_targets = 0;
  sam_hdr_write(file, hdr);

  abg_block_iter iter(bam, params.block_size, params.maxcov, params.hp);

  while(iter.next()){
  	std::cerr << "(" << antimestamp() << "): Loaded " << params.block_size << " target regions\n";
  	realign_process(pool, file, iter.loaded_blocks, seq_block_mutex, params);
  	total_blocks_processed += iter.loaded_blocks.size();
  	std::cerr << "(" << antimestamp() << "): Processed " << total_blocks_processed << " target regions\n";
  }

  iter.destroy();

  /**
  abg_block current_block;
  std::vector<abg_block> loaded_blocks;

  while(sam_read1(bam_inst.fp, bam_inst.header, bam_inst.read) > 0){
  	std::string local_name; //region
  	abg abg_read(bam_inst.read, local_name);
    if(current_block.name.empty()) current_block.name = local_name;
    std::string local_seq = "";
	  //quality string
		uint8_t *q = bam_get_seq(bam_inst.read);
		//update array for read seqeunce
		for(int i=0; i < bam_inst.read->core.l_qseq ; i++) local_seq += seq_nt16_str[bam_seqi(q,i)];

  	//compress if user provided
  	if(params.hp && local_seq.size() > 0) {
  		std::string compressed = "";
  		homopolymer_compressor(local_seq.c_str(), compressed);
  		local_seq = compressed;
  	}
  	//still sequences to load
  	if(local_name == current_block.name) {
  		current_block.seqs.emplace_back(local_seq);
  		current_block.reads.emplace_back(abg_read);
  	}
  	//finished current block
  	else{
  		if(current_block.size() <= params.maxcov) loaded_blocks.emplace_back(current_block);
  		if(loaded_blocks.size() > params.block_size - 1){
  			trimming_process(pool, file, loaded_blocks, seq_block_mutex, params);
    		total_blocks_processed += loaded_blocks.size();
    		std::cerr << "(" << antimestamp() << "): Processed " << total_blocks_processed << " target regions\n";
    		loaded_blocks.clear();
  		}
  		current_block.clear();
  		current_block.name = local_name;
  		current_block.seqs.emplace_back(local_seq);
  		current_block.reads.emplace_back(abg_read);
  	}
  }

  if(current_block.size() <= params.maxcov) loaded_blocks.emplace_back(current_block);

  if(!loaded_blocks.empty()) {
  	trimming_process(pool, file, loaded_blocks, seq_block_mutex, params);
  	total_blocks_processed += loaded_blocks.size();
    std::cerr << "(" << antimestamp() << "): Processed " << total_blocks_processed << " target regions\n";
    loaded_blocks.clear();
  }
  */

	bam_hdr_destroy(hdr);
  sam_close(file); // close bam file
}