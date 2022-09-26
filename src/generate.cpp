#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include <iostream>
#include <vector>
#include <htslib/sam.h>
#include <thread>
#include <chrono>
#include <memory>
#include <map>
#include "generate.hpp"
#include "antimestamp.hpp"
#include "BS_thread_pool.hpp"
#include "anbamtools.hpp"
#include "generate.hpp"
#include "organo_opts.hpp"
#include "tsv_helper.hpp"
#include "anbed.hpp"
#include "abg.hpp"

void duplicates(const std::vector<abg>& abg_reads, const std::vector<uint32_t>& sorted_indeces, std::vector<std::pair<uint32_t,uint32_t>>& duplicate_indeces)
{

	uint32_t last = 0;
	uint32_t i = last + 1;
	for(; i < abg_reads.size(); ++i){
		if(abg_reads[sorted_indeces[last]].name != abg_reads[sorted_indeces[i]].name){
			duplicate_indeces.emplace_back(std::make_pair(last, i));
			last = i;
		}
	}
	if(last < abg_reads.size()) duplicate_indeces.emplace_back(std::make_pair(last, i));
}

void generate(const std::string bam, const std::string bed, const std::string haplotypes, const organo_opts& params)
{
	//setup thread pool
 	BS::thread_pool pool(params.threads);
 	//thread -> index of allocated object-instances
 	std::map<std::thread::id, int> thread2index;
 	//thread index lock
 	std::mutex thread_int_mutex;
 	int thread_int = 0;
 	//run all threads and log thread id -> index
 	pool.parallelize_loop(0, params.threads,[&thread_int_mutex, &thread_int, &pool, &thread2index](const int a, const int b){
 		//sleep or else not all specified threads will be logged
 		std::this_thread::sleep_for(std::chrono::milliseconds(500));
 		thread_int_mutex.lock();
 		thread2index.insert({ std::this_thread::get_id(), thread_int });
 		++thread_int;
 		thread_int_mutex.unlock();
 	}).wait();
 	std::cerr << "(" << antimestamp() << "): Initialised " << params.threads << " additional threads " << std::endl;

 	//initiate bam file openers for each thread
 	std::vector<ansbam> bam_insts(thread2index.size(), ansbam());
 	for(uint32_t i = 0; i < thread2index.size(); ++i) bam_insts[i].init(bam, true);

 	//set lock for outputing sequencing block to stdout
 	std::mutex seq_block_mutex;

 	//load bed file
	BEDmap bedmap;
	std::vector<BED> bed_regions;
	//parse bed file and store results in variable above
	bedparser(bed, bedmap, 0);
	//set total regions observed
    int total_regions = 0;
    int total_contigs = 0;
    //traverse and count total regions
    auto it = bedmap.begin();
    while(it != bedmap.end()){
    	sort(it->second.begin(), it->second.end(), sortbed);
    	uniquebed(it->second);
    	total_regions += it->second.size();
    	++total_contigs;
    	for(const auto& bed : it->second) bed_regions.emplace_back(bed);
    	it = bedmap.erase(it);
    }

    std::cerr << '(' << antimestamp() << "): Found " << total_contigs << " total contigs and " << total_regions << " total regions of interest\n";

    //load haplotype-tagged reads, if provided
    std::unordered_map<std::string, int> read2ha;
    if(haplotypes.size() > 0){
    	load_ha_tagged_reads(haplotypes, read2ha);
    	std::cerr << "(" << antimestamp() << "): Loaded " << read2ha.size() << " total haplotype-tagged reads" << std::endl;
    }

    //set-upt bam file opener to stdout
    samFile* bamstdout = sam_open("-", "wb" );
    bam_hdr_t *hdr = bam_hdr_init();
    hdr->l_text = strlen(init_header);
    hdr->text = strdup(init_header);
    hdr->n_targets = 0;
    sam_hdr_write(bamstdout, hdr);

 	std::cerr << "(" << antimestamp() << "): Processing " << bam << std::endl;
	pool.parallelize_loop(0, bed_regions.size(),
		[&pool, &thread2index, &bam_insts, &bed_regions, &read2ha, &params, &seq_block_mutex, bamstdout](const int a, const int b){
			const std::string sa_tag = "SA";
			const std::string np_tag = "np";
			std::thread::id thread_id = std::this_thread::get_id();
			int thread_index = thread2index.find(thread_id)->second;
			ansbam& bam_inst = bam_insts[thread_index];

			for (int i = a; i < b; ++i){
				auto& region = bed_regions[i];
				std::string region_bed = region.toBEDstring();
				hts_itr_t *iter = sam_itr_querys(bam_inst.idx, bam_inst.header, region_bed.c_str());
				if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " << region_bed << std::endl;
				else{
					std::vector<abg> abg_reads;
					std::vector<std::string> abg_edge_seqs;
					//std::cerr << "(" << antimestamp() << "): Iterating through reads at " << region_bed << std::endl;
					//iterate through each alignment and generate AB-graphs
					while(sam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0){
						if(bam_inst.read->core.qual >= params.min_mapq && !(bam_inst.read->core.flag & BAM_FSECONDARY || bam_inst.read->core.flag & BAM_FSUPPLEMENTARY)){
							std::string edge_seq;
							abg_generate_msg msg;
							abg_generate(bam_inst.read, region.start, region.end, msg, edge_seq);
							//std::cerr << "(" << antimestamp() << "): --Generated " << region_bed << std::endl;
							if(msg.successful) {
								//std::cerr << "(" << antimestamp() << "): --Success " << region_bed << std::endl;
								if(params.fasta){
									char span_tag;
									if(msg.spanning_l && msg.spanning_r) span_tag = 'b';
									else if(msg.spanning_l) span_tag = 'l';
									else if(msg.spanning_r) span_tag = 'r';
									else span_tag = 'n';

									seq_block_mutex.lock();
									std::cout << '>' << region_bed << " SP:Z:" << span_tag << " RN:Z:" << (char*)bam_inst.read->data << '\n';
									std::cout << edge_seq << '\n';
									seq_block_mutex.unlock();
								}
								else{
									abg read;
									read.name = (char*)bam_inst.read->data;
									read.spanning_l = msg.spanning_l;
									read.spanning_r = msg.spanning_r;
								    auto sa = bam_aux_get(bam_inst.read, sa_tag.c_str());
								    read.supplement = sa != NULL ? true : false;
								    auto np = bam_aux_get(bam_inst.read, np_tag.c_str());
								    if(np != NULL) read.np = std::make_shared<int>(bam_aux2i(np));
									auto it_ha = read2ha.find(read.name);
									if(it_ha != read2ha.end()) read.haplotype = std::make_shared<int>(it_ha->second);
									abg_reads.emplace_back(read);
									abg_edge_seqs.emplace_back(edge_seq);
								}
								//std::cerr << "(" << antimestamp() << "): --Added " << region_bed << std::endl;
							}
						
						}
					}

					if(!abg_reads.empty()){
						//std::cerr << "(" << antimestamp() << "): Found " << abg_reads.size() << " reads " << std::endl;
						std::vector<uint32_t> indeces(abg_reads.size());
						for(uint32_t i = 0; i < abg_reads.size(); ++i) indeces[i] = i;
						sort(indeces.begin(), indeces.end(), [&abg_reads](const uint32_t& i, const uint32_t& j) -> bool
				        { 
				        	abg& i_tag = abg_reads[i];
				        	abg& j_tag = abg_reads[j];
				        	//no haplotype info available for both
				            if(!i_tag.haplotype && !j_tag.haplotype) return i_tag.name < j_tag.name;
				            else if(!i_tag.haplotype && j_tag.haplotype) return false;
				            else if(i_tag.haplotype && !j_tag.haplotype) return true;
				            //same haplotype
				            else if(*i_tag.haplotype == *j_tag.haplotype) return i_tag.name < j_tag.name;
				            //different haplotype
				            else return *i_tag.haplotype < *j_tag.haplotype;
				        });
						/**
						//std::cerr << "(" << antimestamp() << "): Sorted " << std::endl;
						std::vector<std::pair<uint32_t,uint32_t>> duplicate_indeces;
						duplicates(abg_reads, indeces, duplicate_indeces);
						//std::cerr << "(" << antimestamp() << "): Mark multimapped" << std::endl;
						std::vector<uint32_t> final_indeces;
						for(auto it = duplicate_indeces.begin(); it != duplicate_indeces.end();++it){
							if((int)it->second - (int)it->first < 1) exit(1);
							if(it->second - it->first == 1) final_indeces.emplace_back(indeces[it->first]);
							//if ambiguous
							else {
								uint32_t primary_i = indeces.size();
								for(uint32_t i = it->first; i < it->second; ++i){
									if(primaries[indeces[i]]){
										primary_i = i;
										break;
									}
								}
								if(primary_i < indeces.size()) final_indeces.emplace_back(indeces[primary_i]);
							}
						}
						*/
						//std::cerr << "(" << antimestamp() << "): Reduced" << std::endl;
						seq_block_mutex.lock();
						write2bam(bamstdout, region_bed, indeces, abg_reads, abg_edge_seqs);
						seq_block_mutex.unlock();
						//std::cerr << "(" << antimestamp() << "): Finished" << std::endl;
					}
				}
				hts_itr_destroy(iter);
			}
 	}).wait();
 	bam_hdr_destroy(hdr);
    sam_close(bamstdout); // close bam file
	 	
	for(uint32_t i = 0; i < thread2index.size(); ++i) bam_insts[i].destroy();

}