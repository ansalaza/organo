#include <htslib/sam.h>
#include <memory>
#include <vector>
#include <limits>
#include <string>
#include <iostream>
#include "anbamtools.hpp"
#include "utils.hpp"
#include "abg.hpp"

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2))

ansbam::ansbam(){};

void ansbam::init(const std::string file, bool _index)
{
	index = _index;
	fp = hts_open(file.c_str(),"r");
	header = sam_hdr_read(fp);
	//load index
	if(index) idx = sam_index_load(fp, file.c_str());
	//init read
	read = bam_init1();
};

void ansbam::destroy()
{
	sam_close(fp);
	fp = nullptr;
	if(index) hts_idx_destroy(idx);
	idx = nullptr;
	bam_hdr_destroy(header);
	header = nullptr;
	bam_destroy1(read);
	read = nullptr;
};

void project_positions(bam1_t* alignment, bool& clipped_l, bool& clipped_r, std::vector<int>& refcoords)
{
	refcoords.resize(alignment->core.l_qseq, -1);
	uint32_t *cigar = bam_get_cigar(alignment);
	//current positions in ref and query
	int rpos = alignment->core.pos;
	int qpos = 0;
	for (uint32_t i = 0; i < alignment->core.n_cigar; ++i) { 
        const int op = bam_cigar_op(cigar[i]);
		const int ol = bam_cigar_oplen(cigar[i]);
		if(op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP){
			if(i == 0) clipped_l = true;
			if(i == alignment->core.n_cigar - 1) clipped_r = true;
			if(op == BAM_CSOFT_CLIP) qpos += ol;
		}
		else if(op == BAM_CMATCH || op ==  7 || op == 8){
			for(int j = 0; j < ol; ++j){
				refcoords[qpos] = rpos;
				++rpos;
				++qpos;
			}
		}
		else if(op == BAM_CINS) qpos += ol;
		else if(op == BAM_CDEL) rpos += ol;
	}
}

void get_breakpoints(const int& start,  const int& end, const bool& clipped_l, const bool& clipped_r, abg_generate_msg& msg, std::vector<int>& projectedcoords, std::unique_ptr<std::pair<int, int>>& ptr)
{
	int leftmost = -1, rightmost = -1, qstart = -1, qend = -1, qstart_dist = -1, qend_dist = -1;

	for(uint32_t i = 0; i < projectedcoords.size(); ++i){
		//current coordinate in reference
		int ccoord = projectedcoords[i];
		if(ccoord > -1){
			//update the left-most index
			if(leftmost == -1) leftmost = i;
			//update the right most index
			if(rightmost == -1 || ccoord > projectedcoords[rightmost]) rightmost = i;
			//distance to start coordinate 
			int cstart_dist = ccoord - start;
			//distance to end coordinate
			int cend_dist = end - ccoord;
			//new closest start coordinate found
			if(cstart_dist >= 0 && (qstart_dist < 0 || cstart_dist < qstart_dist)){
				qstart_dist = cstart_dist;
				qstart = i;
			}
			//new closest end coordinate found
			if(cend_dist >= 0 && (qend_dist < 0 || cend_dist < qend_dist)){
				qend_dist = cend_dist;
				qend = i;
			} 
		}
	}
	if(projectedcoords[rightmost] < start || projectedcoords[leftmost] > end){
		qstart = -1;
		qend = -1;
		msg.successful = false;
		msg.spanning_l = false;
		msg.spanning_r = false;
	}
	//region deleted
	else if(qstart > -1 && qend > -1 && qstart > qend){
		qstart = -1;
		qend = -1;
		msg.successful = true;
		msg.spanning_l = true;
		msg.spanning_r = true;
	} 
	else{
		//readjust if alignment is clipped
		if(clipped_l && qstart > 0) while(projectedcoords[qstart - 1] == -1) --qstart;
		if(clipped_r && qend > 0 && qend < (int)projectedcoords.size() - 1) while(projectedcoords[qend + 1] == -1) ++qend;

		if(leftmost >= 0 && projectedcoords[leftmost] <= start && !(qstart == 0 && clipped_l)) msg.spanning_l = true; else msg.spanning_l = false;
		if(rightmost >= 0 && projectedcoords[rightmost] >= end && !(qend == (int)projectedcoords.size() - 1 && clipped_r)) msg.spanning_r = true; else msg.spanning_r = false;
		msg.successful = true;
	}

	if(msg.successful){
		//spanning both sides
		if(msg.spanning_l && msg.spanning_r) ptr.reset(new std::pair<int,int>(qstart, qend));
		//only left spanning
		else if(msg.spanning_l) ptr.reset(new std::pair<int,int>(qstart, projectedcoords.size()));
		//only right spanning
		else if(msg.spanning_r) ptr.reset(new std::pair<int,int>(0, qend));
		//neither spanning, take whole read
		else ptr.reset(new std::pair<int,int>(0, projectedcoords.size()));
	}
}

void abg_generate(bam1_t* read,	int rstart, int rend, abg_generate_msg& msg, std::string& seq)
{
	std::unique_ptr<std::pair<int,int>> query_ptr;
	std::vector<int> projectedcoords;
	bool clipped_l = false, clipped_r = false;
	project_positions(read, clipped_l, clipped_r, projectedcoords);
	get_breakpoints(rstart, rend, clipped_l, clipped_r, msg, projectedcoords, query_ptr);
	if(msg.successful){
		if((query_ptr->first == -1 && query_ptr->second != -1) || (query_ptr->first != -1 && query_ptr->second == -1)){
			std::cerr << "ERROR: unexpected querty start/end coords found for read " << (char*)read->data << '\n';
			exit(1);
		}
		if(query_ptr->first == -1) seq = ""; 
		else {
			int l_qsubseq = query_ptr->second - query_ptr->first;
			//quality string
			uint8_t *q = bam_get_seq(read);
			//update array for read seqeunce
			for(int i = 0; i  < l_qsubseq; i++) seq += seq_nt16_str[bam_seqi(q, i + query_ptr->first)];
			//seq = read_seq.substr(query_ptr->first, query_ptr->second - query_ptr->first);
		}
	}	
}

void bam_iterator(ansbam& bam_inst)
{
	abg_block current_block;
  	std::vector<abg_block> loaded_blocks;
}

void write2bam(samFile* outfile, const std::string& region_bed, const std::vector<uint32_t>& indeces, const std::vector<abg>& abg_reads, const std::vector<std::string>& edge_seqs)
{
	for(const auto& i : indeces){
		if(i >= edge_seqs.size() || i >= abg_reads.size()){
			std::cerr << "Error at region " << region_bed << ": index " << i << " out of bounds " << abg_reads.size() << ',' << edge_seqs.size() << '\n';
			exit(1);
		}
		const std::string& rseq = edge_seqs[i];
		const abg& abg_read = abg_reads[i];
		//std::cerr << abg_read.name << '\n';
		bam1_t *q = bam_init1();
	    //`q->data` structure: qname-cigar-seq-qual-aux
	    q->l_data = abg_read.name.size() + 1 + (int)(1.5*rseq.size() + (rseq.size() % 2 != 0));
	    q->m_data = q->l_data;
	    //std::cerr << q->l_data << ',' << q->m_data << '\n';
	    q->data = (uint8_t*)realloc(q->data, q->m_data);
	    //std::cerr << "realloc\n";
	    q->core.l_qname = abg_read.name.size() + 1; // +1 includes the tailing '\0'
	    q->core.l_qseq = rseq.size();
	    q->core.n_cigar = 0; // we have no cigar sequence
	    q->core.tid = -1;
	    q->core.pos = -1;
	    q->core.mtid = -1;
	    q->core.mpos = -1;
	    memcpy(q->data, abg_read.name.c_str(), q->core.l_qname); // first set qname
	    //std::cerr << "memcpy\n";
	    uint8_t *s = bam_get_seq(q);
	    for (int i = 0; i < q->core.l_qseq; ++i) bam1_seq_seti(s, i, seq_nt16_table[rseq[i]]);
	    //std::cerr << "seq set\n";
	    s = bam_get_qual(q);
	    for (int i = 0; i < q->core.l_qseq; ++i) s[i] = base_qual;
	    //std::cerr << "qual set\n";
	    std::string tag;
		abg_read.to_bam_aux(region_bed, tag);
	    abg_bam_aux_update(region_bed, abg_read, q);
	    //std::cerr << "aux updated\n";
	    if(sam_write1(outfile, nullptr, q) == -1) std::cerr << "Failed writing at " << abg_read.name << " at region " << region_bed << '\n';
	    bam_destroy1(q);
	    //std::cerr << "bam destroyed\n";
	}
}

abg_block_iter::abg_block_iter(const std::string& file, const uint32_t& b, const uint32_t& m, const bool& h): block_size(b), maxcov(m), hp(h){bam_inst.init(file, false);};

bool abg_block_iter::next()
{
	loaded_blocks.clear();
	while(loaded_blocks.size() < block_size){
		if(sam_read1(bam_inst.fp, bam_inst.header, bam_inst.read) <= 0){
			if(!current_block.empty() && current_block.size() <= maxcov) {
				loaded_blocks.emplace_back(current_block);
				current_block.clear();
			}
			break;
		}
		else {
			std::string local_name; //region
			abg abg_read(bam_inst.read, local_name);
			if(current_block.name.empty()) current_block.name = local_name;
			std::string local_seq = "";
			//quality string
			uint8_t *q = bam_get_seq(bam_inst.read);
			//update array for read seqeunce
			for(int i=0; i < bam_inst.read->core.l_qseq ; i++) local_seq += seq_nt16_str[bam_seqi(q,i)];
			if(abg_read.realigned) local_seq = local_seq.substr(abg_read.realigned->first, abg_read.realigned->second);
			//compress if user provided
			if(hp && local_seq.size() > 0) {
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
				if(current_block.size() <= maxcov) loaded_blocks.emplace_back(current_block);
				current_block.clear();
				current_block.name = local_name;
				current_block.seqs.emplace_back(local_seq);
				current_block.reads.emplace_back(abg_read);
			}
		}
  }
  return !loaded_blocks.empty();
};

void abg_block_iter::destroy(){bam_inst.destroy();};