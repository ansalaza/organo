#include <string>
#include <vector>
#include <memory>
#include <htslib/sam.h>
#include <iostream>
#include <sstream>
#include "abg.hpp"

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2))
std::string NAME_TAG = "R:Z:";
std::string REGION_TAG = "TR:Z:";
std::string SPANNING_TAG = "SP:Z:";
std::string SUPPLEMENT_TAG = "SA:i:";
std::string REALIGNED_TAG = "RA:Z:";
std::string HAPLOTYPE_TAG = "HP:i:";
std::string HIFI_PASSES_TAG = "np:i:";
std::string READ_DEPTH_TAG = "DP:i:";
std::string ALLELE_FREQ_TAG = "AF:f:";
std::string CC_ID_TAG= "CC:i:";

abg::abg()
{
	spanning_l = true;
	spanning_r = true;
	supplement = false;
};
abg::abg(bam1_t* b, std::string& region)
{
	spanning_l = true;
	spanning_r = true;
	supplement = false;

	name = (char*)b->data;
	std::string tag = REGION_TAG + ((char*)bam_get_aux(b) + 3);
	auto aux_ptr = bam_aux_get(b, SPANNING_TAG.substr(0,2).c_str());
	if(aux_ptr != NULL) tag += ' ' + SPANNING_TAG + (char*)(aux_ptr+1);
	aux_ptr = bam_aux_get(b, SUPPLEMENT_TAG.substr(0,2).c_str());
	if(aux_ptr != NULL) tag += ' ' + SUPPLEMENT_TAG + std::to_string(bam_aux2i(aux_ptr));
	aux_ptr = bam_aux_get(b, REALIGNED_TAG.substr(0,2).c_str());
	if(aux_ptr != NULL) tag += ' ' + REALIGNED_TAG + (char*)(aux_ptr+1);
	aux_ptr = bam_aux_get(b, HAPLOTYPE_TAG.substr(0,2).c_str());
	if(aux_ptr != NULL) tag += ' ' + HAPLOTYPE_TAG + std::to_string(bam_aux2i(aux_ptr));
	aux_ptr = bam_aux_get(b, HIFI_PASSES_TAG.substr(0,2).c_str());
	if(aux_ptr != NULL) tag += ' ' + HIFI_PASSES_TAG + std::to_string(bam_aux2i(aux_ptr));
	if(!tag.empty()) parse(tag, region);
};


void abg::parse(std::string& str, std::string& region)
{
	std::string value;
    std::istringstream stream(str);
    while(std::getline(stream, value, ' ')) {
    	//supplement alignment tag
    	if(value.substr(0,5) == SUPPLEMENT_TAG) supplement = value.substr(5) == "1" ? true : false;
    	//cargano-specific realignment tag
    	else if(value.substr(0,5) == REALIGNED_TAG) {
    		std::string value2;
    		std::istringstream stream2(value.substr(5));
    		int l = -1;
    		int r = -1;
    		int i = 0;
    		while(std::getline(stream2, value2, ',')) {
    			if(i == 0) l = std::stoi(value2);
    			else if(i == 1) r = std::stoi(value2);
    			++i;
    		}
    		if(l == -1 || r == -1) {
    			std::cerr << "Unexpected 'RA:Z:' tag:" << value << '\n';
    			exit(1);
    		}
    		else realigned.reset(new std::pair<int,int>(l,r));
    	}
    	//target region tag
    	else if(value.substr(0,5) == REGION_TAG) region = value.substr(5);
    	//haplotype tag
    	else if(value.substr(0,5) == HAPLOTYPE_TAG) haplotype = std::make_shared<int>(std::stoi(value.substr(5)));
    	//hifi number of passes tag
    	else if(value.substr(0,5) == HIFI_PASSES_TAG) np = std::make_shared<int>(std::stoi(value.substr(5)));
    	//spanning tag
    	else if(value.substr(0,5) == SPANNING_TAG){
    		char info = value.substr(5)[0];
    		if(info == 'n'){
    			spanning_l = false;
    			spanning_r = false;
    		}
    		else if(info == 'l'){
    			spanning_l = true;
    			spanning_r = false;
    		}
    		else if(info == 'r'){
    			spanning_l = false;
    			spanning_r = true;
    		}
    		else if(info == 'b'){
    			spanning_l = true;
    			spanning_r = true;
    		}
    	}
    }
}

bool abg::spanning() const {return spanning_l && spanning_r;};

void abg::to_string(std::string& tmp) const
{
	tmp = NAME_TAG + name;

	tmp += ' ' + SPANNING_TAG;
	if(spanning_l && spanning_r) tmp += 'b';
	else if(spanning_l && !spanning_r) tmp += 'l';
	else if(!spanning_l && spanning_r) tmp += 'r';
	else tmp += 'n';

	tmp += ' ' + SUPPLEMENT_TAG;
	tmp += supplement ? "1" : "0";

	if(realigned){
		tmp += ' ' + REALIGNED_TAG;
		tmp += std::to_string(realigned->first) + ',' + std::to_string(realigned->second);
	}

	if(haplotype) tmp += ' ' + HAPLOTYPE_TAG + std::to_string(*haplotype);
	if(np) tmp += ' ' + HIFI_PASSES_TAG + std::to_string(*np);	
}

void abg::to_bam_aux(const std::string& region, std::string& tmp) const
{
	tmp = region;
	tmp += ' ' + SPANNING_TAG;
	if(spanning_l && spanning_r) tmp += 'b';
	else if(spanning_l && !spanning_r) tmp += 'l';
	else if(!spanning_l && spanning_r) tmp += 'r';
	else tmp += 'n';

	tmp += ' ' + SUPPLEMENT_TAG;
	tmp += supplement ? "1" : "0";

	if(realigned){
		tmp += ' ' + REALIGNED_TAG;
		tmp += std::to_string(realigned->first) + ',' + std::to_string(realigned->second);
	}

	if(haplotype) tmp += ' ' + HAPLOTYPE_TAG + std::to_string(*haplotype);
	if(np) tmp += ' ' + HIFI_PASSES_TAG + std::to_string(*np);
}

void abg_bam_aux_update(const std::string& region, const abg& tag, bam1_t *read)
{

	bam_aux_update_str(read, REGION_TAG.substr(0,3).c_str(), region.size() + 1, region.c_str());
	std::string tmp_tag;
	if(tag.spanning()) tmp_tag += 'b';
	else if(tag.spanning_l) tmp_tag += 'l';
	else if(tag.spanning_r) tmp_tag += 'r';
	else tmp_tag += 'n';
	int supplement_tag = tag.supplement ? 1 : 0;
	bam_aux_update_str(read, SPANNING_TAG.substr(0,2).c_str(), tmp_tag.size() + 1, tmp_tag.c_str());
	bam_aux_update_int(read, SUPPLEMENT_TAG.substr(0,2).c_str(), supplement_tag);
	if(tag.realigned){
		tmp_tag = std::to_string(tag.realigned->first) + ',' + std::to_string(tag.realigned->second);
		bam_aux_update_str(read, REALIGNED_TAG.substr(0,2).c_str(), tmp_tag.size() + 1, tmp_tag.c_str());
	}
	if(tag.haplotype) bam_aux_update_int(read, HAPLOTYPE_TAG.substr(0,2).c_str(), *tag.haplotype);
	if(tag.np) bam_aux_update_int(read, HIFI_PASSES_TAG.substr(0,2).c_str(), *tag.np);
}


abg_block::abg_block(){};
std::size_t abg_block::size() const {return reads.size();}
bool abg_block::empty() const {return reads.empty();}
void abg_block::clear(){name.clear(); reads.clear(); seqs.clear();}