#ifndef ABG_HPP
#define ABG_HPP

#include <string>
#include <vector>
#include <memory>
#include <htslib/sam.h>

class abg{
	public:
		std::string name;
		bool spanning_l;
		bool spanning_r;
		bool supplement;
		std::shared_ptr<int> np;
		std::shared_ptr<std::pair<int,int>> realigned;
		std::shared_ptr<int> haplotype;

		abg();
		abg(bam1_t*, std::string&);
		bool spanning() const;
		void parse(std::string&, std::string&);
		void to_string(std::string&) const;
		void to_bam_aux(const std::string&, std::string&) const;
};

void abg_bam_aux_update(
	const std::string& region, 
	const abg& abg_read, 
	bam1_t *read
	);

struct abg_block{
	std::string name;
	std::vector<abg> reads;
	std::vector<std::string> seqs;

	abg_block();
	std::size_t size() const;
	bool empty() const;
	void clear();
};

#endif