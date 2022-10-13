#ifndef ORGANO_OPTS_HPP
#define ORGANO_OPTS_HPP

#include <string>
#include <cstdint>

struct organo_opts {
	//abg_generate args
	uint8_t min_mapq;
	bool assembly;

	//trim args
	double maxerror;
	double minsim;
	int k;

	//merge args
	std::string matrix;
	std::string tsv;
	int maxalleles;
	int mincov;
	int lowcov;
	int maxcov;
	double minaf;

	//general args
	uint32_t threads;
	uint32_t block_size;
	bool group;
	bool haplotype;
	bool hp;
	bool hpd;

	organo_opts(){};
};

#endif