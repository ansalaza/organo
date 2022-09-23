#ifndef ANDISTMATRIX_HPP
#define ANDISTMATRIX_HPP

#include <cstdint>
#include <vector>

class andistmatrix{
	public:
		std::vector<double> distmat;
		andistmatrix(uint32_t);
		void set_dist(uint32_t, uint32_t, double);
		double get_dist(uint32_t, uint32_t) const;
		uint32_t representative(std::vector<uint32_t>&) const;
		uint32_t size()const;

	private:
		uint32_t n;
		uint32_t d_n;
};

#endif
