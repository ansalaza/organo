#ifndef PAIRWISE_ALIGNMENTS_HPP
#define PAIRWISE_ALIGNMENTS_HPP

#include <string>
#include <vector>
#include "organo_opts.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "abg.hpp"
#include "andistmatrix.hpp"

void spanning_alignment(
	wfa::WFAligner& aligner, 
	std::string& query, 
	abg& q_tag, 
	std::string& ref
	);

void trimming_pairwise_alignment(
	wfa::WFAligner& aligner,
	const organo_opts& params, 
	std::vector<abg>& abg_reads,
	std::vector<std::string>& seqs, 
	andistmatrix& distmatrix
	);

void realignment2(
	wfa::WFAligner& aligner_gap, 
	const organo_opts& params, 
	std::vector<abg>& abg_reads, 
	std::vector<std::string>& seqs
	);

void realignment(
	wfa::WFAligner& aligner_edit, 
	wfa::WFAligner& aligner_gap, 
	const organo_opts& params, 
	std::vector<abg>& abg_reads, 
	std::vector<std::string>& seqs
	);

void spanning_aware_pairwise_alignment(
	wfa::WFAligner& aligner, 
	const organo_opts& params, 
	std::vector<abg>& abg_reads, 
	std::vector<std::string>& seqs, 
	std::vector<uint32_t>& seqs_l,
	andistmatrix& distmatrix
	);

void nonspanning_assignment_task(
	const organo_opts& params, 
	andistmatrix& distmatrix, 
	const std::vector<abg>& sequences,
	const std::vector<std::vector<uint32_t>>& ccs, 
	std::vector<std::vector<uint32_t>>& ccs_expanded,
	std::vector<uint32_t>& unassigned
	);

#endif