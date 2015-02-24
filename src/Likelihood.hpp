#pragma once
#include <map>
#include <vector>


#include "LikelihoodParams.hpp"
#include "Peptide.hpp"

class Likelihood {
	typedef std::is_same<std::vector<Peptide>::size_type pepupdate_st;

	std::vector<Peptide> peptides;
	std::vector<double> likelihoods;
	static_assert(pepupdate_st, likelihoods::size_type>, "Vector of likelihood has a different size_type from the Peptides")

	double likelihood;
	// mapping from param to peptides that must be updated
	// Maps an index on the c to a vector of indices on peptides to be updated
	std::map<pepupdate_st, std::vector<pepupdate_st>> onupdate_c;
	// Maps an index on the o to a vector of indices on peptides to be updated
	std::map<std::string, std::map<pepupdate_st, std::vector<pepupdate_st>>> onupdate_o;

}