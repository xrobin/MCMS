#pragma once

#include <Rcpp.h>
#include <string>

/** This is a basic wrapper around Rcpp::S4 that only ensures we have a "Peptides" or "Protein" object
 * It contains a single accessor method: slot, that behaves like that in the S4 class.
 * PeptidesModel has a getProteinModel slot as well.
 */


class ProteinModel  {
	public:
	Rcpp::S4 model;

	ProteinModel(const Rcpp::S4 aModel): model(aModel) {
		const std::string modelClass = std::string(Rcpp::as<Rcpp::CharacterVector>(aModel.attr("class"))[0]);
		if (modelClass != "Protein") {
			Rcpp::stop(std::string("Model of class Protein expected, ") + modelClass + " received");
		}
	}

	SEXP slot(const std::string& name) const {
		return model.slot(name);
	}
};


class PeptidesModel  {
	public:
	Rcpp::S4 model;

	PeptidesModel(const Rcpp::S4 aModel): model(aModel) {
		const std::string modelClass = std::string(Rcpp::as<Rcpp::CharacterVector>(aModel.attr("class"))[0]);
		if (modelClass != "Peptides") {
			Rcpp::stop(std::string("Model of class Peptides expected, ") + modelClass + " received");
		}
	}

	SEXP slot(const std::string& name) const {
		return model.slot(name);
	}

	ProteinModel getProteinModel() const {
		return model.slot("protein");
	}
};