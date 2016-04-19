#ifndef DNest4_Template_MyModel
#define DNest4_Template_MyModel

#include "DNest4/code/DNest4.h"
#include <valarray>
#include <ostream>

class MyModel
{
	private:
		//The bin weights
		std::valarray<double> params;

		//Model visibilities
		std::valarray<double> mvis;
		//const std::valarray<double>& sbmean;

		//Number of bins below resolution
		int belowres;

		//Compute the model visibilities given the bin weights
		void calculate_mvis();

	public:
		// Constructor only gives size of params
		MyModel();

		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
};

#endif
