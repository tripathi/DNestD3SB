#ifndef DNest4_Data
#define DNest4_Data

#include <valarray>

/*
* An object of this class is a dataset
*/
class Data
{
	private:
		// The visibility data
		std::valarray<double> rho;
		std::valarray<double> dvis;
		std::valarray<double> dwgt;

		// Bins and mean SB from image
		std::valarray<double> binedge; //Upper limit of each bin
		std::valarray<double> rsb;
		std::valarray<double> sbmean;

	public:
		// Constructor
		Data();

		// Load data from a file
		void loadvis(const char* filename);
		void loadsb(const char* filename);

		// Access to the data points
		const std::valarray<double>& get_rho() const
		{ return rho; }
		const std::valarray<double>& get_dvis() const
		{ return dvis; }
		const std::valarray<double>& get_dwgt() const
		{ return dwgt; }

		const std::valarray<double>& get_binedge() const
		{ return binedge; }
		const std::valarray<double>& get_rsb() const
		{ return rsb; }
		const std::valarray<double>& get_sbmean() const
		{ return sbmean; }

	private:
		// Static "global" instance
		static Data instance;

	public:
		// Getter for the global instance
		static Data& get_instance()
		{ return instance; }
};

#endif
