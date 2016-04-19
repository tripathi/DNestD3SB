#include "Data.h"
#include <fstream>
#include <vector>

using namespace std;

// The static instance
Data Data::instance;

Data::Data()
{

}

void Data::loadvis(const char* filename)
{
	// Vectors to hold the data
	std::vector<double> _rho;
	std::vector<double> _dvis;
	std::vector<double> _dwgt;

	// Open the file
	fstream fin(filename, ios::in);

	// Temporary vector variables
	double tempu, tempv, tempvis, tempwgt ;
	double arcsec = 180./M_PI*3600.;

	// Read until end of file
	while(fin>>tempu && fin>>tempv && fin>>tempvis && fin>>tempwgt)
	{
	  _rho.push_back(sqrt(tempu*tempu + tempv*tempv)/arcsec);
	  _dvis.push_back(tempvis);
	  _dwgt.push_back(tempwgt);
	}

	// Close the file
	fin.close();

	// Copy the data to the valarrays
	rho = valarray<double>(&_rho[0], _rho.size());
	dvis = valarray<double>(&_dvis[0], _dvis.size());
	dwgt = valarray<double>(&_dwgt[0], _dwgt.size());
}

void Data::loadsb(const char* filename)
{
	// Vectors to hold the data
	std::vector<double> _rsb;
	std::vector<double> _binedge={0.01/140.}; //Placeholder for rin
	std::vector<double> _sbmean;

	// Temporary vector variables
	double tempb, temprsb, tempsbmean ;

	// Open the file
	fstream fin(filename, ios::in);

	// Read until end of file
	while(fin>>tempb && fin>>temprsb && fin>>tempsbmean)
	{
		_binedge.push_back(tempb);
		_rsb.push_back(temprsb);
	  _sbmean.push_back(tempsbmean);
	}

	// Close the file
	fin.close();

	// Copy the data to the valarrays
	binedge = valarray<double>(&_binedge[0], _binedge.size());
	rsb = valarray<double>(&_rsb[0], _rsb.size());
	sbmean = valarray<double>(&_sbmean[0], _sbmean.size());
}
