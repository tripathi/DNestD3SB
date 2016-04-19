#include "MyModel.h"
#include "Data.h"
#include <math.h>
#include "DNest4/code/DNest4.h"
#include <iostream>
#include <algorithm>
#include <gsl/gsl_sf_bessel.h>
#include <cassert>
#include <cmath>


using namespace std;
using namespace DNest4;

MyModel::MyModel()
:params(Data::get_instance().get_rsb().size()),
mvis(Data::get_instance().get_rho().size()),
jinc(Data::get_instance().get_rho().size()*Data::get_instance().get_binedge().size())
//sbmean(Data::get_instance().get_sbmean())
{

}

void MyModel::calculate_jinc()
{
	const auto& rho = Data::get_instance().get_rho(); //Change to u,v as needed
	const auto& rbin = Data::get_instance().get_binedge(); //nbins+1 size
	int n = rbin.size();

	for(size_t j=0; j<rho.size(); j++)
	{
		for (size_t i=0; i<rbin.size(); i++)
		{
			double theta = 2.*M_PI*rbin[i]*rho[j];
			jinc[j*n+i] = gsl_sf_bessel_J1 (theta) / theta; //jinc[i][j]
		}
	}
}
void MyModel::calculate_mvis()
{
	const auto& rho = Data::get_instance().get_rho(); //Change to u,v as needed
	const auto& rbin = Data::get_instance().get_binedge(); //nbins+1 size
	int n = rbin.size();

	//NOW in DATA.cpp - FIX! double rin = 0.001/140.; //Change to a variable
	//rbin[0] = rin;

//	mvis.reserve(rho.size()); //How do I incorporate this with .h init
//	valarray<double> jinc(rho.size());
	valarray<double> wdiff(rbin.size());

	//Need to calculate outer product, difference of weights, and BesselJ
	//double y = gsl_sf_bessel_J1 (x);

	// gsl_blas_dsdot (const gsl_vector_float * x, const gsl_vector_float * y, double * result)

	for (size_t i=0; i<params.size()-1; i++)
	{
		wdiff[i+1] = params[i]-params[i+1];
	}
	wdiff[0] = -params[0];
	wdiff[params.size()] = params[params.size()-1];

	 //theta = 2.*M_PI*rbin*rho;

	for(size_t j=0; j<rho.size(); j++)
	{
		mvis[j]=0;
		for (size_t i=0; i<rbin.size(); i++)
		{
			//fprintf(stderr, "i: %zu Wdiff is %e \n", i, wdiff[i]);
			// double theta = 2.*M_PI*rbin[i]*rho[j];
			// double jinc = gsl_sf_bessel_J1 (theta) / theta; //jinc[i][j]
			mvis[j] += 2.*M_PI*rbin[i]*rbin[i]*jinc[j*n+i]*wdiff[i]; //Need to be sum or dot prod
			//if (!isfinite(mvis[j]))	fprintf(stderr, "i: %zu j: %zu: mvis: %e \n", i, j, mvis[j]);
			//if (j ==22) fprintf(stderr, "%zu: wdiff: %f rbin: %e theta = %e jinc = %e \n", j, wdiff[i], rbin[i], theta, jinc);
		}
	}

	// for (size_t j=0; j<5; j++)
	// {
	// 	fprintf(stderr, " %zu: Rho = %e Vis = %e \n", j, rho[j], mvis[j]);
	// }

}

void MyModel::from_prior(RNG& rng)
{
	const auto& sbmean = Data::get_instance().get_sbmean();
	const auto& rho = Data::get_instance().get_rho(); //Change to u,v as needed
	const auto& rsb = Data::get_instance().get_rsb();
	double res = 1./rho.max();// THIS IS WRONG, buy using it for now *180./M_PI*3600.*3.0e8/340e9;
	fprintf(stderr, "WRONG Res is %e, max rho is %e \n", res, rho.max());
	belowres = 0;
	for (size_t i = 0; i<rsb.size(); i++)
	{
		if (rsb[i] <= res) belowres++;
	}
	fprintf(stderr,"Number below res is %d \n", belowres);
	//belowres = std::count_if(rsb.begin(), rsb.end(), [](int i) {return i < res;});

	for(size_t i=0; i<params.size(); i++)
	{
		if ((int) i>=belowres)
		{
			params[i] = sbmean[i]*10.*rng.rand();
		} else {
			params[i] = sbmean[i] *(1.+2.*rng.rand());
		}
		fprintf(stderr, "i:%zu SB:%f Param:%f Factor: %f \n", i, sbmean[i], params[i], params[i]/sbmean[i]);

		//Also need to manage the cases where below the resolution AND/OR add power law
	}
	calculate_jinc();
	calculate_mvis();
}

double MyModel::perturb(RNG& rng)
{
	const auto& sbmean = Data::get_instance().get_sbmean();
	int which = rng.rand_int(params.size()); //Which params to move
	params[which] += rng.randh();
	if (which >=belowres)
	{
		wrap(params[which], 0., 3.*sbmean[which]); //Does this need to be in a for loop?
	} else {
		wrap(params[which], sbmean[which], 3.*sbmean[which]);
	}

	//fprintf(stderr, "which: %d, param: %f \n", which, params[which]);
	calculate_mvis();

	double logH = 0.; //Leaving as 0 since using uniform priors
	return logH;
}

double MyModel::log_likelihood() const
{
//Enforce positive surface brightness
//(Do I still need to do this? Or is this taken care above?)

//Calculate penalty
	//std::vector<double> priori(params.size());
	int nturns = 0.; //Faster with adjacent_difference?
	double dw1, dw2;

	for (size_t i=1; i<params.size()-1; i++)
	{
		dw1 =params[i]-params[i-1];
		dw2 =params[i+1]-params[i];
		if(dw1*dw2 < 0) nturns++;
	}

	//	std::adjacent_difference(qq.begin(), qq.end(), priori.begin());
	//for(std::valarray<double>::iterator it=priori.begin(); it < priori.end(); ++it)
	//  fprintf(stderr, "%e", *it);

//Chi^2, compared to both component of vis and the weights


// Grab the values from the dataset
	const auto& dvis = Data::get_instance().get_dvis();
	const auto& dwgt = Data::get_instance().get_dwgt();

	for(size_t i=0; i<mvis.size(); i++)
	{
  	if (!isfinite(mvis[i])) fprintf(stderr,"Mvis: %e \n", mvis[i]);
		if (!isfinite(dvis[i])) fprintf(stderr,"Dvis: %e \n", dvis[i]);
		if (!isfinite(dwgt[i])) fprintf(stderr,"Dwgt: %e \n", dwgt[i]);
	}
	// assert(!std::isnan(mvis)));
	// fprintf(stderr,"Mvis isn't nan\n");
	// assert(!std::isnan(dvis)));
	// fprintf(stderr,"Dvis isn't nan\n");
	// assert(!std::isnan(dwgt)));
	double chi2 = (dwgt*pow(mvis-dvis,2.)).sum();
	double prior = nturns * 2.*dvis.size()/params.size();

	//fprintf(stderr, "Turns: %d Prior: %e Chi2: %e \n", nturns, prior, chi2);
	return -0.5*(chi2 +prior);
}

void MyModel::print(std::ostream& out) const
{
  for(size_t i=0; i<params.size(); i++)
  	out<<params[i]<<' ';
}

string MyModel::description() const
{
	return string("Bin weights in order from smallest to largest radius");
}
