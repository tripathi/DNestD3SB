#include <iostream>
#include "DNest4/code/DNest4.h"
#include "MyModel.h"
#include "Data.h"

using namespace DNest4;

int main(int argc, char** argv)
{
	Data::get_instance().loadsb("DATA/fullA22_bcbsb.txt");
	Data::get_instance().loadvis("DATA/fullAvis.txt");

	DNest4::start<MyModel>(argc, argv);

	MyModel test;
	//const auto& rsb = Data::get_instance().get_rsb();
  //fprintf(stderr, "First! %zu \n", b.size());
	//RNG  q;
	//test.from_prior(q);
	// test.log_likelihood();

	return 0;
}
