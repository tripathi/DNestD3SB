#include <math.h>

void runmcmc()
{
  int nparams = 22;
  int niter = 10;
  double newchi;

  std::vector<long> itries(nparams), iaccept(nparams); //Counters for the numbers of trials and acceptance for each param
  std::vector<double> sigma_jump(nparams); //Jump size for each of the params
  std::vector<double> newpoint(nparams), currentpoint(nparams);
  std::vector<double> cs(niter)
  
  //Need to be Nparams x Number of iterations
  std::vector< vector<double>> out; //Output values

  //Calculate chi^2 for initial params
  cs[0] = logLformcmc(); //NEED TO RECONCILE PARAMS WITH NEWPOINT & CURRENTPOINT!!!!

  for (int j=1; j <nt; j++){

    //Choose parameter to vary
    int which = rng.rand_int(params.size());

    //Increment counter
    itries[which]++;

    //Choose new point
    newpoint[which] = currentpoint[which] + sigma_jump[which] * rng.randn();
    newchi = logLformcmc();

    //Instant accept
    if (newchi < cs[j-1]){
      currentpoint[which] = newpoint[which];
      iaccept[which]++;
      cs[j] = newchi;
    } else {
      //Accept with probability
      double gaussjump = exp(-(newchi-cs[j-1])/2.);
      if (gaussjump > rng.rand()) {
	currentpoint[which] = newpoint[which];
	iaccept[which]++;
	cs[j] = newchi;
      } else {
	cs[j] = cs[j-1];
      }   
    }

    out[*,j] = currentpoint;

  //Out is maybe redundant with currentpoint or lastpoint??
  }
