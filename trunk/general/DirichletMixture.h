#ifndef DIRICHLETMIXTURE_H
#define DIRICHLETMIXTURE_H

#include <math.h>
#include <string>
#include <cfloat>
#include "param.h"

using namespace std;

namespace sciphy
{

  class DirichletMixture
  {
  public:
    DirichletMixture(string mixture);
    ~DirichletMixture();
    
    void logProbOfCounts(float* c, int j, double& val);
    void calcPosteriorProb(float* c, double* mix); 
    void calcPosteriorProb_Old(float* c, double* mix);
    void getDirichletPriors(float* priors);
    void logProbOrderedCounts(int* c, double& val);

  private:
    string name;
    int numComp;
    float* coeff;

    /* alpha values are in comp[c][1..alphabet_size] 
       comp[c][0] is the sum of alphas for the component
    */
    float** comp;

    /* we precompute these values which will be reused over and over */
    double* log_coeff;
    double** lgamma_comp;

    void initialize();
    void initLogGammaValues();
  };

}

#endif
  
