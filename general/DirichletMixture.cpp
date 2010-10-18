#include "DirichletMixture.h"
#include "blocks9.h"
#include "recode4.h"
#include <iostream>
#include <vector>

using namespace std;

namespace sciphy{
  
#define FAST_LGAMMA

  //int n_reuse = 0, n_new = 0;

#ifdef FAST_LGAMMA

#define LG_RESOLUTION 0.001
#define LG_MULTIPLY (int) (1/LG_RESOLUTION)
#define LG_MAX_NUM 1000000
#define LG_MAX (int) (LG_MAX_NUM*LG_RESOLUTION)

  int maxCachedValue = -1;

  static float* lgammaCache = NULL;
  
  static float my_lgamma_compute(float v);

  static inline float
  my_lgamma(float v) 
  {
    int k = (int) (v*LG_MULTIPLY);

    if (v >= LG_MAX) {
      return lgamma(v);
    } else if (v < maxCachedValue && lgammaCache[k] < FLT_MAX) {
      //n_reuse++;
      return lgammaCache[k];
    } else {
      //n_new++;
      return my_lgamma_compute(v);
    }
  }


  static float
  my_lgamma_compute(float v)
  {
    int k = (int) (v*LG_MULTIPLY);
    float res = lgamma(v);

    /* grow the cache as necesary */
    if (v >= maxCachedValue) {

      if (!lgammaCache) {
	maxCachedValue = 10;
	int size = maxCachedValue * LG_MULTIPLY;
	lgammaCache = (float *) malloc( sizeof(float) * size);
	if (!lgammaCache) {
	  cerr << "Out of memory error: failed to allocate memory to lgammaCache" << endl;
	  exit(1);
	}
	for (int i=0; i<size; i++) 
	  lgammaCache[i] = FLT_MAX;
	
      } else {
	int begin = maxCachedValue * LG_MULTIPLY;
	maxCachedValue = (int) (v + 1);
	int size = maxCachedValue * LG_MULTIPLY;
	lgammaCache = (float *) realloc(lgammaCache, sizeof (float) * size);
	if (!lgammaCache) {
	  cerr << "Out of memory error: failed to re-allocate memory to lgammaCache" << endl;
	  exit(1);
	}
	for (int i=begin; i<size; i++) 
	  lgammaCache[i] = FLT_MAX;
	
	//cout << "resized cache size to " << maxCachedValue << endl;
      }
    }

    if (lgammaCache[k] == FLT_MAX)
      lgammaCache[k] = res;

    return res;
  }

#else

#define my_lgamma(v) lgamma(v)

#endif

  DirichletMixture::DirichletMixture(string mixtureName)
  {
    name = mixtureName;
    initialize();
  }

  DirichletMixture::~DirichletMixture()
  {
  }

  void
  DirichletMixture::initialize() {
    int i, j;

    if (name == "blocks9")
      numComp = BLOCKS9_NMIX;
    else if (name == "recode4")
      numComp = RECODE4_NMIX;

    coeff = new float[numComp];
    comp = new float*[numComp];

    if (name == "blocks9") {
      for (i=0; i<numComp; i++) {
	coeff[i] = blocks9_coeffs[i];
	comp[i] = new float[ALPHABET_SIZE+1];
	for (j=0; j<=ALPHABET_SIZE; j++)
	  comp[i][j] = blocks9_alphas[i][j];
      }
    } else if (name == "recode4") {

      for (i=0; i<numComp; i++) {
	coeff[i] = recode4_coeffs[i];
	comp[i] = new float[ALPHABET_SIZE+1];
	for (j=0; j<=ALPHABET_SIZE; j++)
	  comp[i][j] = recode4_alphas[i][j];
      }
    }
    
    log_coeff = new double[numComp];
    lgamma_comp = new double*[numComp];

    for (int i=0; i<numComp; i++) {
      log_coeff[i] = log(coeff[i]);

      lgamma_comp[i] = new double[ALPHABET_SIZE+1];
      for (int j=0; j<=ALPHABET_SIZE; j++)
	lgamma_comp[i][j] = lgamma(comp[i][j]);
    }
  }

  /*
   * Return the distribution of amino acids as given by the mixture
   * this should be just the simple background distribution
   * p = q1*p1 + q2*p2 + ...
   */
  void
  DirichletMixture::getDirichletPriors(float* priors) {
    int i, j;
    double tot = 0;
    for (i = 0; i < ALPHABET_SIZE; i++) {
      priors[i] = 0;
      for (j = 0; j < numComp; j++) {
	priors[i] += coeff[j] * comp[j][i+1];
      }
      tot += priors[i];
    }

    for (i = 0; i < ALPHABET_SIZE; i++)
      priors[i] /= tot;
  }

  /*
   * Likehood of vector count c given the jth component
   *           G(n+1)*G(a^j)               G(n_i+a^j_i)
   * P(n|j) = --------------  * Prod    ------------------
   *             G(n+a^j)       i=1..20  G(n_i+1)*G(a^j_i)
   *
   * G() is the Gamma function
   * n=sum{i=1..N} n_i and a^j =sum{i=1..N} a^j_i
   */

  void
  DirichletMixture::logProbOfCounts(float* c, int j, double& val) 
  {
    int i;
    float n=0;

    for (i=0;i<ALPHABET_SIZE;i++) {
      n += c[i];
    }

    val = my_lgamma(n+1) + lgamma_comp[j][0] - my_lgamma(n+comp[j][0]);

    for (i=0; i<ALPHABET_SIZE; i++) {
      if (c[i] > 0)
	val += my_lgamma(c[i]+comp[j][i+1]) - my_lgamma(c[i]+1) - lgamma_comp[j][i+1];
    }
  }

  /* This gives the estimate of the posterior probability of a profile column
   * from the observed count vector c as given in Eq. 15 in Combios 1996.
   *                                  n_i + a^j_i
   *  p_i = sum_{j=1..l} Prob(j|n) ------------------
   *                                sum(n) + sum(a^j)
   */

  void
  DirichletMixture::calcPosteriorProb(float* c, double* mix) {

    /*
     * calculate posterior prob of P(q|c) of jth component given count vector c
     *                   q_j * Prob(n|j)
     * Prob(j|n) = ----------------------------
     *              sum{i=1..K} q_i * Prob(n|i)
     */

    double p, totp=0, totm=0;
    int i, j;
    double x, sumv=0;
    double* likelihood = new double[numComp];
    for (j=0; j<numComp; j++) {
      logProbOfCounts(c, j, x);
      //cout << "logp " << x << " exp(x) " << exp(x) << endl;
      likelihood[j] = exp(x);
      sumv += coeff[j] * likelihood[j];
    }

    double sumc = 0;
    for (i=0; i<ALPHABET_SIZE; i++) 
      sumc += c[i];

    for (i=0; i<ALPHABET_SIZE; i++) {
      totp = 0;
      for (j=0; j<numComp; j++) {
	p = (coeff[j] * likelihood[j]) / sumv;
	totp += p * (c[i] + comp[j][i+1])/(sumc + comp[j][0]);
      }
      mix[i] = totp;
      totm += mix[i];
    }

    // normalize
    for (i=0; i<ALPHABET_SIZE; i++)
      mix[i] /= totm;
    
    delete [] likelihood;
  }

  /*
   * Probability of generating the vector count c by the mixture
   * 
   * by jth component
   *              G(a^j)                 G(n_i+a^j_i)
   * P(n|j) = --------------  * Prod   ---------------
   *             G(n+a^j)       i=1..20   G(a^j_i)
   *
   * G() is the Gamma function
   * n=sum{i=1..N} n_i and a^j =sum{i=1..N} a^j_i
   */

  void
  DirichletMixture::logProbOrderedCounts(int* c, double& val) 
  {
    int i, j;
    double n=0, normalizer=-DBL_MAX, p;

    for (i=0;i<ALPHABET_SIZE;i++) 
      n += c[i];

    // we might have overflow problem when the counts are very large
    double* logp = new double[numComp];

    for (j=0; j<numComp; j++) {
      p = lgamma_comp[j][0] - my_lgamma(n+comp[j][0]);

      for (i=0; i<ALPHABET_SIZE; i++) {
	// if c[i] == 0, the right hand side is 0
	if (c[i] > 0) 
	  p += my_lgamma(c[i]+comp[j][i+1]) - lgamma_comp[j][i+1]; 
      }
      
      logp[j] = log_coeff[j] + p; // log( coeff[j]*Prob(n) )

      if (normalizer < logp[j])
	normalizer = logp[j];
    }

    double x = 0;
    for (j=0; j<numComp; j++) {
      x += exp(logp[j] - normalizer); // P=sum(q1*p1 + q2*p2 + ...)
    }

    val = normalizer + log(x);

    delete [] logp;
  }

}
