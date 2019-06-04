#ifndef THERMR_SIMPLER_ELASTIC_COHERENT_SIGCOH_LEGNDR 
#define THERMR_SIMPLER_ELASTIC_COHERENT_SIGCOH_LEGNDR 

inline auto legndr( double x, std::vector<double>& p, int np ){
 /*--------------------------------------------------------------------
  * Generate Legendre polynomials at x by recursion.
  * Place p(subl) in p(l+1).
  *--------------------------------------------------------------------*/

  // throw exception if x is not in valid range
  if ( x < -1 or x > 1 ){ throw std::exception(); }
  p[0] = 1; p[1] = x;
  for (int i = 0; i < np-1; ++i){
    p[i+2] = x*p[i+1]-p[i]+x*p[i+1]-(x*p[i+1]-p[i])/(i+2);
  }
}

#endif
