
#ifndef THERMR_LEGNDR_HH
#define THERMR_LEGNDR_HH
 
#include <boost/math/special_functions/legendre.hpp>

inline auto legndr( double x, std::vector<double>& p, int np ){
 /*--------------------------------------------------------------------
  * Generate Legendre polynomials at x by recursion.
  * Place p(subl) in p(l+1).
  *--------------------------------------------------------------------
  */
  
  // throw exception if x is not in valid range
  if ( x < -1 or x > 1 ){ throw std::exception(); }


  if (np > 0){
    // If you need the 2nd or higher order, must actually calculate things
    for ( size_t i = 0; int(i) < np + 1 and i < p.size(); ++i ){
      p[i] = boost::math::legendre_p(i,x);
    }
  }
  else {
    p[0] = 1; p[1] = x;
  }

}

#endif
