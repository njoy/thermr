
#ifndef THERMR_LEGNDR_HH
#define THERMR_LEGNDR_HH


template <typename Range, typename Float>
inline auto legndr( const Float& x, Range& p, int np ){
 /*--------------------------------------------------------------------
  * Generate Legendre polynomials at x by recursion.
  * Place p(subl) in p(l+1).
  *--------------------------------------------------------------------
  */
  //if ( x < -1 or x > 1 ){ throw std::exception(); }
  // Defined to be the first two Legendre polynomial values.
  p[0] = 1; p[1] = x;
  if (np > 1){
    // If you need the 2nd or higher order, must actually calculate things
    for ( int i = 0; i < np - 1 and i < int(p.size())-2; ++i ){
      p[i+2] = 2*x*p[i+1] - p[i] - (x*p[i+1]-p[i]) / (i+2);
    }
  }
}
 
#endif
