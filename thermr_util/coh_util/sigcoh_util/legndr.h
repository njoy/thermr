
auto legndr( double x, std::vector<double>& p, int np ){
 /*--------------------------------------------------------------------
  * Generate Legendre polynomials at x by recursion.
  * Place p(subl) in p(l+1).
  *--------------------------------------------------------------------
  */

  p[0] = 1;
  p[1] = x;
   if (np < 2) return;
   for ( int i = 0; i < np - 1; ++i ){
      p[i+2] = x*p[i+1] - p[i] + x*p[i+1] - (x*p[i+1]-p[i]) / (i+2);
   }
   return;
}


