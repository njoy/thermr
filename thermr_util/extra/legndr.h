
void legndr( const double& x, std::vector<double>& p ){
  /*--------------------------------------------------------------------
   * Generate Legendre polynomials at x by recursion.
   * Place p(subl) in p(l+1).
   *--------------------------------------------------------------------
   */

   p[0] = 1;
   p[1] = x;
   if (p.size() < 2) return;
   for ( int i = 0; i < p.size()-1; ++i ){
      p[i+2] = x*p[i+1] + (x*p[i+1] - p[i]) * ( 1 - 1.0 / (i+2) );
   }
}


