
auto legndr( const double& x, std::vector<double>& p ){
  /*--------------------------------------------------------------------
   * Generate Legendre polynomials at x by recursion.
   * Place p(subl) in p(l+1).
   *--------------------------------------------------------------------
   */
   int m1,i;
   double g,h;

   p[0]=1;
   p[1]=x;
   if (p.size() < 2) return;
   m1 = p.size() - 1;
   for ( int i = 0; i < m1; ++i ){
      g=x*p[i+1];
      h=g-p[i];
      p[i+2]=h+g-h/(i+1+1);
   }
   return;
}


