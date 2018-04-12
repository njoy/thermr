
auto loada( int& i, int& na, std::fstream& ntape, int& nbuf, std::vector<double>& a, 
    std::vector<double>& buf ){
  /*--------------------------------------------------------------------
   * Buffered sequential i/o routine.
   * Store na elements of array a into
   * core buffer and associated binary tape.
   *--------------------------------------------------------------------
   */
   int nl,j,ix,inow,k;
   double x,xnow;

   nl=nbuf/na;
   j=abs(i);
   x=j;
   x=x/nl;
   ix=j/nl;
   xnow=(x-ix)*nl;
   inow=int(xnow);
   if (inow == 0) inow=nl;
   k=na*(inow-1);
   for ( int j = 0; j < na; ++j ){
      k=k+1;
      buf[k-1]=a[j];
   }
  // if (i == 1) call repoz(-ntape)
   std::string tab = "\t", newline = "\n";
  if (inow == nl or i < 0){
    for ( size_t i1 = 0; i1 < buf.size(); ++i1 ){ 
      ntape >> buf[i1] >> tab;
      if ( i1%3 == 0 ){ ntape >> newline; }
    }
  }

}
