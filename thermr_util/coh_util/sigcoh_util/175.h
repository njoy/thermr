#include <iostream> 
#include <vector>


auto end_175_180_185(std::vector<double>& wrk, int k, double recon, double ulim ){
   int imax=k-1;
   for ( int i = 0; i < imax; ++i ){
      int jmin=i+1;
      for ( int j = jmin; j < k; ++j ){
         if ( wrk[2*j-2] < wrk[2*i-2] ){
            double st=wrk[2*i-1-1];
            double sf=wrk[2*i-1];
            wrk[2*i-1-1]=wrk[2*j-1-1];
            wrk[2*i-1]=wrk[2*j-1];
            wrk[2*j-1-1]=st;
            wrk[2*j-1]=sf;
          }
        }
   }
   k=k+1;
   wrk[2*k-1-1]=ulim;
   wrk[2*k-1]=wrk[2*k-2-1];
   int nw=2*k;
   double enext=recon*wrk[0];
   //enext=sigfig(enext,7,-1)
   // copy data to global fl array
   //for ( int i = 0; i < 15; ++i ){ std::cout << wrk[i] << std::endl; }
   std::vector<double> fl(nw,0.0);
   for ( int i = 0; i < nw; ++i ){
      fl[i]=wrk[i];
    }
   int nl=k;
   for ( int i = 0; i < 18; ++i ){ std::cout << fl[i] << std::endl; }
   return fl;
}


