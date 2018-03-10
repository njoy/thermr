#include "155.h"

auto do150( int& k, double& tsq, std::vector<double>& wrk, double& eps, double& f, int& nw ){
  std::cout << "150  "  << k<< std::endl;
  for ( int i = 1; i <= k; ++i ){
    //do 155 i=1,k
    if (tsq < wrk[2*i-1-1] or tsq >= (1+eps)*wrk[2*i-1-1]) { 
      if ( i == k ){ 
        do155( k, nw, tsq, f, wrk );
        break;
      }
      else {
        continue; 
      }
    }
    wrk[2*i-1]=wrk[2*i-1]+f;
    //go to 160
    return;
  }
}

