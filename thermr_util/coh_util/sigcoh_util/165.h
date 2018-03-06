
auto do165( int k, std::vector<double>& wrk, double recon, double ulim, double f ){
  for ( int i = 1; i < k; ++i ){
   //if (tsq.lt.wrk(2*i-1).or.tsq.ge.(1+eps)*wrk(2*i-1)) go to 170
    wrk[2*i-1]=wrk[2*i-1]+f;
    return end_175_180_185( wrk, k, recon, ulim );
  }
  // this really shouldn't be here
  return wrk;
}


