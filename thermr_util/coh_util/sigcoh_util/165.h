
auto do165( int k, std::vector<double>& wrk, double recon, double ulim, double f, double tsq, double eps ){
  std::cout << "165" << std::endl;
  for ( int i = 1; i <= k; ++i ){
    std::cout <<"--------- " <<  i << std::endl;
    if (tsq < wrk[2*i-1-1] or tsq >= (1+eps)*wrk[2*i-1-1]) continue; //go to 170
    wrk[2*i-1]=wrk[2*i-1]+f;
    return end_175_180_185( wrk, k, recon, ulim );
  }
  // this really shouldn't be here
  return end_175_180_185( wrk, k, recon, ulim );
}


