
auto do165( int& k, std::vector<double>& wrk, double& recon, double& ulim, double& f, double& tsq, double& eps, int nw ){
  std::cout << "165  "<<k << std::endl;
  for ( int i = 1; i <= k; ++i ){
    if (tsq < wrk[2*i-1-1] or tsq >= (1+eps)*wrk[2*i-1-1]){
      if ( i == k ){ 
//        std::cout << "170" << std::endl;
//        k = k + 1;
//        if ((2*k) > nw) std::cout << "call error('sigcoh','storage exceeded.',' ')" << std::endl;
//        wrk[2*k-1-1] = tsq;
//        wrk[2*k-1] = f;
        return false;
      }
      continue; //go to 170
    }
    wrk[2*i-1]=wrk[2*i-1]+f;
    return true;
  }

  return true;
}


