
auto do155( int& k, int nw, double tsq, double f, std::vector<double>& wrk ){
  //155 continue
  std::cout << "155 " << k  << std::endl;
  k=k+1;
  if ((2*k) > nw) std::cout << "call error('sigcoh','storage exceeded.',' ')" << std::endl;
  wrk[2*k-1-1]=tsq;
  wrk[2*k-1]=f;
}

