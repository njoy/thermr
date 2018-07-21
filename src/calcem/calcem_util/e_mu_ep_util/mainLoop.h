#include <iostream>

auto mainLoop(const int& nep, const int& npage, int& j, int& k, int& ib,
  std::vector<double>& scr, const std::vector<double>& yu, const double& sum, 
  int& nw, int& ncds, int& nb ){
  int iend, istart = 1;
  while (true) {
    std::cout << 595 << std::endl;
    iend=nep;
    if ((iend-istart) >= npage/2) iend=istart+npage/2-1;
    j=k-1;
    ib=istart-1;

    do {
     ///std::cout << 596 << std::endl;
     j=j+2;
     ib=ib+1;
     scr[j-1]=yu[1+2*ib-1];
     scr[j+1-1]=yu[2+2*ib-1]*2/sum;
    }
    while ( ib < iend );
    nw=j+1;
    if (k == 0) {
      std::cout << 597 << std::endl;
      //call moreio(0,0,nscr,scr,nb,nw)
      if (nb == 0) {
        std::cout << 598 << std::endl;
        ncds=ncds+1+(j*(nep+1)+5)/6;
        return std::make_tuple(istart,iend);
      }
      istart=iend+1;
      continue;
    }
    k=0;
    //call tab1io(0,0,nscr,scr,nb,nw)
    if (nb == 0) {
      std::cout << 598 << std::endl;
      ncds=ncds+1+(j*(nep+1)+5)/6;
      return std::make_tuple(istart,iend);
    } 
    istart=iend+1;
    continue;

  }
  return std::make_tuple(istart,iend);

}

