#include <iostream>

auto mainLoop(int& nep, int& npage, int&j, int& k, int& ib,
  std::vector<double>& scr, std::vector<double>& yu, double& sum, int& nw, 
  int& ncds, int& nb ){
  int iend, istart = 1;
  while (true) {
         //595 continue
         std::cout << 595 << std::endl;
         iend=nep;
         if ((iend-istart) >= npage/2) iend=istart+npage/2-1;
         j=k-1;
         ib=istart-1;

        //596 continue
         while (true){
         j=j+2;
         ib=ib+1;
         scr[j-1]=yu[1+2*ib-1];
         scr[j+1-1]=yu[2+2*ib-1]*2/sum;
         if (ib < iend) { continue; }// go to 596
         else { break; }
         }
         nw=j+1;
         if (k == 0) {
           std::cout << "go to 597" << std::endl;
           //597 continue
             //call moreio(0,0,nscr,scr,nb,nw)
             if (nb == 0) {
               std::cout << "go to 598" << std::endl;
               ncds=ncds+1+(j*(nep+1)+5)/6;
               break;
             }
             istart=iend+1;
             std::cout << "go to 595" << std::endl;
             continue;
         }
         k=0;
         //call tab1io(0,0,nscr,scr,nb,nw)
         if (nb == 0) {
           std::cout << "go to 598" << std::endl;
           ncds=ncds+1+(j*(nep+1)+5)/6;
           break;
         } 
         istart=iend+1;
         std::cout << "go to 595" << std::endl;
         continue;

        //598 continue
        // ncds=ncds+1+(j*(nep+1)+5)/6;
       }

}

