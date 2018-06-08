
#include "../../extra/legndr.h"

auto do380( int& jscr, int& i, int& j, const int& nl, int& nll, std::vector<double>& scr, 
  std::vector<double>& x, std::vector<std::vector<double>>& y, const double& em9,
  double& xlast, double& ylast, int& jnz, double& ulast, double& u2last, 
  double& u3last, std::vector<double>& p ){
        std::cout << 380 << std::endl;
        jscr=7+(j-1)*(nl+1);
        scr[jscr-1]=x[i-1];
        if (y[1-1][i-1] >= em9) {
          scr[1+jscr-1]=sigfig(y[1-1][i-1],9,0);
        }
        else { 
          scr[1+jscr-1]=sigfig(y[1-1][i-1],8,0);
        } // endif
        for ( int il = 1; il < nl; ++ il ){
          scr[il+jscr]=sigfig(y[il][i-1],9,0);
          if (scr[il+jscr] > 1) {
            // only warn for big miss, but always fix the overflow
            //  use this same 1+0.0005 value in aceth
            if (scr[il+jscr] > 1+0.0005) {
              std::cout << "call mess('calcem',strng,'')" << std::endl;
            } // endif
            scr[il+jscr-1]=1;
          } // endif
          if (scr[il+jscr] < -1) {
            // only warn for big miss, but always fix the underflow
            if (scr[il+jscr] < -(1+0.0005)) {
              std::cout << "call mess('calcem',strng,'')" << std::endl;
            } // endif
            scr[il+jscr]=-1;
          } // endif
        } // enddo
        xlast=x[i-1];
        ylast=y[1-1][i-1];
        if (ylast != 0) jnz=j;
        ulast=0;
        u2last=0;
        u3last=0;
        nll=3;
        for ( int il = 1; il < nl; ++il ){
          legndr( y[il][i-1], p );
          ulast  = ulast  + p[2-1];
          u2last = u2last + p[3-1];
          u3last = u3last + p[4-1];
        } // enddo
        ulast  = ulast  * y[1-1][i-1]/(nl-1);
        u2last = u2last * y[1-1][i-1]/(nl-1);
        u3last = u3last * y[1-1][i-1]/(nl-1);
        i=i-1;

}



