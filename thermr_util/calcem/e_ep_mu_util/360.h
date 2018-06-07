

auto do360(int& j, const int& jmax, std::vector<double>& xsi, 
    std::vector<double>& x, double& xlast, double& ylast, double& ulast,
    double& u2last, double& u3last, const double& tolmin,
    std::vector<std::vector<double>>& y, std::vector<double>& p2,
    std::vector<double>& p3,
    int& nll, const int& nl, 
    //double& uu, double& u2, double& u3, int& nll, const int& nl, 
    std::vector<double>& p, std::vector<double>& ubar, const int ie, 
    const int i ){
  // 360 continue
  double uu, u2, u3;
  std::cout << 360 << std::endl;
  j=j+1;
  if (j >= jmax) std::cout << "call error('calcem','storage exceeded.',' ')" << "       " << j <<"     " << jmax  << std::endl;
  if (j > 1) xsi[ie-1]=xsi[ie-1]+(x[i-1]-xlast)*(y[0][i-1]+ylast)*0.5;
  if (j > 1) {
    uu=0;
    u2=0;
    u3=0;
    nll=3;
    for ( int il = 1; il < nl; ++il ){
      legndr( y[il][i-1], p );
      uu=uu+p[2-1];
      u2=u2+p[3-1];
      u3=u3+p[4-1];
    } // enddo
    uu=uu/(nl-1);
    uu=uu*y[0][i-1];
    u2=u2/(nl-1);
    u3=u3/(nl-1);
    u2=u2*y[0][i-1];
    u3=u3*y[0][i-1];
    ubar[ie-1]=ubar[ie-1]+0.5*(x[i-1]-xlast)*(uu+ulast);
    p2[ie-1]=p2[ie-1]+0.5*(x[i-1]-xlast)*(u2+u2last);
    p3[ie-1]=p3[ie-1]+0.5*(x[i-1]-xlast)*(u3+u3last);
  } // endif
  // if (j != 3 or xsi[ie-1] >= tolmin) go to 380
  if (j == 3 and xsi[ie-1] < tolmin){ j = 2; }
  // j=2

}



