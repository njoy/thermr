
#include "terp.h"

auto iel( int mat, int itemp, int iold, int inew, int ne, int nex, 
    std::vector<double>& tempr, std::vector<double>& fl, int za, 
    std::vector<double>& scr, double awr, double emax, int natom, int nbin,
    std::vector<double>& esi 
    ){
  /*-------------------------------------------------------------------
   * Compute the elastic scattering from polyethylene or hydrogen
   * in zirconium hydride using the incoherent approximation.
   * The calcem energy grid is used.  If mat == 11 or 12, built-in
   * parameters from the ENDF/B-III evaluations are used.
   * If mat == 20, material parameters have been read from ENDF6.
   *-------------------------------------------------------------------
   */
   int idis,iex,iet,iu,ix,nj,nr,np,ip,ir,ltt,nb,nw,matdp;
   int n,nup,nup1,isave,nne;
   double dwa,c1,e,c2,u,rc2,temp,tt1,ttn,tnxt,math,mth,mfh,mtref,sb;
   double x1,r1x1,xsec,x2,unow;
   std::vector<double> ex(20), ej(20);
   std::vector<double> xie;
   int nupmax = 10;
   std::vector<double> tmp { 296,400,500,600,700,800,1000,1200 },
     dwh { 8.4795, 9.0854, 9.8196, 10.676, 11.625, 12.643, 14.822, 17.125 },
     dwz { 1.9957, 2.6546, 3.2946, 4.5835, 5.2302, 6.5260, 7.8236 };
   
   double c11a=162.88, c11b=296, c11c=34.957, c11d=350, c11e=40.282, c12a=81.44, c13a=6.3366, up=1.1, dn=.9;

   // initialize
   temp=tempr[itemp];
   nj=nex+1;

   // select material parameters
   if (mat == 11) {
      sb=c11a;
      dwa=c11c+(temp-c11b)*(c11e-c11c)/(c11d-c11b);
   }
   else if (mat == 12) {
      dwa=terp(tmp,dwh,8,temp,3);
      sb=c12a;
   }
   else if (mat == 13) {
      dwa=terp(tmp,dwz,8,temp,3);
      sb=c13a;
   }
   else if (mat == 20) {
      sb=fl[1];
      nr=round(fl[5-1]);
      np=round(fl[6-1]);
      if (np == 1) {
         tt1=fl[7+2*nr-1];
         if (abs(temp-tt1) > temp/10){
           std::cout << "call error('iel',&'bad temperature for debye-waller factor',' ') " << std::endl;
         }
         dwa=fl[8+2*nr];
      }
      else {
         tt1=fl[7+2*nr-1];
         ttn=fl[5+2*nr+2*np-1];
         if (temp < dn*tt1 or temp > up*ttn) { 
           std::cout << "call error('iel',&'bad temperature for debye-waller factor',' ')" << std::endl;
         }
         if (tt1 > temp) fl[7+2*nr-1]=temp;
         if (ttn < temp) fl[5+2*nr+2*np]=temp;
         ip=2;
         ir=1;
         //call terpa(dwa,temp,tnxt,idis,fl,ip,ir)
      }
   } 
   else {
     std::cout << "call error('iel','unknown material identifier.',' ')" << std::endl;
    }
   c1=sb/(2*natom);

   // check on calcem energy grid
   nne=0;
   while (true) {
     if (esi[nne+1-1] == 0) break;
     nne=nne+1;
   }
   //allocate(xie(nne))

   // write head and tab2 records for mf6
   // in lanl format
   ltt=6;
   math=matdp;
   mfh=6;
   mth=mtref+1;
   scr[1]=za;
   scr[2]=awr;
   scr[3]=0;
   scr[4]=ltt; // temporary flag for incoherent inelastic data
   scr[5]=0;
   scr[6]=0;
   nw=6;
   //call contio(0,0,nscr,scr,nb,nw)
   scr[1]=1;
   scr[2]=1;
   scr[3]=-2; // special flag for incoherent inelastic data
   scr[4]=1; // law=1 for incoherent inelastic data
   scr[5]=1;
   scr[6]=2;
   scr[7]=2;
   scr[8]=2;
   scr[9]=1.e-5;
   scr[10]=1;
   scr[11]=emax;
   scr[12]=1;
   nw=12;
   // call tab1io(0,0,nscr,scr,nb,nw)
   scr[1]=temp;
   scr[2]=0;
   scr[3]=3; // LANG=3 is a special flag for equally probable cosines
   scr[5]=1;
   scr[6]=nne;
   scr[7]=nne;
   scr[8]=2;
   nw=8;
   //call tab2io(0,0,nscr,scr,nb,nw)
   n=nbin;
   nup=n;
   if (nup > nupmax) nup=nupmax;
   nup1=nup+1;
   nw=n+8;

   /*
   // compute cross sections on calcem energy grid
   // and equi-probable angles.
   do iex=1,nne
      e=esi(iex)
      c2=2*e*dwa
      u=-1
      rc2=1/c2
      x1=exp(-2*c2)
      r1x1=1/(1-x1)
      xsec=c1*rc2*(1-x1)
      scr[1)=0
      scr[2)=e
      scr[3)=0
      scr[4)=0
      scr[5)=n+2
      scr[6)=n+2
      scr[7)=e
      scr[8)=1
      do iu=1,n
         x2=exp(-c2*(1-u))
         unow=1+rc2*log((1-x1)/n+x2)
         scr[8+iu)=n*rc2*(exp(-c2*(1-unow))*(c2*unow-1)-&
           x2*(c2*u-1))*r1x1
         u=unow
      enddo
      xie(iex)=xsec
      if (iprint == 2) write(nsyso,'(4x,1p,2e12.4,0p,10f9.4)')&
        e,xsec,(scr[8+iu),iu=1,nup)
      if (iprint == 2 and n > nup)&
        write(nsyso,'(28x,10f9.4)') (scr[8+iu),iu=nup1,n)
      call listio(0,0,nscr,scr,nb,nw)
   enddo
   do iex=1,ne
      call finda(iex,ex,nex,iold,bufo,nbuf)
      do ix=1,nex
         ej(ix)=ex(ix)
      enddo
      ej(nj)=terp(esi,xie,nne,ex(1),3)
      if (iex == ne) ej(nj)=0
      iet=iex
      if (iex == ne) iet=-iet
      call loada(iet,ej,nj,inew,bufn,nbuf)
   enddo
   isave=iold
   iold=inew
   inew=isave
   nex=nj
   call asend(0,nscr)
   ncdse=5+ne*((n+11)/6)
   return
   end subroutine iel

   */
}

