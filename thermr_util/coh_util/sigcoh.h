#include <iostream>
#include <vector>

auto sigcoh( double e, double enext, std::vector<double> s, int nl, int lat, 
  double temp, double emax, int natom ){
 /*-------------------------------------------------------------------
  * Compute the first nl Legendre components of the coherent scatter-
  * ing at energy e from lattice type lat.  Here enext is the next
  * Bragg edge.  Initialize if e=0.  A list of reciprocal lattice
  * shells and weights is precomputed and stored for use at all e.
  * Long, closely-spaced shells are grouped together to speed up the
  * calculation.
  *       lat=1  graphite
  *       lat=2  be
  *       lat=3  beo
  *       lat=10 read from endf6
  * nl returns no. of Bragg edges on initialization call.
  *-------------------------------------------------------------------
  */
  int nord, nw, k, i1m, i1, l1, i2m, i2, l2, i3m, i3, l3, l, i, imax, jmin, j,
    il, lmax, last;
  double amne, econ, tsqx, a, c, amsc, scoh, wal2, wint, x, w1, w2, w3, tsq, 
    tau, w, f, st, sf, blast, re, t2, ulim, phi, elim, u, twopis, c1, c2, 
    recon, scon;
  std::vector<double> p(6);
  int nd = 10;
  std::vector<double> dwf1 { 2.1997, 2.7448, 3.2912, 3.8510, 4.4210, 4.9969, 
    6.1624, 7.3387, 9.6287, 11.992 },
    dwf2 { 3.16663, 3.88842, 4.62944, 5.40517, 6.19880, 7.0042, 8.63665, 10.2865, 0, 0 },
    dwf3 {  2.153, 2.6374, 3.1348, 3.6513, 4.1798, 4.7164, 5.8052, 6.9068, 0, 0 },
    tmp { 29, 40, 50, 600, 700, 800, 1000, 1200, 1600, 2000 };

  
   double gr1 = 2.4573e-8, gr2 = 6.700e-8, gr3 = 12.011e0, gr4 = 5.50e0, 
     be1 = 2.2856e-8, be2 = 3.5832e-8, be3 = 9.01e0, be4 = 7.53e0, 
     beo1 = 2.695e-8, beo2 = 4.39e-8, beo3 = 12.5e0, beo4 = 1.e0, 
     sqrt3 = 1.732050808e0, cw = 0.658173e-15, half = 0.5e0, eps = 0.05e0, 
     zero = 0, hbar = 1.05457266e-27, amu = 1.6605402e-24, amassn = 1.008664904,
     ev = 1.60217733e-12;

   // save k,recon,scon


   // initialize.
   // if (e > 0) go to 210
   twopis = (4*M_PI*M_PI);
   amne = amassn * amu;
   econ = ev * 8 * ( amne / hbar ) / hbar;
   recon = 1 / econ;
   tsqx = econ / 20;
   nord = 2;
   // if (lat == 10) go to 200
   if (lat == 1){
      // graphite constants
      a = gr1;
      c = gr2;
      amsc = gr3;
      scoh = gr4 / natom;
      // wal2 = terp(tmp,dwf1,nd,temp,nord);
   }
   else if (lat == 2) {
      // beryllium constants
      a = be1;
      c = be2;
      amsc = be3;
      scoh = be4/natom;
      // wal2 = terp(tmp,dwf2,nd,temp,nord);
   }
   else if (lat == 3){
      // beryllium oxide constants
      a = beo1;
      c = beo2;
      amsc = beo3;
      scoh = beo4/natom;
      // wal2 = terp(tmp,dwf3,nd,temp,nord);
   } 
   else {
      std::cout << "OH NO! Error over here. Illegal lat value" << std::endl;
      // call error('sigcoh','illegal lat.',' ')
   }
   c1 = 4 / ( 3 * a * a );
   c2 = 1 / ( c * c );
   scon = scoh * ( 16 * M_PI*M_PI )/( 2 * a * a * c * sqrt3 * econ );
   wint = cw * amsc * wal2;
   t2 = hbar / ( 2 *amu * amsc );
   ulim = econ * emax;
   nw = 10000;
   std::vector<double> wrk(nw);


   // compute and sort lattice factors.
   phi = ulim / twopis;
   i1m = a * std::pow(phi,0.5);
   i1m = i1m + 1;
   k = 0;
   /*
   do 185 i1=1,i1m
   l1=i1-1
   i2m=int(half*(l1+sqrt(3*(a*a*phi-l1*l1))))
   i2m=i2m+1
   do 180 i2=i1,i2m
   l2=i2-1
   x=phi-c1*(l1*l1+l2*l2-l1*l2)
   i3m=0
   if (x.gt.zero) i3m=int(c*sqrt(x))
   i3m=i3m+1
   do 175 i3=1,i3m
   l3=i3-1
   w1=2
   if (l1.eq.l2) w1=1
   w2=2
   if (l1.eq.0.or.l2.eq.0) w2=1
   if (l1.eq.0.and.l2.eq.0) w2=half
   w3=2
   if (l3.eq.0) w3=1

   tsq=tausq(l1,l2,l3,c1,c2,twopis)
   if (tsq.le.zero.or.tsq.gt.ulim) go to 160
   tau=sqrt(tsq)
   w=exp(-tsq*t2*wint)*w1*w2*w3/tau
   f=w*form(lat,l1,l2,l3)
   if (k.gt.0.and.tsq.gt.tsqx) go to 150
   k=k+1
   if ((2*k).gt.nw) call error('sigcoh','storage exceeded.',' ')
   wrk(2*k-1)=tsq
   wrk(2*k)=f
   go to 160
  150 continue
   do 155 i=1,k
   if (tsq.lt.wrk(2*i-1).or.tsq.ge.(1+eps)*wrk(2*i-1)) go to 155
   wrk(2*i)=wrk(2*i)+f
   go to 160
  155 continue
   k=k+1
   if ((2*k).gt.nw) call error('sigcoh','storage exceeded.',' ')
   wrk(2*k-1)=tsq
   wrk(2*k)=f
  160 continue
   tsq=tausq(l1,-l2,l3,c1,c2,twopis)
   if (tsq.le.zero.or.tsq.gt.ulim) go to 175
   tau=sqrt(tsq)
   w=exp(-tsq*t2*wint)*w1*w2*w3/tau
   f=w*form(lat,l1,-l2,l3)
   if (k.gt.0.and.tsq.gt.tsqx) go to 165
   k=k+1
   if ((2*k).gt.nw) call error('sigcoh','storage exceeded.',' ')
   wrk(2*k-1)=tsq
   wrk(2*k)=f
   go to 175
  165 continue
   do 170 i=1,k
   if (tsq.lt.wrk(2*i-1).or.tsq.ge.(1+eps)*wrk(2*i-1)) go to 170
   wrk(2*i)=wrk(2*i)+f
   go to 175
  170 continue
   k=k+1
   if ((2*k).gt.nw) call error('sigcoh','storage exceeded.',' ')
   wrk(2*k-1)=tsq
   wrk(2*k)=f
  175 continue
  180 continue
  185 continue
   imax=k-1
   do i=1,imax
      jmin=i+1
      do j=jmin,k
         if (wrk(2*j-1).lt.wrk(2*i-1)) then
            st=wrk(2*i-1)
            sf=wrk(2*i)
            wrk(2*i-1)=wrk(2*j-1)
            wrk(2*i)=wrk(2*j)
            wrk(2*j-1)=st
            wrk(2*j)=sf
         endif
      enddo
   enddo
   k=k+1
   wrk(2*k-1)=ulim
   wrk(2*k)=wrk(2*k-2)
   nw=2*k
   enext=recon*wrk(1)
   enext=sigfig(enext,7,-1)
   ! copy data to global fl array
   allocate(fl(nw))
   do i=1,nw
      fl(i)=wrk(i)
   enddo
   deallocate(wrk)
   nl=k
   return

   !--bragg parameters already read from endf6
  200 continue
   k=int(fl(6))
   nl=k
   blast=0
   scon=1
   enext=fl(9)
   enext=sigfig(enext,7,-1)
   do i=1,nl
      l=1+2*(i-1)
      fl(l)=fl(l+8)*econ
      fl(l+1)=fl(l+9)-blast
      blast=fl(l+9)
   enddo
   return

   !--compute cross sections at this energy
  210 continue
   re=1/e
   do il=1,nl
      s(il)=0
   enddo
   last=0
   do i=1,k
      tsq=fl(2*i-1)
      elim=tsq*recon
      if (elim.ge.e) exit
      f=fl(2*i)
      if (e.gt.emax) f=0
      u=1-2*elim*re
      lmax=nl-1
      call legndr(u,p,lmax)
      do il=1,nl
         s(il)=s(il)+f*p(il)
      enddo
      if (i.eq.k) last=1
   enddo
   do il=1,nl
      s(il)=s(il)*scon*re
   enddo
   if (last.eq.1) elim=emax
   if (elim.gt.emax) elim=emax
   enext=sigfig(elim,7,-1)
   if (e.gt.sigfig(enext,7,-1)) enext=sigfig(elim,7,+1)
   return

   contains

      real(kr) function tausq(m1,m2,m3,c1,c2,twopis)
      integer::m1,m2,m3
      real(kr)::c1,c2,twopis
      tausq=(c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*twopis
      return
      end function tausq

   end subroutine sigcoh

   */
}

