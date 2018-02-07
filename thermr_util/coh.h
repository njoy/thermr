#include <vector>
#include <tuple>

auto coh( int lat, int inew, int ne, int nex, double temp ){
 /*-------------------------------------------------------------------
  * Compute the coherent scattering cross sections for a crystalline
  * material.  The cross section is computed on an energy grid
  * chosen adaptively to represent the actual function within a
  * given tolerance with linear-linear interpolation.
  *-------------------------------------------------------------------
  */
  int nl, imax, nx, nj, i, nlt, nlt1, nb, nw, nbragg, ix, j, iex, il, isave, ltt;
  int nlmax = 6;
  double e, enext, en, xm, ym, test;
  std::vector<double> ej(20), ex(20), x(5), y(5), z(5), b(12);
  double half = 0.5, small = 1.0e-10, tolmin = 1.0e-6, eps = 3.0e-5; 

   // initialize.
   nl = 1;

   // temp = tempr(itemp) --> Just pass in the actual temperature 
   imax = 20;
   if (nl > nlmax) { 
     std::cout << "too many legendre orders" << std::endl; 
     return; 
   }

   nx = nl + 1;
   nj = nex + nl;
   for ( int i = 0; i < nl; ++i ){
      ex[nex+i] = 0;
   }
   nlt = 5;
   if (ne < nlt) nlt = ne;
   nlt1 = nlt - 1;

   //allocate(stk(nx,imax)) --> Use eigen? 

   // determine the energy grid adaptively and
   // store the cross sections in a scratch file.
   nbragg = nl;
   e = 0;
   //sigcoh(e,enext,s,nbragg,lat,temp,emax,natom) --> appears later in thermr
   ix = 1;
   j = 0;
   iex = 0;

  // 100 continue
   iex = iex + 1;
   // finda(iex,ex,nex,iold,bufo,nbuf) --> from util
   
   // if (ex(1).gt.enext*(1+small)) go to 105

   x[1] = ex[1];
   y[1] = ex[2];
   z[1] = ex[nex];
   j = j + 1;
   // call loada(j,ex,nj,inew,bufn,nbuf) --> from util
   
   // go to 100
   
  // 105 continue
   ix = ix + 1;
   x[ix] = ex[1];
   y[ix] = ex[2];
   z[ix] = ex[nex];
   //if (ix.lt.nlt) go to 100
   // prime stack with first bragg edge
   e = enext;
   // sigcoh(e,enext,s,nl,lat,temp,emax,natom) --> elsewhere in thermr
   
   // stk(1,1)=e
   
   for ( int i1 = 0; i1 < nl; ++i1 ){
      // stk(1+il,1)=s(il)
   } 
   i = 1;
   // add next bragg edge to stack
   
  // 120 continue
   e = enext;
   // sigcoh(e,enext,s,nl,lat,temp,emax,natom) --> from thermr
   // upstk(e,s,nl,nx,i,stk) --> from thermr
   // make sure input grid points are included
   
   /*
  125 continue
   
   do 130 ix=1,nlt
   if (x(ix).gt.stk(1,i)*(1+small)) go to 135
  130 continue
   go to 140
  135 continue
   if (x(ix).ge.stk(1,i-1)*(1-small)) go to 140
   e=x(ix)
   call sigcoh(e,en,s,nl,lat,temp,emax,natom)
   call upstk(e,s,nl,nx,i,stk)
   ! compare linear approximation to true function.
  140 continue
   if (i.eq.imax) go to 160
   xm=half*(stk(1,i-1)+stk(1,i))
   xm=sigfig(xm,7,0)
   if (stk(1,i-1)-stk(1,i).lt.eps*xm) go to 160
   call sigcoh(xm,en,s,nl,lat,temp,emax,natom)
   do 150 il=1,nl
   call terp1(stk(1,i),stk(1+il,i),&
     stk(1,i-1),stk(1+il,i-1),xm,ym,2)
   test=tol*abs(s(il))
   if (test.lt.tolmin) test=tolmin
   if (abs(s(il)-ym).gt.test) go to 210
  150 continue
   ! all components pass.  save top point in stack and continue.
  160 continue
   j=j+1
   ej(1)=stk(1,i)
  170 continue
   if (ej(1).le.x(3)*(1+small).or.iex.eq.ne) go to 190
   do ix=1,nlt1
      x(ix)=x(ix+1)
      y(ix)=y(ix+1)
      z(ix)=z(ix+1)
   enddo
   iex=iex+1
   call finda(iex,ex,nex,iold,bufo,nbuf)
   x(nlt)=ex(1)
   y(nlt)=ex(2)
   z(nlt)=ex(nex)
   if (iex.eq.ne) nlt=nlt-1
   if (iex.eq.ne) nlt1=nlt1-1
   go to 170
  190 continue
   ej(2)=terp(x,y,nlt,ej(1),nlt1)
   ej(nex)=terp(x,z,nlt,ej(1),nlt1)
   if (stk(1,i).gt.emax*(1+small)) ej(nex)=0
   do il=1,nl
      ej(nex+il)=stk(1+il,i)
      if (stk(1,i).gt.emax*(1+small)) ej(nex+il)=0
   enddo
   if (stk(1,i).gt.emax*(1+small)) go to 230
   call loada(j,ej,nj,inew,bufn,nbuf)
   i=i-1
   if (i.gt.1) go to 125
   go to 120
   ! test fails.  add point to stack and continue.
  210 continue
   call upstk(xm,s,nl,nx,i,stk)
   go to 125
   ! linearization complete.  save last point.
  230 continue
   ne=j
   j=-j
   call loada(j,ej,nj,inew,bufn,nbuf)
   isave=iold
   iold=inew
   inew=isave
   nex=nj
   ltt=7
   b(1)=za
   b(2)=awr
   b(3)=0
   b(4)=ltt ! temporary flag
   b(5)=0
   b(6)=0
   math=matdp
   mfh=6
   mth=mtref+1
   call contio(0,0,nscr,b,nb,nw)
   b(1)=1
   b(2)=1
   b(3)=-nbragg ! use lip field for nbragg
   b(4)=0 ! law=0 for coherent elastic data
   b(5)=1
   b(6)=2
   b(7)=2
   b(8)=2
   b(9)=1.e-5_kr
   b(10)=1
   b(11)=emax
   b(12)=1
   nw=12
   call tab1io(0,0,nscr,b,nb,nw)
   ncdse=3
   call asend(0,nscr)
   deallocate(stk)
   return
   end subroutine coh
   */
}

