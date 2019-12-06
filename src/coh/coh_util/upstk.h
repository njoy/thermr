//#include <Eigen/Dense>
#include <iostream>
#include <range/v3/all.hpp>
template <typename F, typename A, typename R>
auto upstkOLD( const F& e, const A& s, R& stk, int nl, int& i ){
 /*-------------------------------------------------------------------
  * Update the linearization stack with energy e and cross
  * sections s.  Here, i is the current index to the stack in stk,
  * nl is the number of legendre orders in s, and nx is the
  * cycle length in the stack.
  *-------------------------------------------------------------------
  */

  //std::cout <<"----   " <<  stk[0][i-1] << std::endl;
  stk[0][i]=stk[0][i-1];
  stk[0][i-1]=e;
  //std::cout <<"----   " <<  stk[0][0] << std::endl;
  //std::cout <<"----   " <<  stk[1][0] << std::endl;
  //std::cout <<"----   " <<  stk[0][1] << std::endl;
  //std::cout <<"----   " <<  stk[1][1] << std::endl;
  //std::cout <<"----   " <<  stk[0][2] << std::endl;
  //std::cout <<"----   " <<  stk[1][2] << std::endl;


  for ( int j = 0; j < nl; ++ j ){
    stk[1+j][i]=stk[1+j][i-1];
    stk[1+j][i-1]=s[j];
  }
  ++i;
  return;
  std::cout << stk.size()<<s.size() << e << std::endl;

}


template <typename F, typename A, typename R>
auto upstk( const F& e, const A& s, R& stk2, int nl, int& i, int nx ){
 /*-------------------------------------------------------------------
  * Update the linearization stack with energy e and cross
  * sections s.  Here, i is the current index to the stack in stk,
  * nl is the number of legendre orders in s, and nx is the
  * cycle length in the stack.
  *-------------------------------------------------------------------
  */

  int index = 0+i*nx;
  stk2[0+(i  )*nx] = stk2[0+(i-1)*nx];
  stk2[0+(i-1)*nx] = e;

  for ( int j = 0; j < nl; ++ j ){
    stk2[(1+j)+(i)*nx]=stk2[(1+j)+(i-1)*nx];
    stk2[(1+j)+(i-1)*nx]=s[j];
  }
  ++i;
  return;
  //std::cout << stk2.size()<<s.size() << e << nl << std::endl;

}


