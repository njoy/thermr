#include <iostream>
#include <range/v3/all.hpp>
template <typename F, typename A, typename R>
void upstkOLD( const F& e, const A& s, R& stk, int nl, int& i ){
 /*-------------------------------------------------------------------
  * Update the linearization stack with energy e and cross
  * sections s.  Here, i is the current index to the stack in stk,
  * nl is the number of legendre orders in s, and nx is the
  * cycle length in the stack.
  *-------------------------------------------------------------------
  */

  stk[0][i]=stk[0][i-1];
  stk[0][i-1]=e;

  for ( int j = 0; j < nl; ++ j ){
    stk[1+j][i]=stk[1+j][i-1];
    stk[1+j][i-1]=s[j];
  }
  ++i;
}


template <typename F, typename A, typename R>
void upstk( const F& e, const A& s, R& stk2, int nl, int& i, int nx ){
 /*-------------------------------------------------------------------
  * Update the linearization stack with energy e and cross
  * sections s.  Here, i is the current index to the stack in stk,
  * nl is the number of legendre orders in s, and nx is the
  * cycle length in the stack.
  *-------------------------------------------------------------------
  */

  stk2[0+(i  )*nx] = stk2[0+(i-1)*nx];
  stk2[0+(i-1)*nx] = e;

  for ( int j = 0; j < nl; ++ j ){
    stk2[(1+j)+(i)*nx]=stk2[(1+j)+(i-1)*nx];
    stk2[(1+j)+(i-1)*nx]=s[j];
  }
  ++i;

}


