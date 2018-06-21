
#include <Eigen/Dense>
#include <iostream>

auto upstk( const double& e, const std::vector<double>& s, 
    Eigen::MatrixXd& stk, int nl, int& i ){
 /*-------------------------------------------------------------------
  * Update the linearization stack with energy e and cross
  * sections s.  Here, i is the current index to the stack in stk,
  * nl is the number of legendre orders in s, and nx is the
  * cycle length in the stack.
  *-------------------------------------------------------------------
  */

  stk(0,i)   = stk(0,i-1);
  stk(0,i-1) = e;
  for ( int j = 1; j <= nl; ++ j ){
    stk(j,i)   = stk(j,i-1);
    stk(j,i-1) = s[j-1];
  }

  ++i;
}


