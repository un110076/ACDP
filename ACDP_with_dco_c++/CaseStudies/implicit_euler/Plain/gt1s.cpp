/*
Adjoint Code Design Patterns with dco/c++
case study: implicit Euler (tangent driver)
author: Uwe Naumann (2018)
*/

#include "dco.hpp"
#include "implicit_euler.hpp"

int main(int, char* v[]) {
  int n=std::stoi(v[1]);
  typedef dco::gt1s<double>::type DCO_T;
  const double x0=1, T=1, eps=1e-15;            
  std::vector<DCO_T> p(n,1); 
  for (int i=0;i<n;i++) {
    DCO_T x=x0;
    dco::derivative(p[i])=1.0;
    implicit_euler(x,p,T,eps);
    std::cout << "dx/dp[" << i << "]="  
              << dco::derivative(x) << std::endl;
    dco::derivative(p[i])=0.0;
  }
  return 0;
}
