/*
Adjoint Code Design Patterns with dco/c++
case study: implicit Euler (forward finite differences driver)
author: Uwe Naumann (2018)
*/

#include <iostream>
#include <vector>
#include<cfloat>

#include "implicit_euler.hpp"

int main(int, char* v[]) {
  int n=std::stoi(v[1]);
  const double x0=1, T=1, eps=1e-15;            
  std::vector<double> p(n,1);
  double x=x0; 
  implicit_euler(x,p,T,eps);
  for (int i=0;i<n;i++) {
    double h=sqrt(DBL_EPSILON);
    double xp=x0;
    p[i]+=h;
    implicit_euler(xp,p,T,eps);
    std::cout << "dx/dp[" << i << "]=" << (xp-x)/h << std::endl;
    p[i]-=h;
  }
  return 0;
}

