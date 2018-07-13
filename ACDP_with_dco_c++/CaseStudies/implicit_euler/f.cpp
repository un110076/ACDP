/*
Adjoint Code Design Patterns with dco/c++
case study: implicit Euler (primal driver)
author: Uwe Naumann (2018)
*/

#include <iostream>
#include <vector>

#include "implicit_euler.hpp"

int main(int, char*v[]) {
  int n=std::stoi(v[1]);
  double x=1; 
  const double T=1, eps=1e-15;            
  std::vector<double> p(n,1); 
  implicit_euler(x,p,T,eps);
  std::cout << "x=" << x << std::endl;
  return 0;
}
