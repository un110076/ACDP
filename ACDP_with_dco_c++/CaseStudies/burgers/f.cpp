/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include <iostream>
#include <vector>

#include "diffusion.hpp"

int main(int c, char* v[]){
  if (c!=3) throw;
  int n=std::stoi(v[1]), m=std::stoi(v[2]);
  VT<double> y(VT<double>::Zero(n)); 
  const double pi=3.141592653589793;
  for (int i=1;i<n-1;i++) y[i]=sin((2*pi*i)/n);
  double d=1e-3;
  euler(m,d,y);
  std::cout << y << std::endl;
  return 0;
}

