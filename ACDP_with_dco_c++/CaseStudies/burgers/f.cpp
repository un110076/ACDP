/*
Adjoint Code Design Patterns with dco/c++
case study: Burgers equation
author: Uwe Naumann (2018)
*/

#include <iostream>
#include <vector>

#include "burgers.hpp"

int main(int c, char* v[]){
  if (c!=3) throw;
  int n=std::stoi(v[1]), m=std::stoi(v[2]);
  VT<double> y(VT<double>::Zero(n)); 
  for (int i=1;i<n-1;i++) y[i]=sin((2*pi*i)/n);
  euler(m,d,y);
  std::cout << y << std::endl;
  return 0;
}

