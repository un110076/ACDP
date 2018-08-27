/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include <iostream>
#include <vector>

#include "diffusion.hpp"

int main(int c, char* v[]){
  if (c!=4) throw;
  int n=std::stoi(v[1]), m=std::stoi(v[2]), ncs=std::stoi(v[3]);
  VT<double> y(VT<double>::Ones(n)); 
  double yl=0, yr=0;
  euler(m,ncs,yl,yr,y);
  std::cout << y << std::endl;
  return 0;
}

