/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include "dco.hpp"

#include "../diffusion.hpp"

int main(int c, char* v[]){
  typedef dco::gt1s<double>::type DCO_T;
  if (c!=4) throw;
  int n=stoi(v[1]), m=stoi(v[2]), ncs=stoi(v[3]);
  double yl=0, yr=0;
  for(int i=0;i<n;i++) { 
    VT<DCO_T> y(VT<DCO_T>::Ones(n)); 
    dco::derivative(y[i])=1;
    euler(m,ncs,yl,yr,y);
    cout << dco::derivative(y[(n-1)/2]) << endl;
  }
  return 0;
}
