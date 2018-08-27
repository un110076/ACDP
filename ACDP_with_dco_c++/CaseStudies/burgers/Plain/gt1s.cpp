/*
Adjoint Code Design Patterns with dco/c++
case study: Burgers equation
author: Uwe Naumann (2018)
*/

#include "dco.hpp"

#include "../diffusion.hpp"

int main(int c, char* v[]){
  typedef dco::gt1s<double>::type DCO_T;
  if (c!=3) throw;
  int n=stoi(v[1]), m=stoi(v[2]);
  double d=1e-2;
  for(int i=0;i<n;i++) { 
    VT<DCO_T> y(VT<DCO_T>::Zero(n));
    const double pi=3.141592653589793;
    for (int i=1;i<n-1;i++) y[i]=sin((2*pi*i)/n);
    dco::derivative(y[i])=1;
    euler(m,d,y);
    cout << dco::derivative(y[(n-1)/2]) << endl;
  }
  return 0;
}
