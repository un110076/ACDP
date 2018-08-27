/*
Adjoint Code Design Patterns with dco/c++
case study: Burgers equation
author: Uwe Naumann (2018)
*/

#include <iostream>
#include <cfloat>
#include <cmath>
using namespace std;

#include "diffusion.hpp"

int main(int c, char* v[]){
  if (c!=3) throw;
  int n=stoi(v[1]), m=stoi(v[2]);
  VT<double> y(VT<double>::Zero(n));
  const double pi=3.141592653589793;
  for (int i=1;i<n-1;i++) y[i]=sin((2*pi*i)/n);
  double d=1e-2;
  for(int j=0;j<n;j++) {
    VT<double> yph(y), ymh(y);
    double h=(y[j]<10) ? sqrt(DBL_EPSILON) : sqrt(DBL_EPSILON)*abs(y[j]); 
    ymh[j]-=h;
    euler(m,d,ymh);
    yph[j]+=h;
    euler(m,d,yph);
    cout << j << " " << (yph[(n-1)/2]-ymh[(n-1)/2])/(2*h) << endl;
  }
  return 0;
}
