/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include <iostream>
#include <cfloat>
#include <cmath>
using namespace std;

#include "diffusion.hpp"

int main(int c, char* v[]){
  if (c!=4) throw;
  int n=stoi(v[1]), m=stoi(v[2]), ncs=stoi(v[3]);
  VT<double> y(VT<double>::Ones(n));
  double yl=0., yr=0.;
  for(int j=0;j<n;j++) {
    VT<double> yph(y), ymh(y);
    double h=(y[j]<10) ? sqrt(DBL_EPSILON) : sqrt(DBL_EPSILON)*abs(y[j]); 
    ymh[j]-=h;
    euler(m,ncs,yl,yr,ymh);
    yph[j]+=h;
    euler(m,ncs,yl,yr,yph);
    cout << j << " " << (yph[(n-1)/2]-ymh[(n-1)/2])/(2*h) << endl;
  }
  return 0;
}
