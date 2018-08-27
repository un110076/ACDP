/*
Adjoint Code Design Patterns with dco/c++
case study: Burgers equation
author: Uwe Naumann (2018)
*/

#include "dco.hpp"

#include "../diffusion.hpp"

int main(int c, char* v[]){
  typedef dco::ga1sm<double> DCO_AM;
  typedef DCO_AM::type DCO_A;
  typedef DCO_AM::tape_t DCO_TT;
  if (c!=3) throw;
  int n=stoi(v[1]), m=stoi(v[2]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Zero(n));
  const double pi=3.141592653589793;
  for (int i=1;i<n-1;i++) y_indep[i]=sin((2*pi*i)/n);
  double d=1e-2;
  DCO_TT* tape=DCO_TT::create();
  for (int i=0;i<n;i++) tape->register_variable(y_indep[i]);
  y=y_indep;
  euler(m,d,y);  
  tape->register_output_variable(y[(n-1)/2]);
  dco::derivative(y[(n-1)/2])=1.;
  tape->interpret_adjoint();
  for(int i=0;i<n;i++)
    cout << dco::derivative(y_indep[i]) << endl;
  DCO_TT::remove(tape);
  return 0;
}

