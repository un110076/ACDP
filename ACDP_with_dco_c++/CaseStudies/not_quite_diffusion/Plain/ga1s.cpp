/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include "dco.hpp"

#include "../diffusion.hpp"

int main(int c, char* v[]){
  typedef dco::ga1sm<double> DCO_AM;
  typedef DCO_AM::type DCO_A;
  typedef DCO_AM::tape_t DCO_TT;
  if (c!=4) throw;
  int n=stoi(v[1]), m=stoi(v[2]), ncs=stoi(v[3]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Ones(n));
  double yl=0, yr=0;
  DCO_TT* tape=DCO_TT::create();
  for (int i=0;i<n;i++) tape->register_variable(y_indep[i]);
  y=y_indep;
  euler(m,ncs,yl,yr,y);  
  tape->register_output_variable(y[(n-1)/2]);
  dco::derivative(y[(n-1)/2])=1.;
  tape->interpret_adjoint();
  for(int i=0;i<n;i++)
    cout << dco::derivative(y_indep[i]) << endl;
  DCO_TT::remove(tape);
  return 0;
}

