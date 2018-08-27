/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include "ACDP_EarlyTangentPreaccumulation.hpp"

#include "../diffusion.hpp"

struct EarlyTangentPreaccumulation_Target : ACDP_PrimalBase<DCO_T> {

  int base,n,nsteps,ncs;
  const DCO_B &d;


  EarlyTangentPreaccumulation_Target(
    const int &n, const int &base, const int &nsteps, 
    const int &ncs, const double &d
  ) : base(base), n(n), nsteps(nsteps), ncs(ncs), d(d) {}

  void evaluate_primal() {
    VT<DCO_T> y(n);
    for (int i=0;i<n;i++) y(i)=input_value(i);
    for (int i=base;i<min(base+ncs,nsteps);i++) {
      VT<DCO_T> y_prev=y;
      newton(nsteps,d,y_prev,y);
    }
    for (int i=0;i<n;i++) output_value(i)=y(i);
  }

};

inline void euler(
    const int &nsteps, const int &ncs, const double& d, VT<DCO_A>& y) {
  int n=y.size();
  DCO_TT* context_tape=dco::tape(y[0]);
  for (int base=0;base<nsteps;base+=ncs) {
      ACDP_EarlyTangentPreaccumulation *e= new ACDP_EarlyTangentPreaccumulation(
        context_tape, 
        new EarlyTangentPreaccumulation_Target(n,base,nsteps,ncs,d),
        n,n
      );
      for (int i=0;i<n;i++) {
        e->register_input(y(i));
        e->register_output(y(i));
      }
      context_tape->register_acdp(e);
      e->link_target();
      e->evaluate_augmented_primal();
  }
}

int main(int, char* v[]){
  int n=stoi(v[1]), nsteps=stoi(v[2]), ncs=stoi(v[3]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Zero(n));
  const double pi=3.141592653589793;
  for (int i=1;i<n-1;i++) y_indep[i]=sin((2*pi*i)/n);
  double d=1e-2;
  DCO_TT* context_tape=DCO_TT::create();
  for (int i=0;i<n;i++) context_tape->register_variable(y_indep(i));
  y=y_indep;
  euler(nsteps,ncs,d,y);  
  context_tape->register_output_variable(y((n-1)/2));
  dco::derivative(y((n-1)/2))=1.;
  context_tape->interpret_adjoint();
  for(int i=0;i<n;i++)
    cout << dco::derivative(y_indep[i]) << endl;
  DCO_TT::remove(context_tape);
  return 0;
}
