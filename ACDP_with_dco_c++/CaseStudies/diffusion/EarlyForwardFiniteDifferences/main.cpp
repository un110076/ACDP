/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include "ACDP_EarlyForwardFiniteDifferences.hpp"

#include "../diffusion.hpp"

struct EarlyForwardFiniteDifferences_Target : ACDP_PrimalBase<DCO_B> {

  int base,n,nsteps,ncs;
  const DCO_B &yl,&yr;


  EarlyForwardFiniteDifferences_Target(
    const int &n, const int &base, const int &nsteps, 
    const int &ncs, const double &yl, const double &yr
  ) : base(base), n(n), nsteps(nsteps), ncs(ncs), yl(yl), yr(yr) {}

  void evaluate_primal() {
    VT<DCO_B> y(n);
    for (int i=0;i<n;i++) y(i)=input_value(i);
    for (int i=base;i<min(base+ncs,nsteps);i++) {
      VT<DCO_B> y_prev=y;
      newton(nsteps,y_prev,yl,yr,y);
    }
    for (int i=0;i<n;i++) output_value(i)=y(i);
  }

};

template <>
inline void euler(
    const int &nsteps, const int &ncs, const double& yl, 
    const double& yr, VT<DCO_A>& y) {
  int n=y.size();
  DCO_TT* context_tape=dco::tape(y[0]);
  for (int base=0;base<nsteps;base+=ncs) {
      ACDP_EarlyForwardFiniteDifferences *e= new ACDP_EarlyForwardFiniteDifferences(
        context_tape, 
        new EarlyForwardFiniteDifferences_Target(n,base,nsteps,ncs,yl,yr),
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
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Ones(n));
  double yl=0, yr=0;
  DCO_TT* context_tape=DCO_TT::create();
  for (int i=0;i<n;i++) context_tape->register_variable(y_indep(i));
  y=y_indep;
  euler(nsteps,ncs,yl,yr,y);  
  context_tape->register_output_variable(y((n-1)/2));
  dco::derivative(y((n-1)/2))=1.;
  context_tape->interpret_adjoint();
  for(int i=0;i<n;i++)
    cout << dco::derivative(y_indep[i]) << endl;
  DCO_TT::remove(context_tape);
  return 0;
}

