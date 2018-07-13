/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include "ACDP_StaticEvolutionRecomputeAll.hpp"

#include "../diffusion.hpp"

struct StaticEvolutionRecomputeAll_Target : ACDP_EvolutionStep {

  int n;
  DCO_B yl,yr;
  VT<DCO_A> y_indep,y;
  DCO_TT *local_tape;

  StaticEvolutionRecomputeAll_Target(
    DCO_TT *context_tape, const int &nsteps, const int &n, 
    const DCO_B &yl, const DCO_B &yr
  ) : ACDP_EvolutionStep(context_tape, nsteps), n(n), yl(yl), yr(yr),
      y_indep(VT<DCO_A>::Zero(n)), y(VT<DCO_A>::Zero(n)) {}

  void evaluate_primal() {
    VT<DCO_B> y(n);
    for (int i=0;i<n;i++) y(i)=input_value(i);
    VT<DCO_B> y_prev=y;
    newton(nsteps,y_prev,yl,yr,y);
    for (int i=0;i<n;i++) output_value(i)=y(i);
  }

  void evaluate_augmented_primal() {
    for (int i=0;i<n;i++) y_indep(i)=input_value(i);
    local_tape=DCO_TT::create();
    for (int i=0;i<n;i++) local_tape->register_variable(y_indep(i));
    y=y_indep;
    VT<DCO_A> y_prev=y;
    newton(nsteps,y_prev,yl,yr,y);
    for (int i=0;i<n;i++) output_value(i)=dco::value(y(i));
  }

  void evaluate_adjoint() {
    for (int i=0;i<n;i++) {
      dco::derivative(y(i))+=output_adjoint(i);
      output_adjoint(i)=0;
    }
    local_tape->interpret_adjoint();
    for (int i=0;i<n;i++) {
      input_adjoint(i)+=dco::derivative(y_indep[i]);
      dco::derivative(y_indep[i])=0;
    }
    DCO_TT::remove(local_tape); 
  }

};

template <>
inline void euler(
    const int &nsteps, const int&, const double &yl, const double &yr, 
    VT<DCO_A>& y
) {
  int n=y.size();
  DCO_TT* context_tape=dco::tape(y[0]);
  ACDP_StaticEvolutionRecomputeAll *e=new ACDP_StaticEvolutionRecomputeAll(
    context_tape,
    new StaticEvolutionRecomputeAll_Target(context_tape,nsteps,n,yl,yr),
    n,n,nsteps
  );
  for(int i=0;i<n;i++) {
    e->register_input(y(i));
    e->register_output(y(i));
  }
  e->link();
  context_tape->register_acdp(e);
  e->evaluate_augmented_primal();
}

int main(int, char* v[]){
  int n=stoi(v[1]), nsteps=stoi(v[2]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Ones(n));
  DCO_B yl=0,yr=0;
  DCO_TT* context_tape=DCO_TT::create();
  for(int i=0;i<n;i++) 
    context_tape->register_variable(y_indep(i));
  DCO_TPT tpos=context_tape->get_position();
  y=y_indep;
  euler(nsteps,1,yl,yr,y);
  context_tape->register_output_variable(y((n-1)/2));
  dco::derivative(y((n-1)/2))=1.;
  context_tape->interpret_adjoint_and_reset_to(tpos);
  for(int i=0;i<n;i++) 
    std::cout << dco::derivative(y_indep[i]) << std::endl;
  DCO_TT::remove(context_tape);
  return 0;
}
