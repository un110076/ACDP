/*
Adjoint Code Design Patterns with dco/c++
case study: implicit Euler 
author: Uwe Naumann (2018)
*/

#include "../implicit_euler.hpp"

#include "ACDP_EarlyTangentPreaccumulation.hpp"

struct EarlyTangentPreaccumulation_Target : ACDP_PrimalBase<DCO_T> {

  // passive context
  const double dt,eps; const int step;

  EarlyTangentPreaccumulation_Target(
    const double &dt, const double &eps, const int &i
  ) : dt(dt), eps(eps), step(i+1) {}

  void evaluate_primal() {
    DCO_T x=input_value(0), x_prev=input_value(1), p=input_value(2);
    newton(x,x_prev,p,step,dt,eps); 
    output_value(0)=x;
  }
};

void implicit_euler(DCO_A &x, std::vector<DCO_A> &p, const double &T, const double &eps) {
  int n=p.size(); double dt=T/n;
  for (int i=0;i<n;i++) {
    DCO_A x_prev=x;
    DCO_TT* context_tape=dco::tape(x);
    ACDP_EarlyTangentPreaccumulation *e= new ACDP_EarlyTangentPreaccumulation(
      context_tape, 
      new EarlyTangentPreaccumulation_Target(dt,eps,i),3,1
    );
    e->register_input(x); e->register_input(x_prev); e->register_input(p[i]);
    e->register_output(x); 
    context_tape->register_acdp(e);
    e->link_target(); 
    e->evaluate_augmented_primal();
  }
}

int main(int, char* v[]) {
  int n=std::stoi(v[1]);
  const double x0=1, T=1, eps=1e-15;            
  std::vector<DCO_A> p(n,1); 
  DCO_A x=x0;  
  DCO_TT* context_tape=DCO_TT::create();
  context_tape->register_variable(x); 
  context_tape->register_variable(p);
  DCO_TPT tp=context_tape->get_position();
  implicit_euler(x,p,T,eps);
  dco::derivative(x)=1;
  context_tape->interpret_adjoint_and_reset_to(tp);
  for (int i=0;i<n;i++)
    std::cout << "dx/dp[" << i << "]="  << dco::derivative(p[i]) << std::endl;
  DCO_TT::remove(context_tape);
  return 0;
}
