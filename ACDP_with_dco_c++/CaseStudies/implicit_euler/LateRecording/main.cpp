/*
Adjoint Code Design Pattern with dco/c++
case study: implicit Euler
author: Uwe Naumann (2018)
*/

#include "../implicit_euler.hpp"

#include "ACDP_LateRecording.hpp"

struct LateRecording_Target : ACDP_ArgCP {
  // passive context
  const double dt,eps; const int n,ncs,base;
  // local active variables
  DCO_A x, x_indep; std::vector<DCO_A> p;
  // local tape 
  DCO_TT *local_tape; 

  LateRecording_Target(
      DCO_TT *context_tape, 
      const double &dt, //  time step
      const double &eps, //  accuracy of Newton
      const int &n, // number of time steps
      const int &ncs, // number of consecutive time steps per chunk
      const int &base // number of first time step in chunk
  ) :
    ACDP_ArgCP(context_tape), dt(dt), eps(eps), n(n), ncs(ncs), 
    base(base), p(std::vector<DCO_A>(n,0)) {}

  void evaluate_primal() {
    for (int j=base;j<std::min(base+ncs,n);j++) {
      DCO_B x_prev=input_value(0);
      newton(input_value(0),x_prev,input_value(j+1),j+1,dt,eps);
    }
    output_value(0)=input_value(0);
  }

  void evaluate_augmented_primal() {
    local_tape=DCO_TT::create();
    x_indep=input_value(0);
    local_tape->register_variable(x_indep);
    for (unsigned int j=1;j<input_count();j++) p[j-1]=input_value(j);
    local_tape->register_variable(p);
    x=x_indep;
    for (int j=base;j<std::min(base+ncs,n);j++) {
      DCO_A x_prev=x;
      newton(x,x_prev,p[j],j+1,dt,eps);
    }
    output_value(0)=dco::value(x);
  }

  void evaluate_adjoint() {
    dco::derivative(x)+=output_adjoint(0); output_adjoint(0)=0;
    local_tape->interpret_adjoint();
    input_adjoint(0)+=dco::derivative(x_indep);
    dco::derivative(x_indep)=0;
    for (unsigned int j=1;j<input_count();j++) {
      input_adjoint(j)+=dco::derivative(p[j-1]);
      dco::derivative(p[j-1])=0;
    }
    DCO_TT::remove(local_tape);
  }

};

void implicit_euler(
    DCO_A& x, 
    std::vector<DCO_A>& p, 
    const double& T, 
    const double& eps, 
    const int &ncs
) {
  int n=p.size(); double dt=T/n;
  for (int base=0;base<n;base+=ncs) {
    DCO_TT* context_tape=dco::tape(x);
    ACDP_LateRecording *e= new ACDP_LateRecording(
      context_tape, 
      n+1, // number of active inputs
      1, // number of active outputs
      new LateRecording_Target(context_tape,dt,eps,n,ncs,base) // target
    );
    e->register_input(x);
    for (int k=0;k<n;k++) e->register_input(p[k]);
    DCO_A y; // avoid i/o aliasing
    e->register_output(y);
    e->link();
    context_tape->register_acdp(e);
    e->evaluate_augmented_primal();
    x=y;
  }
}

int main(int, char* v[]) {
  int n=std::stoi(v[1]), ncs=std::stoi(v[2]);
  const double x0=1, T=1, eps=1e-15;            
  std::vector<DCO_A> p(n,1); 
  DCO_A x=x0;  
  DCO_TT* context_tape=DCO_TT::create();
  context_tape->register_variable(x);
  context_tape->register_variable(p);
  DCO_TPT tp=context_tape->get_position();
  implicit_euler(x,p,T,eps,ncs);
  dco::derivative(x)=1;
  context_tape->interpret_adjoint_and_reset_to(tp);
  for (int i=0;i<n;i++)
    std::cout << "dx/dp[" << i << "]="  
              << dco::derivative(p[i]) << std::endl;
  DCO_TT::remove(context_tape);
  return 0;
}
