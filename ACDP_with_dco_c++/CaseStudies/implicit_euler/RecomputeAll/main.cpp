/*
Uwe Naumann hereby disclaims all copyright interest in the 
Adjoint Code Design Patterns (ACDP) software.

Uwe Naumann, Aachen, Germany, 11 March 2019
*/

/*
ACDP is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ACDP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "../implicit_euler.hpp"

#include "ACDP_StaticEvolutionRecomputeAll.hpp"

struct StaticEvolutionRecomputeAll_Target : ACDP_EvolutionStep {
  // passive context
  const double dt, eps; const int n;
  // required in evaluate_augmented_primal and evaluate_adjoint
  DCO_A x, x_indep; std::vector<DCO_A> p;
  DCO_TPT tp;

  StaticEvolutionRecomputeAll_Target(
    DCO_TT *context_tape, 
    const double &dt, const double &eps, const int &n
  ) : ACDP_EvolutionStep(context_tape,n), dt(dt), eps(eps), n(n),
      p(std::vector<DCO_A>(n,0)) {}

  void evaluate_primal() {
    DCO_B& x=input_value(0);
    DCO_B x_prev=x;
    DCO_B& p=input_value(step+1);
    newton(x,x_prev,p,step+1,dt,eps);
    output_value(0)=x;
  }

  void evaluate_augmented_primal() {
    x_indep=input(0);
    DCO_TT *context_tape=dco::tape(x_indep); 
    tp=context_tape->get_position();
    context_tape->register_variable(x_indep);
    p[step]=input_value(step+1);
    context_tape->register_variable(p[step]);
    x=x_indep;
    DCO_A x_prev=x;
    newton(x,x_prev,p[step],step+1,dt,eps);
    output_value(0)=dco::value(x);
  }

  void evaluate_adjoint() {
    DCO_TT *context_tape=dco::tape(x); 
    dco::derivative(x)+=output_adjoint(0); output_adjoint(0)=0;
    context_tape->interpret_adjoint_to(tp);
    input_adjoint(0)+=dco::derivative(x_indep);
    dco::derivative(x_indep)=0;
    input_adjoint(step+1)+=dco::derivative(p[step]);
    dco::derivative(p[step])=0;
    context_tape->reset_to(tp);
  }

};

void implicit_euler(DCO_A &x, std::vector<DCO_A> &p, const double &T, const double &eps) {
  int n=p.size(); double dt=T/n;
  DCO_TT* context_tape=dco::tape(x);
  ACDP_StaticEvolutionRecomputeAll *e=new ACDP_StaticEvolutionRecomputeAll(
    context_tape, 
    new StaticEvolutionRecomputeAll_Target(context_tape,dt,eps,n),n+1,1,n
  );
  e->register_input(x);
  for (auto& v : p) e->register_input(v);
  e->register_output(x); 
  context_tape->register_acdp(e);
  e->link();
  e->evaluate_augmented_primal();
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
    std::cout << "dx/dp[" << i << "]="  
              << dco::derivative(p[i]) << std::endl;
  DCO_TT::remove(context_tape);
  return 0;
}
