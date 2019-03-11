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

#include "ACDP_EarlyForwardFiniteDifferences.hpp"

struct EarlyForwardFiniteDifferences_Target : ACDP_PrimalBase<DCO_B> {

  // passive context
  const double dt,eps; const int step;

  EarlyForwardFiniteDifferences_Target(
    const double &dt, const double &eps, const int &i
  ) : dt(dt), eps(eps), step(i+1) {}

  void evaluate_primal() {
    DCO_B x=input_value(0), x_prev=input_value(1), p=input_value(2);
    newton(x,x_prev,p,step,dt,eps); 
    output_value(0)=x;
  }
};

void implicit_euler(DCO_A &x, std::vector<DCO_A> &p, const double &T, const double &eps) {
  int n=p.size(); double dt=T/n;
  for (int i=0;i<n;i++) {
    DCO_A x_prev=x;
    DCO_TT* context_tape=dco::tape(x);
    ACDP_EarlyForwardFiniteDifferences *e= new ACDP_EarlyForwardFiniteDifferences(
      context_tape, 
      new EarlyForwardFiniteDifferences_Target(dt,eps,i),3,1
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
