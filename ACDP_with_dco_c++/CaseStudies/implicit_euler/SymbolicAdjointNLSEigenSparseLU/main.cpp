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

#include "ACDP_SymbolicAdjointNLSEigenSparseLU.hpp"

template <typename T, int N=Eigen::Dynamic>
using VT=ACDP_SymbolicAdjointNLSEigenSparseLU::VT<T,N>;

template <typename T>
using MT=ACDP_SymbolicAdjointNLSEigenSparseLU::MT<T>;

template<typename T>
struct SymbolicAdjointNLSEigenSparseLU_Target : ACDP_PrimalBase<T> {

  const double dt;
  const int i;

  SymbolicAdjointNLSEigenSparseLU_Target(const double& dt, const int i) : dt(dt), i(i) {}

  void evaluate_primal() {
    const T* x=&(this->input_value(0)); 
    const T* x_prev=&(this->input_value(1)); 
    const T* p=&(this->input_value(2));
    T r=f(*x,*x_prev,*p,i,dt);
    this->output_value(0)=r;
  }

};

void implicit_euler(DCO_A& x, VT<DCO_A>& p, const double& T, const double &eps) {
  int n=p.size();
  double dt=T/n;
  for (int i=0;i<n;i++) {
    DCO_A x_prev=x;
    DCO_TT *tape=dco::tape(x);
    SymbolicAdjointNLSEigenSparseLU_Target<DCO_B> *f=new SymbolicAdjointNLSEigenSparseLU_Target<DCO_B>(dt,i+1);
    SymbolicAdjointNLSEigenSparseLU_Target<DCO_T> *f_t=new SymbolicAdjointNLSEigenSparseLU_Target<DCO_T>(dt,i+1);
    SymbolicAdjointNLSEigenSparseLU_Target<DCO_A> *f_a=new SymbolicAdjointNLSEigenSparseLU_Target<DCO_A>(dt,i+1);
    ACDP_SymbolicAdjointNLSEigenSparseLU *e=
      new ACDP_SymbolicAdjointNLSEigenSparseLU(tape,1,2,eps,f,f_t,f_a);
    e->register_input(x);
    e->register_input(x_prev);
    e->register_input(p[i]);
    e->register_output(x);
    e->link_target();
    tape->register_acdp(e);
    e->evaluate_augmented_primal();
  }
}

int main(int, char* v[]) {
  int n=std::stoi(v[1]);
  const double x0=1, T=1, eps=1e-15;            
  VT<DCO_A> p(VT<DCO_A>::Ones(n)); 
  DCO_A x=x0;  
  DCO_TT* tape=DCO_TT::create();
  tape->register_variable(x); 
  for (int i=0;i<n;i++) tape->register_variable(p[i]);
  DCO_TPT tpos=tape->get_position();
  implicit_euler(x,p,T,eps);
  tape->register_output_variable(x);
  dco::derivative(x)=1;
  tape->interpret_adjoint_and_reset_to(tpos);
  for (int i=0;i<n;i++)
    std::cout << "dx/dp[" << i << "]="  << dco::derivative(p[i]) << std::endl;
  DCO_TT::remove(tape);
  return 0;
}
