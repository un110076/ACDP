/*
Adjoint Code Design Patterns with dco/c++
case study: implicit Euler (adjoint driver)
author: Uwe Naumann (2018)
*/

#include "dco.hpp"
typedef dco::ga1sm<double> DCO_AM;
typedef DCO_AM::type DCO_A;
typedef DCO_AM::tape_t DCO_TT;

#include "implicit_euler.hpp"

int main(int, char* v[]) {
  int n=std::stoi(v[1]);
  const double T=1, eps=1e-15;            
  std::vector<DCO_A> p(n,1); DCO_A x=1;  
  DCO_TT* tape=DCO_TT::create();
  tape->register_variable(p);
  implicit_euler(x,p,T,eps);
  dco::derivative(x)=1;
  tape->interpret_adjoint();
  for (int i=0;i<n;i++)
    std::cout << "dx/dp[" << i << "]="  
              << dco::derivative(p[i]) << std::endl;
  DCO_TT::remove(tape);
  return 0;
}
