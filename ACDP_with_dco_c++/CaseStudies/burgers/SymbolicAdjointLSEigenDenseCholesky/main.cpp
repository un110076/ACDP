/*
Adjoint Code Design Pattern with dco/c++
case study: Burgers equation
author: Uwe Naumann (2018)
*/

#include "ACDP_SymbolicAdjointLSEigenDenseCholesky.hpp"

#include "burgers.hpp"

inline void newton(
    DCO_TT* context_tape,
    const int& m, const double& d, const VT<DCO_A>& y_prev, VT<DCO_A>& y
) {
  int n=y.size();
  const double eps=1e-15;
  MT<DCO_A> A(n,n); A.reserve(3*n-2);
  VT<DCO_A> r(VT<DCO_A>::Zero(n));
  f(m,d,y_prev,y,r);
  do {
    dfdy(m,d,y,A);
    ACDP_SymbolicAdjointLSEigenDenseCholesky::DMT<DCO_A> A_dense(A); 
    ACDP_SymbolicAdjointLSEigenDenseCholesky *e= new ACDP_SymbolicAdjointLSEigenDenseCholesky(context_tape,n);
    e->register_A(A_dense);
    e->register_b(r);
    context_tape->register_acdp(e);
    e->evaluate_augmented_primal();
    y-=r;
    f(m,d,y_prev,y,r);
  } while (r.norm()>eps);
}

inline void euler(
    DCO_TT* context_tape,
    const int& m, const double& d, VT<DCO_A>& y
) {
  for (int j=0;j<m;j++) {
    VT<DCO_A> y_prev=y;
    newton(context_tape,m,d,y_prev,y);
  }
}

int main(int, char* v[]){
  int n=stoi(v[1]), nsteps=stoi(v[2]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Zero(n));
  const double pi=3.141592653589793;
  for (int i=1;i<n-1;i++) y_indep[i]=sin((2*pi*i)/n);
  double d=1e-2;
  DCO_TT* context_tape=DCO_TT::create();
  for (int i=0;i<n;i++) context_tape->register_variable(y_indep(i));
  y=y_indep;
  euler(context_tape,nsteps,d,y);  
  context_tape->register_output_variable(y((n-1)/2));
  dco::derivative(y((n-1)/2))=1.;
  context_tape->interpret_adjoint();
  for(int i=0;i<n;i++)
    cout << dco::derivative(y_indep[i]) << endl;
  DCO_TT::remove(context_tape);
  return 0;
}
