/*
Adjoint Code Design Pattern with dco/c++
case study: diffusion
author: Uwe Naumann
*/

#include "ACDP_SymbolicAdjointNLSEigenSparseLU.hpp"

#include "../diffusion.hpp"

template <typename T, int N=Eigen::Dynamic>
using VT=ACDP_SymbolicAdjointNLSEigenSparseLU::VT<T,N>;

template<typename T>
struct SymbolicAdjointNLSEigenSparseLU_Target : ACDP_PrimalBase<T> {
  int nsteps; double yl, yr;

  SymbolicAdjointNLSEigenSparseLU_Target(const int nsteps, const double &yl, const double &yr) : nsteps(nsteps), yl(yl), yr(yr) {}

  void evaluate_primal() {
    int ns=this->input_count()/2;
    int np=this->input_count()/2;
    VT<T> r(ns), y(ns), y_prev(np);
    for (int i=0;i<ns;i++) y[i]=this->input_value(i);
    for (int i=0;i<np;i++) y_prev[i]=this->input_value(ns+i);
    f(nsteps,y,yl,yr,y_prev,r);
    for (int i=0;i<ns;i++) this->output_value(i)=r[i];
  }

};

// implicit Euler integration with Newton extracted
inline void euler(DCO_TT* tape, const int nsteps, const double &eps, const double &yl, const double &yr, VT<DCO_A>& y) {
    int n=y.size();
    for (int i=0;i<nsteps;i++) {
      VT<DCO_A> y_prev=y;
      SymbolicAdjointNLSEigenSparseLU_Target<DCO_B> *f=new SymbolicAdjointNLSEigenSparseLU_Target<DCO_B>(nsteps,yl,yr); 
      SymbolicAdjointNLSEigenSparseLU_Target<DCO_T> *f_t= new SymbolicAdjointNLSEigenSparseLU_Target<DCO_T>(nsteps,yl,yr);
      SymbolicAdjointNLSEigenSparseLU_Target<DCO_A> *f_a= new SymbolicAdjointNLSEigenSparseLU_Target<DCO_A>(nsteps,yl,yr);
      ACDP_SymbolicAdjointNLSEigenSparseLU *e= 
        new ACDP_SymbolicAdjointNLSEigenSparseLU(tape,n,n,eps,f,f_t,f_a);
      for (int k=0;k<n;k++) e->register_input(y(k));
      for (int k=0;k<n;k++) e->register_input(y_prev(k));
      for (int k=0;k<n;k++) e->register_output(y(k));
      tape->register_acdp(e);
      e->link_target();
      e->evaluate_augmented_primal();
    }
}

int main(int, char* v[]){
  int n=stoi(v[1]), nsteps=stoi(v[2]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Ones(n));
  double yl=0, yr=0;
  const double &eps=1e-15;
  DCO_TT* tape=DCO_TT::create();
  for (int k=0;k<n;k++) tape->register_variable(y_indep(k));
  DCO_TPT tape_begin=tape->get_position();
  y=y_indep;
  euler(tape,nsteps,eps,yl,yr,y);
  dco::derivative(y((n-1)/2))=1.;
  tape->interpret_adjoint_and_reset_to(tape_begin);
  for(int i=0;i<n;i++) cout << dco::derivative(y_indep(i)) << endl;
  DCO_TT::remove(tape);
  return 0;
}
