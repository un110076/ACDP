/*
Adjoint Code Design Pattern with dco/c++
pattern: SymbolicAdjointNLSEigenSparseLU
author: Uwe Naumann (2018)
*/

#pragma once
#include "../Base/ACDP_Base.hpp"
#include "Eigen/SparseLU"

struct ACDP_SymbolicAdjointNLSEigenSparseLU : ACDP_AdjointBase {

  template <typename T, int N=Eigen::Dynamic>
  using VT=Eigen::Matrix<T,N,1>;

  template <typename T>
  using MT=Eigen::SparseMatrix<T>;

  int ns,np;
  double eps;
  std::unique_ptr<ACDP_PrimalBase<DCO_B>> target;
  std::unique_ptr<ACDP_PrimalBase<DCO_T>> target_t;
  std::unique_ptr<ACDP_PrimalBase<DCO_A>> target_a;
  VT<DCO_B> x,p,r;
  VT<DCO_T> x_t,p_t,r_t;
  VT<DCO_A> x_a,p_a,r_a;
  MT<DCO_B> A;

  ACDP_SymbolicAdjointNLSEigenSparseLU(
    DCO_TT *context_tape, 
    const int& ns, const int& np, const double& eps,
    ACDP_PrimalBase<DCO_B>* T,
    ACDP_PrimalBase<DCO_T>* T_t,
    ACDP_PrimalBase<DCO_A>* T_a
    ) : ACDP_AdjointBase(context_tape), ns(ns), np(np), eps(eps),
        target(T), target_t(T_t), target_a(T_a), 
        x(VT<DCO_B>::Zero(ns)), p(VT<DCO_B>::Zero(np)), r(VT<DCO_B>::Zero(ns)),
        x_t(VT<DCO_T>::Zero(ns)), p_t(VT<DCO_T>::Zero(np)), r_t(VT<DCO_T>::Zero(ns)),
        x_a(VT<DCO_A>::Zero(ns)), p_a(VT<DCO_A>::Zero(np)), r_a(VT<DCO_A>::Zero(ns)),
        A(MT<DCO_B>(ns,ns)) {}

  void link_target() {
    if (!target->linked) {
      for (int i=0;i<ns;i++) target->register_input(x(i));
      for (int i=0;i<np;i++) target->register_input(p(i));
      for (int i=0;i<ns;i++) target->register_output(r(i));
      target->linked=true;
    }
    if (!target_t->linked) {
      for (int i=0;i<ns;i++) target_t->register_input(x_t(i));
      for (int i=0;i<np;i++) target_t->register_input(p_t(i));
      for (int i=0;i<ns;i++) target_t->register_output(r_t(i));
      target_t->linked=true;
    }
    if (!target_a->linked) {
      for (int i=0;i<ns;i++) target_a->register_input(x_a(i));
      for (int i=0;i<np;i++) target_a->register_input(p_a(i));
      for (int i=0;i<ns;i++) target_a->register_output(r_a(i));
      target_a->linked=true;
    }
  }

  void evaluate_augmented_primal() {
    for (int i=0;i<ns+np;i++) target->input_value(i)=input_value(i);
    target->evaluate_primal();
    for (int j=0;j<np;j++) p_t(j)=p(j);
    do {
      std::vector<Eigen::Triplet<DCO_B>> entries;
      for (int i=0;i<ns;i++) {
        for (int j=0;j<ns;j++) x_t(j)=x(j);
        dco::derivative(x_t(i))=1;
        target_t->evaluate_primal();
        for (int j=0;j<ns;j++) 
          if (dco::derivative(r_t(j))!=0) 
            entries.push_back(Eigen::Triplet<DCO_B>(j,i,dco::derivative(r_t(j))));
      } 
      A.setFromTriplets(entries.begin(),entries.end());
      Eigen::SparseLU<MT<DCO_B>> solver;
      solver.analyzePattern(A); 
      solver.factorize(A); 
      r=solver.solve(r);
      x-=r;
      target->evaluate_primal();
    } while (r.norm()>eps);
    for (int i=0;i<ns;i++) output_value(i)=x(i);
  }

  void evaluate_adjoint() {
    for (int i=0;i<ns;i++) { r(i)=-output_adjoint(i); output_adjoint(i)=0; }
    Eigen::SparseLU<MT<DCO_B>> solver;
    solver.analyzePattern(A.transpose());  
    solver.factorize(A.transpose()); 
    r=solver.solve(r);
    for (int i=0;i<ns;i++) x_a(i)=x(i);
    for (int i=0;i<np;i++) p_a(i)=p(i);
    DCO_TT* tape = DCO_TT::create();
    for (int i=0;i<np;i++) tape->register_variable(p_a(i));
    target_a->evaluate_primal();
    for (int i=0;i<ns;i++) tape->register_output_variable(r_a(i));
    for (int i=0;i<ns;i++) dco::derivative(r_a(i))=r(i);
    tape->interpret_adjoint();
    for (int i=0;i<np;i++) input_adjoint(i+ns)+=dco::derivative(p_a(i));
    DCO_TT::remove(tape);
  }
};

