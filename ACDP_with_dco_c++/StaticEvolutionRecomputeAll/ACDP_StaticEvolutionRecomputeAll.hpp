/*
Adjoint Code Design Pattern with dco/c++
pattern: StaticEvolutionRecomputeAll
author: Uwe Naumann (2018)
*/

#pragma once

#include "../EvolutionStep/ACDP_EvolutionStep.hpp"

struct ACDP_StaticEvolutionRecomputeAll : ACDP_AdjointBase {
  std::unique_ptr<ACDP_EvolutionStep> tgt;
  int n,m,base,ns;

  ACDP_StaticEvolutionRecomputeAll(
    DCO_TT* context_tape, ACDP_EvolutionStep* tgt, int n, int m, int ns
  ) : ACDP_AdjointBase(context_tape), tgt(tgt), n(n), m(m), base(0), ns(ns) {}

  void link() {
    if (!tgt->linked) {
      for (int i=0;i<n;i++)
        tgt->register_input(input(i));
      for (int i=0;i<m;i++)
        tgt->register_output(output(i));
      tgt->linked=true;
    }
  }

  void evaluate_primal() { 
    for (int i=0;i<ns;i++) {
      tgt->step=base+i; 
      tgt->evaluate_primal(); 
    }
  }

  void evaluate_augmented_primal() { 
    tgt->push_args();
    for (int i=0;i<ns;i++) {
      tgt->step=base+i; 
      tgt->evaluate_primal(); 
    }
  }
  void evaluate_adjoint() { 
    for (int t=ns;t>0;t--) {
      tgt->top_args();
      for (int i=0;i<t-1;i++) {
        tgt->step=base+i; 
        tgt->evaluate_primal();
      }
      tgt->step=base+t-1; 
      tgt->evaluate_augmented_primal();
      tgt->evaluate_adjoint();
    }
    tgt->pop_args();
  }
};
