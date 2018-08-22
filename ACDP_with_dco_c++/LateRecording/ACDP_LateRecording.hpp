/*
Adjoint Code Design Pattern with dco/c++
pattern: LateRecording
author: Uwe Naumann (2018)
*/

#pragma once
#include "../ArgCP/ACDP_ArgCP.hpp"

struct ACDP_LateRecording : ACDP_AdjointBase {
  // passive context
  int n,m;
  // target code
  std::unique_ptr<ACDP_ArgCP> tgt;

  ACDP_LateRecording(
    DCO_TT* context_tape, // context tape
    int n, // number of active inputs
    int m, // number of active outputs
    ACDP_ArgCP* tgt // target code
  ) : ACDP_AdjointBase(context_tape), n(n), m(m), tgt(tgt) {}

  void link_target() {
    if (!tgt->linked) {
      for (int i=0;i<n;i++)
        tgt->register_input(input(i));
      for (int i=0;i<m;i++)
        tgt->register_output(output(i));
      tgt->linked=true;
    }
  }

  void evaluate_augmented_primal() { 
    tgt->push_args();
    tgt->evaluate_primal(); 
  }

  void evaluate_adjoint() { 
    tgt->top_args(); tgt->pop_args();
    tgt->evaluate_augmented_primal(); 
    tgt->evaluate_adjoint();
  }
};

