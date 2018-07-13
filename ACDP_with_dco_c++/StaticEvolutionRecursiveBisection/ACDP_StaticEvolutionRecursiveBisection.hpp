/*
Adjoint Code Design Pattern with dco/c++
pattern: StaticEvolutionRecursiveBisection
author: Uwe Naumann (2018)
*/

#pragma once
#include "../StaticEvolutionRecomputeAll/ACDP_StaticEvolutionRecomputeAll.hpp"

struct ACDP_StaticEvolutionRecursiveBisection : ACDP_ArgCP {
  std::unique_ptr<ACDP_StaticEvolutionRecomputeAll> tgt;
  int n,m,base,ns,nc;
  DCO_TT* context_tape;

  ACDP_StaticEvolutionRecursiveBisection(
    DCO_TT* context_tape, ACDP_EvolutionStep* tgt, int n, int m, int ns, int nc
  ) : ACDP_ArgCP(context_tape), tgt(new ACDP_StaticEvolutionRecomputeAll(context_tape,tgt,n,m,ns)), n(n), m(m), base(0), ns(ns), nc(nc), context_tape(context_tape) {}

  int split(int ns, int) { return ns/2; }

  void link() {
    if (!tgt->linked) {
      for (int i=0;i<n;i++)
        tgt->register_input(input(i));
      for (int i=0;i<m;i++)
        tgt->register_output(output(i));
      tgt->linked=true;
    }
  }

  void evaluate_augmented_primal() { 
    tgt->link();
    push_args();
    tgt->base=base; tgt->ns=ns; 
    tgt->evaluate_primal(); 
  }

  void evaluate_adjoint() { 
    top_args(); pop_args();
    if (ns>2&&nc-1) {
      int base_all=base;
      int ns_all=ns;
      int split_pos=ns=split(ns_all,nc); 
      evaluate_augmented_primal();
      nc--; ns=ns_all-ns;
      push_args(); 
      base+=split_pos;
      evaluate_adjoint();
      nc++; ns=split_pos;
      base=base_all;
      evaluate_adjoint();
      ns=ns_all;
    } else {
      tgt->base=base; 
      tgt->ns=ns; 
      tgt->tgt->push_args(); 
      tgt->evaluate_adjoint();
    }
  }
};
