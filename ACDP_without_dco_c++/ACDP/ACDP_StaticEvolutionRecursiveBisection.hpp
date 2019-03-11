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

#pragma once
#include "ACDP_StaticEvolutionRecomputeAll.hpp"

struct ACDP_StaticEvolutionRecursiveBisection : ACDP_ArgCP {
  std::unique_ptr<ACDP_StaticEvolutionRecomputeAll> tgt;
  int n,m,base,ns,nc;

  ACDP_StaticEvolutionRecursiveBisection(
    ACDP_EvolutionStep* tgt, int n, int m, int ns, int nc
  ) : tgt(new ACDP_StaticEvolutionRecomputeAll(tgt,n,m,ns)), n(n), m(m), base(0), ns(ns), nc(nc) {}

  int split(int ns, int) { return ns/2; }

  void link_target() {
    if (!tgt->linked) {
      for (int i=0;i<n;i++)
        tgt->register_input(input(i));
      for (int i=0;i<m;i++)
        tgt->register_output(output(i));
      tgt->linked=true;
    }
  }

  void evaluate_primal() {
    tgt->link_target();
    for (int i=0;i<ns;i++) {
      tgt->base=base; tgt->ns=ns; 
      tgt->evaluate_primal();
    }
  }

  void evaluate_augmented_primal() { 
    tgt->link_target();
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
