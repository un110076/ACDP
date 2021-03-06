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

#include "ACDP_EvolutionStep.hpp"

struct ACDP_StaticEvolutionRecomputeAll : ACDP_AdjointBase {
  std::unique_ptr<ACDP_EvolutionStep> tgt;
  int n,m,base,ns;

  ACDP_StaticEvolutionRecomputeAll(
    ACDP_EvolutionStep* tgt, int n, int m, int ns
  ) : tgt(tgt), n(n), m(m), base(0), ns(ns) {}

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
    for (int i=0;i<ns;i++) {
      tgt->current_step=base+i; 
      tgt->evaluate_primal(); 
    }
  }

  void evaluate_augmented_primal() { 
    tgt->push_args();
    for (int i=0;i<ns;i++) {
      tgt->current_step=base+i; 
      tgt->evaluate_primal(); 
    }
  }
  void evaluate_adjoint() { 
    for (int t=ns;t>0;t--) {
      tgt->top_args();
      for (int i=0;i<t-1;i++) {
        tgt->current_step=base+i; 
        tgt->evaluate_primal();
      }
      tgt->current_step=base+t-1; 
      tgt->evaluate_augmented_primal();
      tgt->evaluate_adjoint();
    }
    tgt->pop_args();
  }
};
