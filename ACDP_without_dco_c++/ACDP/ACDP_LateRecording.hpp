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

#include<vector>
#include<stack>
#include<cmath>
#include<iostream>
#include<memory>
using namespace std;

#include "ACDP_ArgCP.hpp"

struct ACDP_LateRecording : ACDP_AdjointBase {
  unique_ptr<ACDP_ArgCP> tgt;
  ACDP_LateRecording(ACDP_ArgCP* tgt) : tgt(tgt) {}

  void link_target() {
    if (!tgt->linked) {
      for (int i=0;i<input_count();i++)
        tgt->register_input(input(i));
      for (int i=0;i<output_count();i++)
        tgt->register_output(output(i));
      tgt->linked=true;
    }
  }

  void evaluate_primal() {
    tgt->evaluate_primal();
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

