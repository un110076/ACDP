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
#include "dco.hpp"

typedef double DCO_B;
typedef dco::gt1s<DCO_B>::type DCO_T;
typedef dco::ga1sm<DCO_B> DCO_AM;
typedef DCO_AM::type DCO_A;
typedef DCO_AM::tape_t DCO_TT;
typedef DCO_TT::position_t DCO_TPT;

template<typename T>
struct ACDP_PrimalBase : dco::ACDP::PrimalBase<T> {
  bool linked;
  ACDP_PrimalBase() : linked(false) {}
};

struct ACDP_AdjointBase : dco::ACDP::AdjointBase<DCO_TT,DCO_A> {
  bool linked;
  ACDP_AdjointBase(DCO_TT* tape) 
    : dco::ACDP::AdjointBase<DCO_TT,DCO_A>(tape), linked(false) {}
};
