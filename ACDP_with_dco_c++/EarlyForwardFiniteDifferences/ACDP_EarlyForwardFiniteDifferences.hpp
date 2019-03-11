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

#include<cfloat>
#include "../Base/ACDP_Base.hpp"

struct ACDP_EarlyForwardFiniteDifferences : public ACDP_AdjointBase {
  std::unique_ptr<ACDP_PrimalBase<DCO_B>> tgt;
  int n,m;
  std::vector<DCO_B> jac,x,y;

  ACDP_EarlyForwardFiniteDifferences(
    DCO_TT *context_tape, ACDP_PrimalBase<DCO_B>* tgt, int n, int m
  ) :
    ACDP_AdjointBase(context_tape), tgt(tgt), n(n), m(m), jac(std::vector<DCO_B>(m*n,0)), x(std::vector<DCO_B>(n,0)), y(std::vector<DCO_B>(m,0)) {}
  
  void link_target() {
    if (!tgt->linked) {
      for (auto& v : x) tgt->register_input(v);
      for (auto& v : y) tgt->register_output(v);
      tgt->linked=true;
    }
  }

  void evaluate_augmented_primal() {
    for (int i=0;i<n;i++) tgt->input_value(i)=input_value(i);
    tgt->evaluate_primal();
    std::vector<DCO_B> y_s(m,0);
    for (int i=0;i<m;i++) y_s[i]=tgt->output_value(i);
    for (int i=0;i<n;i++) {
      for (int j=0;j<n;j++) tgt->input_value(j)=input_value(j);
      DCO_B h=sqrt(DBL_EPSILON); 
      h*=tgt->input_value(i)>1 ? abs(tgt->input_value(i)):1;
      tgt->input_value(i)+=h;
      tgt->evaluate_primal();
      for (int j=0;j<m;j++) 
        jac[j*n+i]=(tgt->output_value(j)-y_s[j])/h;
    }
    for (int j=0;j<m;j++) output_value(j)=y_s[j];
  }

 void evaluate_adjoint() {
    std::vector<DCO_B> tmp(m);
    for (int j=0;j<m;j++) {
      tmp[j]=output_adjoint(j); output_adjoint(j)=0;
    }
    for (int i=0;i<n;i++)
      for (int j=0;j<m;j++) {
        input_adjoint(i)+=jac[j*n+i]*tmp[j];
      }
  }

};
