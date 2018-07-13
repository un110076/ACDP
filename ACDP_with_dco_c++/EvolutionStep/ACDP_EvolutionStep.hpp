/*
Adjoint Code Design Pattern with dco/c++
pattern: EvolutionStep
author: Uwe Naumann (2018)
*/

#pragma once

#include "../ArgCP/ACDP_ArgCP.hpp"

struct ACDP_EvolutionStep : ACDP_ArgCP {
  int nsteps,step; 
  ACDP_EvolutionStep(
    DCO_TT* context_tape, 
    int nsteps
  ) : ACDP_ArgCP(context_tape), nsteps(nsteps) {}
};
