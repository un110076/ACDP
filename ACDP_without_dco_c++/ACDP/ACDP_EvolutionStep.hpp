/*
Adjoint Code Design Pattern without dco/c++
pattern: EvolutionStep
author: Uwe Naumann (2018)
*/

#pragma once

#include "ACDP_ArgCP.hpp"

struct ACDP_EvolutionStep : ACDP_ArgCP {
  int number_steps, current_step; 
  ACDP_EvolutionStep(const int &ns) : number_steps(ns), current_step(-1) {} 
};
