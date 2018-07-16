#include<stack>
#include<cmath>
#include<iostream>
#include<memory>

#include "ACDP_StaticEvolutionRecursiveBisection.hpp"

inline void primal_euler(double &x, const int nsteps, const double& dt) {
  for (int i=0;i<nsteps;i++) x+=dt*sin(x);
}

inline void augmented_primal_euler(double &x, const int nsteps, const double& dt, std::stack<double>& tbr) {
  for (int j=0;j<nsteps;j++) {
    tbr.push(1+dt*cos(x)); // recording partial derivatives
    x+=dt*sin(x);
  }
}

inline void adjoint_euler(double &x_a, const int nsteps, std::stack<double>& tbr) {
  for (int j=nsteps;j>=1;j--) {
    x_a*=tbr.top(); tbr.pop();
  }
}

struct StaticEvolutionRecursiveBisectionTarget : ACDP_EvolutionStep {

  // passive context
  const double& dt;
  const int& nsteps;
  // tbr 
  std::stack<double> tbr;

  StaticEvolutionRecursiveBisectionTarget(const double& dt, const int& nsteps) :
    ACDP_EvolutionStep(nsteps), dt(dt), nsteps(nsteps) {}

  void evaluate_primal() {
    double x=input_value(0);
    primal_euler(x,1,dt);
    output_value(0)=x;
  }

  void evaluate_augmented_primal() {
    double x=input_value(0);
    augmented_primal_euler(x,1,dt,tbr);
    output_value(0)=x;
  }

  void evaluate_adjoint() {
    double x_a=0;
    x_a+=output_adjoint(0); output_adjoint(0)=0;
    adjoint_euler(x_a,1,tbr);
    input_adjoint(0)+=x_a; x_a=0;
  }
};

int main(int c, char* v[]) {
  if (c!=3) throw;
  int nsteps=std::stoi(v[1]);
  int ncp=std::stoi(v[2]);
  AT x(AT(1,1));
  double dt=1.0/(3*nsteps);
  std::stack<double> tbr; 

  // augmented primal
  augmented_primal_euler(x.first,nsteps,dt,tbr);
  std::unique_ptr<ACDP_StaticEvolutionRecursiveBisection> e(
    new ACDP_StaticEvolutionRecursiveBisection(new StaticEvolutionRecursiveBisectionTarget(dt,nsteps),1,1,nsteps,ncp)
  );
  e->register_input(x); 
  e->register_output(x); 
  e->link_target();
  e->evaluate_augmented_primal();
  augmented_primal_euler(x.first,nsteps,dt,tbr);
  std::cout << "x=" <<  x.first << std::endl;

  // adjoint
  adjoint_euler(x.second,nsteps,tbr);
  e->evaluate_adjoint();
  adjoint_euler(x.second,nsteps,tbr);
  std::cout << "x_a=" <<  x.second << std::endl;
  return 0;
}
