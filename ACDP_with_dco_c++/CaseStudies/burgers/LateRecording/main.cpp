/*
Adjoint Code Design Patterns with dco/c++
case study: Burgers equation
author: Uwe Naumann (2018)
*/

#include "ACDP_LateRecording.hpp"

#include "burgers.hpp"

struct LateRecording_Target : ACDP_ArgCP {

  int k,m,ncs;
  const DCO_B &d;
  VT<DCO_A> y_indep,y;
  DCO_TT *local_tape;

  LateRecording_Target(DCO_TT *context_tape, const int k, const int m, const int ncs, const DCO_B& d) :
    ACDP_ArgCP(context_tape), k(k), m(m), ncs(ncs), d(d) {}

  void evaluate_primal() {
    int n=input_count();
    VT<DCO_B> y(n);
    for (int i=0;i<n;i++) y(i)=input_value(i);
    for (int i=k;i<min(k+ncs,m);i++) {
      VT<DCO_B> y_prev=y;
      newton(m,d,y_prev,y);
    }
    for (int i=0;i<n;i++)
      output_value(i)=y(i);
  }

  void evaluate_augmented_primal() {
    int n=input_count();
    y_indep=VT<DCO_A>::Zero(n);
    local_tape=DCO_TT::create();
    for (int i=0;i<n;i++) {
      y_indep(i)=input_value(i);
      local_tape->register_variable(y_indep(i));
    }
    y=y_indep;
    for (int i=k;i<min(k+ncs,m);i++) {
      VT<DCO_A> y_prev=y;
      newton(m,d,y_prev,y);
    }
    for (int i=0;i<n;i++)
      output_value(i)=dco::value(y(i));
  }

  void evaluate_adjoint() {
    int n=input_count();
    for (int i=0;i<n;i++) {
      dco::derivative(y(i))+=output_adjoint(i);
      output_adjoint(i)=0;
    }
    local_tape->interpret_adjoint();
    for (int i=0;i<n;i++) {
      input_adjoint(i)+=dco::derivative(y_indep(i));
      dco::derivative(y_indep(i))=0;
    }
    DCO_TT::remove(local_tape);
  }

};

// implicit Euler integration
inline void euler(const int& m, const int& ncs, const DCO_B& d, VT<DCO_A>& y) {
  DCO_TT* tape=dco::tape(y[0]);
  for (int k=0;k<m;k+=ncs) {
      ACDP_LateRecording *e= new ACDP_LateRecording(
        tape,y.size(),y.size(),
        new LateRecording_Target(tape,k,m,ncs,d)
      );
      for (int i=0;i<y.size();i++) {
        e->register_input(y(i));
        e->register_output(y(i));
      }
      e->link_target();
      tape->register_acdp(e);
      e->evaluate_augmented_primal();
  }
}

int main(int c, char* v[]){
  if (c!=4) throw;
  int n=stoi(v[1]), m=stoi(v[2]), ncs=stoi(v[3]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Zero(n));
  const double pi=3.141592653589793;
  for (int i=1;i<n-1;i++) y_indep[i]=sin((2*pi*i)/n);
  double d=1e-2;
  DCO_TT* tape=DCO_TT::create();
  for (int i=0;i<n;i++) tape->register_variable(y_indep[i]);
  y=y_indep;
  euler(m,ncs,d,y);  
  tape->register_output_variable(y[(n-1)/2]);
  dco::derivative(y[(n-1)/2])=1.;
  tape->interpret_adjoint();
  for(int i=0;i<n;i++)
    cout << dco::derivative(y_indep[i]) << endl;
  DCO_TT::remove(tape);
  return 0;
}
