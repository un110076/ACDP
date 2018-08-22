EIGEN_DIR=$(HOME)/Software/Eigen
DCO_DIR=$(HOME)/Software/dco

TARGETS=$(addsuffix .exe, $(basename $(wildcard *.cpp)))
CPPC=g++
CPPC_FLAGS=-std=c++14 -O3 -Wall -Wextra -pedantic -march=native
DEPS=$(wildcard *.hpp)

all: f.exe fd.exe 
	./f.exe 10 100 1 >f.out
	./fd.exe 10 100 2 >fd.out
	cd Plain&&$(MAKE);
	cd LateRecording&&$(MAKE);
	cd RecomputeAll&&$(MAKE);
	cd RecursiveBisection&&$(MAKE);
	cd EarlyForwardFiniteDifferences&&$(MAKE);
	cd EarlyTangentPreaccumulation&&$(MAKE);
	cd SymbolicAdjointLSEigenDenseCholesky&&$(MAKE);
	cd SymbolicAdjointNLSEigenSparseLU&&$(MAKE);

%.exe : %.cpp $(DEPS)
	$(CPPC) $(CPPC_FLAGS) -I$(EIGEN_DIR) -I$(DCO_DIR) $< -o $@

clean :
	cd Plain&&$(MAKE) clean;
	cd LateRecording&&$(MAKE) clean;
	cd RecursiveBisection&&$(MAKE) clean;
	cd RecomputeAll&&$(MAKE) clean;
	cd EarlyForwardFiniteDifferences&&$(MAKE) clean;
	cd EarlyTangentPreaccumulation&&$(MAKE) clean;
	cd SymbolicAdjointLSEigenDenseCholesky&&$(MAKE) clean;
	cd SymbolicAdjointNLSEigenSparseLU&&$(MAKE) clean;
	rm -f *.exe *.out *.ps *.pdf

.PHONY: all clean print
