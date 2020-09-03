CPPC=g++ 
CPPC_FLAGS=-std=c++14 -O3 -march=native -Wall -Wextra -pedantic -DDCO_DISABLE_AUTO_WARNING
DEPS=../implicit_euler.hpp
BASE_DIR=../../..
EIGEN_DIR=$(HOME)/Software/Eigen
DCO_DIR=$(HOME)/Software/dco
PATTERN_DIR=$(BASE_DIR)/SymbolicAdjointNLSEigenSparseLU

all : main.exe 
	./main.exe 1000 10 >main.out
	diff main.out ../Plain/ga1s.out

%.exe : %.cpp $(DEPS)
	$(CPPC) $(CPPC_FLAGS) -I$(DCO_DIR) -I$(EIGEN_DIR) -I$(PATTERN_DIR) $< -o $@

clean:
	rm -fr main.exe main.out

.PHONY: all clean


