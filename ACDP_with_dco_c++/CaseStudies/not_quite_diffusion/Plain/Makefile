EIGEN_DIR=$(HOME)/Software/Eigen
DCO_DIR=$(HOME)/Software/dco

TARGETS=$(addsuffix .exe, $(basename $(wildcard *.cpp)))
CPPC=g++
CPPC_FLAGS=-std=c++14 -O3 -Wall -Wextra -pedantic -DDCO_DISABLE_AUTO_WARNING -march=native
DEPS=$(wildcard *.hpp)

INC=-I$(DCO_DIR) -I$(EIGEN_DIR) -I..

all: gt1s.exe ga1s.exe 
	./gt1s.exe 10 100 5 >gt1s.out
	./ga1s.exe 10 100 6 >ga1s.out
	diff gt1s.out ga1s.out

%.exe : %.cpp $(DEPS)
	$(CPPC) $(CPPC_FLAGS) $(INC) $< -o $@

clean :
	rm -f *.exe *.out

.PHONY: all clean
