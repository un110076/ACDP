BASE_DIR=$(HOME)/Documents/gitlab/publications/ACDP/code
DCO_DIR=$(HOME)/Software/dco

TARGETS=$(addsuffix .exe, $(basename $(wildcard *.cpp)))

CPPC=g++ 
CPPC_FLAGS=-std=c++14 -O3 -Wall -Wextra -pedantic -DDCO_DISABLE_AUTO_WARNING

DEPS=$(wildcard ../*.hpp)
INC=-I$(DCO_DIR) -I..

all : $(TARGETS) $(DEPS)
	./gt1s.exe 1000 >gt1s.out
	./ga1s.exe 1000 >ga1s.out
	diff gt1s.out ga1s.out

%.exe : %.cpp $(DEPS)
	$(CPPC) $(CPPC_FLAGS) $(INC) $< -o $@

clean:
	rm -fr *.exe *.out

.PHONY: all clean


