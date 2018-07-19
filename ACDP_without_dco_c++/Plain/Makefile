CPPC=g++ 
CPPC_FLAGS=-std=c++14 -O3 -march=native -Wall -Wextra -pedantic 
DEPS=$(wildcard, ../ACDP/*.hpp)
INCS=-I../ACDP/

all : main.exe 
	./main.exe 10

%.exe : %.cpp $(DEPS)
	$(CPPC) $(CPPC_FLAGS) $(INCS) $< -o $@

clean:
	rm -fr main.exe

.PHONY: all clean


