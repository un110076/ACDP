CPPC=g++ 
CPPC_FLAGS=-std=c++14 -O3 -march=native -Wall -Wextra -pedantic 
DEPS=$(wildcard, *.hpp)

all : main.exe 
	./main.exe 10

%.exe : %.cpp $(DEPS)
	$(CPPC) $(CPPC_FLAGS) $< -o $@

clean:
	rm -fr main.exe 

.PHONY: all clean


