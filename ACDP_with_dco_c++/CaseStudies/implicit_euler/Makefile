TARGETS=$(addsuffix .exe, $(basename $(wildcard *.cpp)))
CPPC=g++ 
CPPC_FLAGS=-std=c++14 -O3 -Wall -Wextra -pedantic -march=native
DEPS=$(wildcard *.hpp)

all : $(TARGETS) $(DEPS)
	./f.exe 1000 >f.out
	./fd.exe 1000 >fd.out
	cd Plain&&$(MAKE)
	cd LateRecording&&$(MAKE)
	cd EarlyTangentPreaccumulation&&$(MAKE)
	cd EarlyForwardFiniteDifferences&&$(MAKE)
	cd SymbolicAdjointNLSEigenSparseLU&&$(MAKE)
	cd RecomputeAll&&$(MAKE)
	cd RecursiveBisection&&$(MAKE)

%.exe : %.cpp $(DEPS)
	$(CPPC) $(CPPC_FLAGS) $< -o $@

clean:
	rm -fr *.exe *.out
	cd Plain&&$(MAKE) clean
	cd LateRecording&&$(MAKE) clean
	cd EarlyTangentPreaccumulation&&$(MAKE) clean
	cd EarlyForwardFiniteDifferences&&$(MAKE) clean
	cd SymbolicAdjointNLSEigenSparseLU&&$(MAKE) clean
	cd RecomputeAll&&$(MAKE) clean
	cd RecursiveBisection&&$(MAKE) clean

.PHONY: all clean


