all : 
	cd Plain && $(MAKE)
	cd LateRecording && $(MAKE)
	cd StaticEvolutionRecomputeAll && $(MAKE)
	cd StaticEvolutionRecursiveBisection && $(MAKE)
	doxygen

clean:
	cd Plain && $(MAKE) clean
	cd LateRecording && $(MAKE) clean
	cd StaticEvolutionRecomputeAll && $(MAKE) clean
	cd StaticEvolutionRecursiveBisection && $(MAKE) clean
	rm -fr html

.PHONY: all clean


