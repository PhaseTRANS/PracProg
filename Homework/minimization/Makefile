CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

all: Out.txt fig_higgs.svg

fig_higgs.svg : higgs_fit.dat Makefile
	echo 'set terminal svg; \
	set output "fig_higgs.svg"; \
	set title "Wavefunction"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"higgs.dat" using 1:2 with points title "Data", \
		"higgs_fit.dat" using 1:4 with lines dt 1 title "Fit"' |\
	gnuplot

# build output
Out.txt: main
	./$< > $@

main : main.o

main.o : main.cpp matrix.hpp

.PHONY: clean
clean:
	rm Out.txt main main.o
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)