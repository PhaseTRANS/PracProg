CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

all: Out.txt fig_wave.svg

fig_wave.svg : wave.dat Makefile
	echo 'set terminal svg; \
	set output "fig_wave.svg"; \
	set title "Wavefunction"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"wave.dat" using 1:2 with lines dt 1 title "Numerical f(r)", \
		"wave.dat" using 1:3 with lines dt 2 title "Exact f(r)"' |\
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