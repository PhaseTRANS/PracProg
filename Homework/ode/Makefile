CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

all: Out.txt fig_ode_result.svg fig_orbits.svg

fig_orbits.svg : orbits.dat Makefile
	echo 'set terminal svg; \
	set output "fig_orbits.svg"; \
	set title "Orbit solutions"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"orbits.dat" index 0 using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with lines title "eps=0.0, u'\''=0.0",  \
		"orbits.dat" index 1 using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with lines title "eps=0.0, u'\''=-0.5", \
		"orbits.dat" index 2 using (1/$$2)*cos($$1):(1/$$2)*sin($$1) with lines dt 2 title "eps=0.01, u'\''=-0.5"' \
	| gnuplot

fig_ode_result.svg : ode_result.dat Makefile
	echo 'set terminal svg; \
	set output "fig_ode_result.svg"; \
	set title "ODE solution"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"ode_result.dat" using 1:2 with lines dt 1 title "theta(t)", \
		"ode_result.dat" using 1:3 with lines dt 1 title "omega(t)"' |\
	gnuplot

# build output
Out.txt: main
	./$< > $@

main : main.o

main.o : main.cpp matrix.hpp

.PHONY: clean
clean:
	rm Out.txt main main.o ode_result.dat fig_ode_result.svg fig_orbits.svg
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)