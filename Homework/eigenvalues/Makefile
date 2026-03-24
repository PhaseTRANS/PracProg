CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

all: Out.txt fig_dr.svg fig_rmax.svg fig_waves.svg fig_diag_times.svg wavefunctions.dat diag_times.dat

# plot diagolaization times
fig_diag_times.svg : diag_times.dat Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set tics in; \
	set xlabel "Matrix size N"; \
	set ylabel "elapsed time"; \
	a=0.1; b=1; \
	f(x)=a*x*x*x+b; \
	fit f(x) "diag_times.dat" using 1:2 via a,b; \
	plot \
		"$<" index 0 using 1:2 with linespoints title "times", \
		f(x) title "Fit to N^3";' |\
	gnuplot

# plot wavefunctions
fig_waves.svg : wavefunctions.dat Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set tics in; \
	set xlabel "r"; \
	set ylabel "|Psi|^2"; \
	plot \
		"$<" index 0 using 1:2 with linespoints title "k=0", \
		"$<" index 0 using 1:3 with linespoints title "k=1", \
		"$<" index 0 using 1:4 with linespoints title "k=2";' |\
	gnuplot

# plot change of dr parameter
fig_dr.svg : out.dr.data Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set tics in; \
	set xlabel "dr"; \
	set ylabel "E0"; \
	plot \
		"$<" index 0 using 2:3 with linespoints title "dr";' |\
	gnuplot

# plot change of rmax parameter
fig_rmax.svg : out.rmax.data Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set tics in; \
	set xlabel "rmax"; \
	set ylabel "E0"; \
	plot \
		"$<" index 0 using 1:3 with linespoints title "rmax";' |\
	gnuplot

# build wavefuntion dataset
wavefunctions.dat: main
	./main -rmax 10 -dr 0.05 -wf wavefunctions.dat

# build diagonalization times dataset
diag_times.dat : main
	./main -diagf diag_times.dat

# build dr dataset
out.dr.data: main
	for r in $$(seq 0.1 0.1 1); do \
		./$< -rmax 10 -dr $$r 2>> $@; \
	done

# build rmax dataset
out.rmax.data: main
	for r in $$(seq 1 10); do \
		./$< -rmax $$r -dr 0.3 2>> $@; \
	done

# build output for Part A
Out.txt: main
	./$< > $@
	./$< -rmax 10 -dr 0.3 > $@

main : main.o

main.o : main.cpp matrix.hpp

.PHONY: clean
clean:
	rm Out.txt main main.o out.rmax.data out.dr.data wavefunctions.dat diag_times.dat fig_waves.svg fig_diag_times.svg fig_dr.svg fig_rmax.svg fit.log 
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)