CXX = g++ -std=c++23
CXXFLAGS += -O
LDLIBS += -lstdc++ -lm -lpthread

fig_times.svg: out.times.data Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set datafile separator "	"; \
	set tics in; \
	set xlabel "Matrix size N"; \
	set ylabel "elapsed time"; \
	a=0.1; b=1; \
	f(x)=a*x*x*x+b; \
	fit f(x) "out.times.data" using 1:2 via a,b; \
	plot \
		"$<" index 0 using 1:2 with linespoints title "times", \
		f(x) title "Fit to N^3";' |\
	gnuplot

out.times.data : main
	>$@
	for N in $$(seq 10 2 500); \
	do \
		time --format "$$N\t%e" --output $@ --append \
		./$< $$N 1>out 2>err;\
	done

Out.txt: main
	./$< > $@

main : main.o

main.o : main.cpp matrix.hpp

.PHONY: clean
clean:
	rm out out.times.data log.fig_times.gpi fig_times.svg fit.log err main main.o Out.txt
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)