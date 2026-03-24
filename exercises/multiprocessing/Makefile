CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

DATA_FILES = erf.dat 
PLOT_FILES = fig_erf.svg 
LOG_FILES = log.fig_erf.gpi
OBJECTS = main.o
EXECUTABLES = main

N=1e9
TIME = time --append --output $@ --format "$$nthreads %e %U"
out.times : main Makefile
	>$@
	>log.main
	for nthreads in 1 2 3 4 5 6 7 8 9 10; \
		do $(TIME) ./main -terms $(N) -threads $$nthreads >> log.main;\
	done

fig_times.svg: out.times Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set tics in; \
	set xlabel "threads"; \
	set ylabel "elapsed time"; \
	a=1; b=10;\
	f(x)=a+b/x;\
	fit f(x) "out.times" using 1:2 via a,b; \
	plot \
		"$<" index 0 using 1:2 with linespoints title "times", \
		f(x) title "fit to 1/n";' |\
	tee log.fig_times.gpi | gnuplot

.PHONY: all plots
all: fig_times.svg

plots: all

.PHONY: clean
clean:
	$(RM) $(EXECUTABLES) $(OBJECTS) $(DATA_FILES) $(PLOT_FILES) $(LOG_FILES) Out.txt
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)