CXX = g++ -std=c++23
CXXFLAGS += -Wall -Wextra -Werror
LDLIBS += -lstdc++ -lm 

DATA_FILES = erf.dat gamma.dat lngamma.dat
PLOT_FILES = fig_erf.svg fig_gamma.svg fig_lngamma.svg
LOG_FILES = log.fig_erf.gpi log.fig_gamma.gpi log.fig_lngamma.gpi log.main
OBJECTS = main.o
EXECUTABLES = main

# Plot for error function
fig_erf.svg: erf.dat Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set tics out; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot [][-2:2] \
		"$<" index 0 using 1:2 with points pt 7 ps 0.5 title "Tabulated",\
		"$<" index 1 using 1:3 with lines title "erf(x)";' | \
	sed -e 's|\$$@|$@|g' -e 's|\$$<|$<|g' | \
	tee log.fig_erf.gpi | gnuplot

# Plot for gamma function
fig_gamma.svg: gamma.dat Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set tics out; \
	set xlabel "x"; \
	set ylabel "Gamma(x)"; \
	plot [][0:10] \
		"$<" index 0 using 1:2 with points pt 7 ps 0.5 title "Tabulated",\
		"$<" index 1 using 1:3 with lines title "gamma(x)";' | \
	sed -e 's|\$$@|$@|g' -e 's|\$$<|$<|g' | \
	tee log.fig_gamma.gpi | gnuplot

fig_lngamma.svg: lngamma.dat Makefile
	echo 'set terminal svg; \
	set output "$@"; \
	set tics out; \
	set xlabel "x"; \
	set ylabel "log(gamma(x))"; \
	plot [][0:10] \
		"$<" index 0 using 1:2 with points pt 7 ps 0.5 title "Tabulated",\
		"$<" index 1 using 1:3 with lines title "lngamma(x)";' | \
	sed -e 's|\$$@|$@|g' -e 's|\$$<|$<|g' | \
	tee log.fig_lngamma.gpi | gnuplot

erf.dat gamma.dat lngamma.dat: main
	./$< -xmin -5 -xmax 5 -dx 0.5 -dx2 0.01 2>log.main

.PHONY: all plots
all: fig_erf.svg fig_gamma.svg fig_lngamma.svg

plots: all

.PHONY: clean
clean:
clean:
	$(RM) $(EXECUTABLES) $(OBJECTS) $(DATA_FILES) $(PLOT_FILES) $(LOG_FILES) Out.txt
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)