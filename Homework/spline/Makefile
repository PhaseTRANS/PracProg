CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

all: Out.txt fig_linterp.svg fig_qinterp.svg fig_cinterp.svg

fig_cinterp.svg : tabulated.dat cspline.dat Makefile
	echo 'set terminal svg; \
	set output "fig_cinterp.svg"; \
	set title "Spline Interpolation"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"tabulated.dat" using 1:2 with points title "Data", \
		"" smooth csplines with lines dt 3 title "Gnuplot csplines", \
		"cspline.dat" using 1:2 with lines dt 2 title "Cubic Interp.", \
		"cspline.dat" using 1:3 with lines title "Integrated"' |\
	gnuplot

fig_qinterp.svg : tabulated.dat qspline.dat Makefile
	echo 'set terminal svg; \
	set output "fig_qinterp.svg"; \
	set title "Spline Interpolation"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"tabulated.dat" using 1:2 with points title "Data", \
		"qspline.dat" using 1:2 with lines title "Quadratic Interp.", \
		"qspline.dat" using 1:3 with lines title "Integrated"' |\
	gnuplot

fig_linterp.svg : tabulated.dat lspline.dat Makefile
	echo 'set terminal svg; \
	set output "fig_linterp.svg"; \
	set title "Spline Interpolation"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"tabulated.dat" using 1:2 with points title "Data", \
		"lspline.dat" using 1:2 with lines title "Linear Interp.", \
		"lspline.dat" using 1:3 with lines title "Integrated"' |\
	gnuplot

# build output
Out.txt: main
	./$< > $@

main : main.o

main.o : main.cpp matrix.hpp

.PHONY: clean
clean:
	rm Out.txt main main.o tabulated.dat lspline.dat qspline.dat cspline.dat fig_linterp.svg fig_qinterp.svg fig_cinterp.svg
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)