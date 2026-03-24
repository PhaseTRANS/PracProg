CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

all: Out.txt fig_fit.svg 

# plot fit and data
fig_fit.svg : decay.dat fit.dat Makefile
	echo 'set terminal svg; \
	set output "fig_fit.svg"; \
	set title "Least Squares Fit"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"decay.dat" using 1:2:3 with yerrorbars title "data", \
		"fit.dat" using 1:2 with lines title "fit", \
		"fit.dat" using 1:3:4 with filledcurves fill transparent solid 0.3 title "uncertainty"' |\
	gnuplot

# build output
Out.txt: main
	./$< > $@

main : main.o

main.o : main.cpp matrix.hpp

.PHONY: clean
clean:
	rm Out.txt main main.o decay.dat fit.dat fig_fit.svg 
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)