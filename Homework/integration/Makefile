CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

all: Out.txt fig_erf.svg fig_acc.svg

fig_acc.svg : acc.dat Makefile
	echo 'set terminal svg; \
	set output "fig_acc.svg"; \
	set title "Accuracy plot of erf(z)"; \
	set xlabel "accuracy"; \
	set ylabel "abs(erf(1.0) - exact)"; \
	set logscale xy; \
	plot \
		"acc.dat" using 1:3 with lines dt 1 title "accuracy"' |\
	gnuplot

fig_erf.svg : erf.dat Makefile
	echo 'set terminal svg; \
	set output "fig_erf.svg"; \
	set title "erf(z)"; \
	set xlabel "z"; \
	set ylabel "erf(z)"; \
	plot \
		"erf.dat" using 1:2 with lines dt 1 title "erf(z)"' |\
	gnuplot

# build output
Out.txt: main
	./$< > $@

main : main.o

main.o : main.cpp matrix.hpp

.PHONY: clean
clean:
	rm Out.txt main main.o erf.dat fig_erf.svg
test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)