CXX = g++ -std=c++23
CXXFLAGS += -Wall -O
LDLIBS += -lstdc++ -lm -lpthread

all: Out.txt fig_network.svg

fig_network.svg : network.dat Makefile
	echo 'set terminal svg; \
	set output "fig_network.svg"; \
	set title "Wavefunction"; \
	set xlabel "x"; \
	set ylabel "y"; \
	plot \
		"network.dat" using 7:8 with points dt 2 title "Data", \
		"network.dat" using 1:2 with lines dt 1 title "nn(x)", \
		"network.dat" using 1:3 with lines dt 1 title "dnn(x)", \
		"network.dat" using 1:4 with lines dt 1 title "ddnn(x)", \
		"network.dat" using 1:5 with lines dt 1 title "NN(x)", \
		"network.dat" using 1:6 with lines title "g(x)"' |\
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