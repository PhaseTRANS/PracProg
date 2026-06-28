CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra -pedantic

all: spline

spline: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o spline

clean:
	rm -f spline out.svg *.o