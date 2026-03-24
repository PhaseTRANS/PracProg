CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O2
LDLIBS = -lm

# run: main
# 	./main $(ARGS)

Out.txt: main
	./main > Out.txt

main: main.o sfuns.o
	$(CXX) $^ $(LDLIBS) -o main

main.o: main.cpp sfuns.h
	$(CXX) $(CXXFLAGS) -c main.cpp

sfuns.o: sfuns.cpp sfuns.h
	$(CXX) $(CXXFLAGS) -c sfuns.cpp

.PHONY: clean
clean:
	$(RM) -f main main.o sfuns.o Out.txt

test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)
