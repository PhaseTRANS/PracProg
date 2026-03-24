CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm

Out.txt : main
	./main > Out.txt

main : main.o
	$(CXX) main.o -o main $(LDLIBS)

main.o : main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

.PHONEY:clean
clean:
	$(RM) --force main main.o Out.txt

test:
	@echo $(CXX)
	@echo $(CXX)
	@echo $(RM)