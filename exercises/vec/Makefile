CXX = g++ -std=c++23
CXXFLAGS += -Wall -Wextra -Werror
LDLIBS += -lstdc++ -lm 

run : Out.txt
	cat Out.txt

Out.txt: main
	./$< > $@

main: main.o

main.o : main.cpp vec.h

.PHONY: clean
clean:
	$(RM) -f main main.o Out.txt

test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)
