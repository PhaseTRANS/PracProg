CXX = g++ -std=c++23
CXXFLAGS += -Wall -Wextra -Werror
LDLIBS += -lstdc++ -lm 

Out.txt: main
	./$< > $@

main: main.o

.PHONY: clean
clean:
	$(RM) -f main main.o Out.txt

test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)
