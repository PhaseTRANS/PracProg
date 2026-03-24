CXX = g++ -std=c++23
CXXFLAGS += -Wall -Wextra -Werror
LDLIBS += -lstdc++ -lm 

# run : Out.txt
# 	cat Out.txt

Out.txt: main
# 	echo 1 2 3 4 5 | ./$< > $@  # for part 2
# 	./$< > $@
# 	./main --input infile.txt --output outfile.txt


main: main.o

main.o : main.cpp

.PHONY: clean
clean:
	$(RM) -f main main.o Out.txt

test:
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(RM)
