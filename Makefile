# -*- Makefile -*-

## sherlock
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall

all:main.out
.PHONY: all

main.out:main.cpp src/*.hpp
	$(CXX) $(CXXFLAGS)	main.cpp -o main.out
clean:
	@echo "clean up..."
	rm -f *.o
	rm -f *.out
.PHONY: clean