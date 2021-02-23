# -*- Makefile -*-

## sherlock
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall

all:main.out
.PHONY: all

main.out:main.cpp src/*.hpp utils/*.hpp
	$(CXX) $(CXXFLAGS)	main.cpp -o build/main.out

MD.out:mainMD.cpp srcMD/*.cpp srcMD/*.hpp utils/*.hpp
	$(CXX) $(CXXFLAGS) mainMD.cpp srcMD/*.cpp -o build/MD.out

clean:
	@echo "clean up..."
	rm -rf build/*.o
	rm -rf build/*.out
.PHONY: clean