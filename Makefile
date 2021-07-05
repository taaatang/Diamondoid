# -*- Makefile -*-

## sherlock
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall

all:main.out, MD.out, test.out
.PHONY: all

main.out:main.cpp src/*.cpp src/*.hpp utils/*.hpp
	$(CXX) $(CXXFLAGS)	main.cpp src/*.cpp -o build/main.out

MD.out:mainMD.cpp srcMD/*.cpp srcMD/*.hpp utils/*.hpp
	$(CXX) $(CXXFLAGS) mainMD.cpp srcMD/*.cpp -o build/MD.out

test.out:test.cpp src/*.cpp src/*.hpp utils/*.hpp
	$(CXX) $(CXXFLAGS)	test.cpp src/*.cpp -o build/test.out
clean:
	@echo "clean up..."
	rm -rf build/*.o
	rm -rf build/*.out
.PHONY: clean
