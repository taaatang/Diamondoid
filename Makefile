# -*- Makefile -*-

CXX = g++-13
CXXFLAGS = -std=c++20 -O3 -Wall
BUILD_DIR = build
PROGRAM = $(BUILD_DIR)/diamondMC

all:$(BUILD_DIR) $(PROGRAM)
.PHONY: all

$(BUILD_DIR):
	mkdir -p $@

SRC_FILES = $(wildcard src/*.cpp)
SRC_HEADERS = $(wildcard src/*.hpp)
$(PROGRAM):main.cpp $(SRC_FILES) $(SRC_HEADERS)
	@echo "Building  : " $@
	$(CXX) $(CXXFLAGS) -o $@ main.cpp $(SRC_FILES)

clean:
	@echo "clean up..."
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BUILD_DIR)/*.out
	rm -f $(PROGRAM)
.PHONY: clean
