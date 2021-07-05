# -*- Makefile -*-

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
BUILD_DIR = build
PROGRAM = $(BUILD_DIR)/diamondMC
TEST = $(BUILD_DIR)/test
#all:main.out, MD.out, test.out
all:$(BUILD_DIR) $(PROGRAM)
.PHONY: all

$(BUILD_DIR):
	mkdir -p $@

SRC_FILES = $(wildcard src/*.cpp utils/*.hpp)
SRC_HEADERS = $(wildcard src/*.hpp utils/*.cpp)
$(PROGRAM):main.cpp $(SRC_FILES) $(SRC_HEADERS)
	@echo "Building  : " $@
	$(CXX) $(CXXFLAGS) -o $@ $^

$(TEST):main.cpp $(SRC_FILES) $(SRC_HEADERS)
	@echo "Building  : " $@
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	@echo "clean up..."
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BUILD_DIR)/*.out
	rm -f $(PROGRAM)
	rm -f $(TEST)
.PHONY: clean
