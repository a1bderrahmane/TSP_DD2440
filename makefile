# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# Executable name
TARGET = tsp_solver

# Source files
SOURCES = main.cpp \
          GreedyTSPsolver.cpp \
          Christofides.cpp \
          Implementation2.cpp \
          Implementation3.cpp \
          Implementation4.cpp

# Header files (for dependency tracking)
HEADERS = ITSPsolver.hpp \
          GreedyTSPsolver.hpp \
          Christofides.hpp \
          Implementation2.hpp \
          Implementation3.hpp \
          Implementation4.hpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

IMPL ?= 0

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile each .cpp file into .o
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the program with a given implementation (default: 0)
run: $(TARGET)
	./$(TARGET) $(IMPL)

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Rebuild everything
rebuild: clean all

.PHONY: all clean rebuild run