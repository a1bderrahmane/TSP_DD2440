# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# Executable name
TARGET = tsp_solver

# Source files
SOURCES = main.cpp \
          GreedyTSPsolver.cpp \
          LKH.cpp \

# Object files (one .o per .cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# If you want to track headers, you can list them here or just leave it empty
HEADERS =

# Default target
all: $(TARGET)

# Link step
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile step
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run the program
run: $(TARGET)
	./$(TARGET) $(IMPL)

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Rebuild everything
rebuild: clean all

.PHONY: all clean rebuild run
