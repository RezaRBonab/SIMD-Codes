# Compiler
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -msse4.1

# Source files
SRCS = $(wildcard *.cpp)
HDRS = $(wildcard *.h *.hpp)

# Build directory
BUILD_DIR = build

# Object files
OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)

# Executable name
TARGET = main

# Build target
all: $(BUILD_DIR) $(TARGET)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJS) $(TARGET)
	rm -rf $(BUILD_DIR)

