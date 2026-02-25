# Makefile for Scouting JEC Map Creator
# Requires ROOT installation

# Compiler
CXX = g++

# ROOT configuration
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

# Compiler flags
CXXFLAGS = -O3 -Wall -std=c++17 $(ROOTCFLAGS)
LDFLAGS = $(ROOTLIBS)

# Target executable
TARGET = create_scouting_jec_map

# Source files
SOURCES = create_scouting_jec_map.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $(TARGET)
	@echo "Done! Executable: ./$(TARGET)"

# Compile
%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	@echo "Cleaning..."
	rm -f $(OBJECTS) $(TARGET)
	@echo "Clean done."

# Help
help:
	@echo "Scouting JEC Map Creator - Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make           - Compile the executable"
	@echo "  make clean     - Remove compiled files"
	@echo "  make help      - Show this help message"
	@echo ""
	@echo "Usage after compilation:"
	@echo "  ./create_scouting_jec_map Jet file1.root file2.root ..."
	@echo "  ./create_scouting_jec_map FatJet file1.root file2.root ..."

.PHONY: all clean help
