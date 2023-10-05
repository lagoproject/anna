################################################################################
# File:        Makefile
# Description: Builds the ANNA toolsuite, a scientific tool for analyzing 
#              LAGO astrophysical data. Provides targets for building applications 
#              (dump, example, raw, sol, ...) and cleaning build artifacts.
# Author:      Hernán Asorey
# Email:       asoreyh@gmail.com
# Date:        2012-2023
# 
# Copyright:   [2012] [The LAGO Collaboration]
# License:     BSD-3-Clause
# See the LICENSE file in the project root for full license information.
################################################################################

# Compiler and compile options
CXX = g++
CXXFLAGS = -Wall -Iinclude

# Source and build directories
SRC_DIR = src
BUILD_DIR = build

# Tests
TESTS = check-env

# Applications
APPS = dump example raw sol

# Default target
all: $(TESTS) $(APPS)

# Compile each application
$(APPS): %: $(BUILD_DIR)/%.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule for compiling object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(APPS)

.PHONY: all clean

check-env:
ifndef LAGO_ANNA
$(warning LAGO_ANNA is not set. I will define it to $(PWD), and modify the .bashrc)
$(shell $(PWD)/lago-anna.sh)
LAGO_ANNA=$(PWD)
else
$(info Environment variable LAGO_ANNA is set to $(LAGO_ANNA))
endif