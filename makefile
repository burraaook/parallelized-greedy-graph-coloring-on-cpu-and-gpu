# Compiler and flags
CXX = g++
NVCC = nvcc
CXXFLAGS = -std=c++17 -g #-Wall -Wextra -Werror -Wpedantic
NVCCFLAGS = -g

# Source files
CPUSRC = cpu/src/graph.cpp cpu/src/bitcolor_cpu.cpp
CUFILES = gpu/bitcolor_gpu.cu

# Object files
CPUOBJS = $(CPUSRC:cpu/src/%.cpp=cpu/obj/%.o)
CUOBJS = $(CUFILES:gpu/%.cu=gpu/%.o)

# Headers
CPPHDR = cpu/lib/graph.hpp cpu/lib/dct.hpp

# Targets
CPUTRGT = bitcolor_cpu
GPUTRGT = bitcolor_gpu

all: cpu_dirs gpu_dirs $(CPUTRGT) $(GPUTRGT)

# Directories
cpu_dirs:
	mkdir -p cpu/obj

gpu_dirs:
	mkdir -p gpu/output

# CPU build rules
cpu/obj/%.o: cpu/src/%.cpp $(CPPHDR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(CPUTRGT): $(CPUOBJS)
	$(CXX) $(CXXFLAGS) -o $(CPUTRGT) $(CPUOBJS)

# GPU build rules
gpu/%.o: gpu/%.cu
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

$(GPUTRGT): $(CUOBJS)
	$(NVCC) $(NVCCFLAGS) -o $(GPUTRGT) $(CUOBJS)

# Clean
clean:
	rm -f cpu/obj/*.o gpu/*.o
	rm -f $(CPUTRGT) $(GPUTRGT)
	rm -rf cpu/obj

.PHONY: all clean cpu_dirs gpu_dirs
