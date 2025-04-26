CC = g++
NVCC = nvcc
CFLAGS = -O3 -std=c++14
NVCCFLAGS = -O3 -std=c++14

all: cpu_kmer gpu_kmer gpu_kmer_optimized

cpu_kmer: cpu_kmer.cpp
	$(CC) $(CFLAGS) -o cpu_kmer cpu_kmer.cpp

gpu_kmer: gpu_kmer.cu
	$(NVCC) $(NVCCFLAGS) -o gpu_kmer gpu_kmer.cu

gpu_kmer_optimized: gpu_kmer_optimized.cu
	$(NVCC) $(NVCCFLAGS) -o gpu_kmer_optimized gpu_kmer_optimized.cu

clean:
	rm -f cpu_kmer gpu_kmer gpu_kmer_optimized
