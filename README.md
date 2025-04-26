# CUDA K-mer Counting

this project implements k-mer counting using CUDA for GPU acceleration. it compares the performance of CPU and GPU implementations on different datasets.

## project structure

```
.
├── data/                          # input datasets
│   ├── ncbi_dataset/              # e. coli genome
│   └── ncbi_dataset_saccharomyces_cerevisiae/  # yeast genome
├── src/                           # source code
│   ├── cpu_kmer.c                 # CPU implementation
│   ├── gpu_kmer.cu                # naive GPU implementation
│   └── gpu_kmer_optimized.cu      # optimized GPU implementation
├── scripts/                       # utility scripts
│   ├── run_experiments.sh         # runs experiments and generates results
│   ├── run_ncu.sh                 # runs NVIDIA Nsight Compute profiling
│   ├── generate_timing_csv.py     # extracts timing data into CSV
│   └── plot_results.py            # generates performance plots
└── results/                       # output directory for results
    ├── ecoli/                     # e. coli dataset results
    ├── yeast/                     # yeast dataset results 
    └── plots/                     # performance visualization plots
```

## building the project

to build the project, run:

```bash
cd src
make
```

this will compile all implementations:
- `cpu_kmer`: CPU version
- `gpu_kmer`: naive GPU version
- `gpu_kmer_optimized`: optimized GPU version

## running experiments

to run experiments and generate results:

```bash
cd scripts
./run_experiments.sh
```

this script will:
1. run all implementations on both datasets for k=4,6,8,10,12
2. save output files to the dataset-specific directories (ecoli/ and yeast/)
3. collect timing information for performance analysis

## performance analysis

the project includes several scripts for analyzing and visualizing performance:

```bash
# extract timing data into a CSV file
cd scripts
./generate_timing_csv.py

# generate performance plots
./plot_results.py

# run NVIDIA profiling (requires NVIDIA Nsight Compute)
./run_ncu.sh
```

### performance visualization

the `plot_results.py` script generates several plots to help understand performance:

1. **runtime comparison**: compares execution times between CPU, naive GPU, and optimized GPU implementations for k=4, k=8, and k=12.

2. **gpu speedup**: visualizes the speedup factor (how many times faster the GPU is compared to CPU) for k=4, k=8, and k=12.

3. **performance vs k-mer length**: shows how performance scales with k-mer length for both datasets:
   - separate graphs for e. coli and yeast datasets
   - each graph shows both regular and logarithmic scales
   - CPU, naive GPU, and optimized GPU implementations are compared side-by-side

all plots are saved to the `results/plots/` directory.

## GPU profiling

for detailed GPU performance analysis, the project includes `run_ncu.sh` script that uses NVIDIA Nsight Compute to collect low-level metrics:

- memory throughput
- occupancy
- active threads per warp
- cache hit rates (L1 and L2)
- instructions per cycle (IPC)

these metrics help identify performance bottlenecks and opportunities for optimization.

## requirements

- CUDA toolkit
- gcc
- python3 with matplotlib and pandas
- make
- NVIDIA Nsight Compute (for profiling)

## datasets

the project uses two datasets:
1. e. coli genome (4.6 million bases)
2. yeast genome (12.1 million bases)

these datasets are used to evaluate how the implementations scale with different sequence lengths.
