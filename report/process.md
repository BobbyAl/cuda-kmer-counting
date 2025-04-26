# CUDA K-mer Counting Project

## Data
The dataset was sourced from The National Center for Biotechnology Information (The National Library of Medicine: https://www.ncbi.nlm.nih.gov/). The direct link for the dataset is https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/

There are two files from this NCBI dataset: GCF_000005845.2 & GCA_000005845.2.
GCA (GenBank Assembly) refers to the assembly sumbitted by GenBank, which has raw submissions from researchers, that's a little less curated, but more frequently updated. 

GCF (NCBI RefSeq Assembly) refers to the reference sequences that are curated, manually reviewed and considered the gold standard reference. 

For this project, we will use the GCF dataset due it's higher quality and non-redundancy. 

## Parsing
There wasn't much need for parsing, as the dataset provided only had a header, which we skipped past. 

## K-mer Counting Implementations

### CPU Implementation
- Uses a simple hash table (std::unordered_map) to count k-mers
- Processes the sequence sequentially

### GPU Implementations

#### Naive GPU Version
- Uses global memory for kmer counting
- Processes sequence in chunks to handle memory constraints
- Uses atomic operations for counting
- Features:
  - Chunked processing (1 million bases per chunk)
  - 2D grid structure (x-dimension for parallel processing, y-dimension for chunks)

#### Optimized GPU Version
- Uses shared memory for frequently accessed kmers
- Implements memory coalescing
- More efficient memory access patterns
- Features:
  - Shared memory optimization
  - Better thread utilization
  - Reduced global memory access

## Performance Characteristics

Key observations:
1. GPU implementations are significantly faster than CPU
2. Optimized GPU version consistently outperforms naive version
3. Performance gap increases with larger k values
4. Memory usage becomes a limiting factor for naive GPU with k ≥ 12


## Python Visualization and Analysis

The project includes Python scripts for visualizing and analyzing the performance results. The main visualization script (`plot_results.py`) generates plots to help understand the performance characteristics of the GPU implementation compared to the CPU baseline.

### Metrics

1. **Execution Time**: The raw time taken by each implementation to count k-mers, measured in seconds. This is plotted on both linear and logarithmic scales to show the full range of performance differences.

2. **Speedup Factor**: This is calculated as the ratio of CPU execution time to GPU execution time. For example, if the CPU takes 1 second and the GPU takes 0.1 seconds, the speedup factor is 10x. This metric shows how much faster the GPU implementation is compared to the CPU baseline.

### Visualization Plots

1. **Performance Comparison (2 subplots)**:
   - Execution times for both CPU and GPU implementations across different k-mer sizes
   - Direct comparison between E. coli and yeast datasets
   - Speedup factors showing GPU performance improvement

2. **Log Scale Performance Comparison**:
   - Shows the same data as the main comparison but on a logarithmic scale
   - Particularly useful for visualizing large differences in performance
   - Makes it easier to see the relationship between k-mer size and execution time