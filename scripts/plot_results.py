#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# create plots directory if it doesn't exist
os.makedirs('../results/plots', exist_ok=True)

# read timing results from csv file
results_file = "../results/timing_results.csv"

if not os.path.exists(results_file):
    print(f"Results file {results_file} not found.")
    exit(1)

data = pd.read_csv(results_file)

# convert seconds to milliseconds
data['ecoli_cpu'] = data['ecoli_cpu'] * 1000
data['ecoli_gpu_naive'] = data['ecoli_gpu_naive'] * 1000
data['ecoli_gpu_optimized'] = data['ecoli_gpu_optimized'] * 1000
data['yeast_cpu'] = data['yeast_cpu'] * 1000
data['yeast_gpu_naive'] = data['yeast_gpu_naive'] * 1000
data['yeast_gpu_optimized'] = data['yeast_gpu_optimized'] * 1000

# calculate speedup ratios
data['ecoli_gpu_naive_speedup'] = data['ecoli_cpu'] / data['ecoli_gpu_naive']
data['ecoli_gpu_optimized_speedup'] = data['ecoli_cpu'] / data['ecoli_gpu_optimized']
data['ecoli_optimization_speedup'] = data['ecoli_gpu_naive'] / data['ecoli_gpu_optimized']

data['yeast_gpu_naive_speedup'] = data['yeast_cpu'] / data['yeast_gpu_naive']
data['yeast_gpu_optimized_speedup'] = data['yeast_cpu'] / data['yeast_gpu_optimized']
data['yeast_optimization_speedup'] = data['yeast_gpu_naive'] / data['yeast_gpu_optimized']

# focus on k=4, k=8, k=12 only
focus_k_values = [4, 8, 12]
x = np.arange(2)  # 2 datasets
width = 0.25

# 1. runtime comparison - modified for k=4, k=8, k=12 only
plt.figure(figsize=(15, 6))

for i, k in enumerate(focus_k_values):
    plt.subplot(1, 3, i+1)
    
    # find index for this k value
    k_idx = data.index[data['k'] == k].tolist()[0]
    
    # get times for this k value
    cpu_times = [data['ecoli_cpu'].iloc[k_idx], data['yeast_cpu'].iloc[k_idx]]
    gpu_naive_times = [data['ecoli_gpu_naive'].iloc[k_idx], data['yeast_gpu_naive'].iloc[k_idx]]
    gpu_opt_times = [data['ecoli_gpu_optimized'].iloc[k_idx], data['yeast_gpu_optimized'].iloc[k_idx]]
    
    # plot bars
    plt.bar(x - width, cpu_times, width, label='CPU', color='red')
    plt.bar(x, gpu_naive_times, width, label='Naive GPU', color='blue')
    plt.bar(x + width, gpu_opt_times, width, label='Optimized GPU', color='green')
    
    plt.xlabel('Dataset')
    plt.ylabel('Execution Time (ms)')
    plt.title(f'Runtime Comparison (k={k})')
    plt.xticks(x, ['E. coli', 'Yeast'])
    # only add legend to first subplot
    if i == 0:
        plt.legend()
    plt.grid(True, axis='y')
    plt.yscale('log')  # use log scale for better visibility

plt.tight_layout()
plt.savefig('../results/plots/runtime_comparison.png', dpi=200)
plt.close()

# 2. gpu speedup - k=4, k=8, k=12 only
plt.figure(figsize=(15, 6))
width = 0.35  # wider bars for better visibility

for i, k in enumerate(focus_k_values):
    plt.subplot(1, 3, i+1)  # 1 row, 3 columns
    
    # find the correct index for this k value
    k_idx = data.index[data['k'] == k].tolist()[0]
    
    # get speedups for this k value
    naive_speedups = [data['ecoli_gpu_naive_speedup'].iloc[k_idx], data['yeast_gpu_naive_speedup'].iloc[k_idx]]
    opt_speedups = [data['ecoli_gpu_optimized_speedup'].iloc[k_idx], data['yeast_gpu_optimized_speedup'].iloc[k_idx]]
    
    # plot bars
    plt.bar(x - width/2, naive_speedups, width, label='Naive GPU', color='blue')
    plt.bar(x + width/2, opt_speedups, width, label='Optimized GPU', color='green')
    
    # add text labels above bars
    for j, v in enumerate(naive_speedups):
        plt.text(j - width/2, v + 1, f'{v:.1f}x', ha='center', va='bottom')
    
    for j, v in enumerate(opt_speedups):
        plt.text(j + width/2, v + 1, f'{v:.1f}x', ha='center', va='bottom')
    
    plt.xlabel('Dataset')
    plt.ylabel('Speedup Factor (x times)')
    plt.title(f'GPU Speedup over CPU (k={k})')
    plt.xticks(x, ['E. coli', 'Yeast'])
    plt.ylim(0, max(max(naive_speedups), max(opt_speedups)) * 1.2)  # add headroom for labels
    
    # only show legend on first subplot
    if i == 0:
        plt.legend()
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig('../results/plots/gpu_speedup.png', dpi=200)
plt.close()

# 3. performance vs k-mer size - both datasets, 4 subplots (2x2)
# using a smaller figure size and simpler markers
plt.figure(figsize=(12, 10))

# e. coli - regular scale
plt.subplot(2, 2, 1)
plt.plot(data['k'], data['ecoli_cpu'], '-', marker='.', label='CPU', color='red')
plt.plot(data['k'], data['ecoli_gpu_naive'], '-', marker='.', label='GPU Naive', color='blue')
plt.plot(data['k'], data['ecoli_gpu_optimized'], '-', marker='.', label='GPU Optimized', color='green')

plt.xlabel('k-mer size (k)')
plt.ylabel('Execution Time (ms)')
plt.title('E. coli Performance')
plt.legend()
plt.grid(True)

# add vertical lines for k=4, k=8, k=12
for k in focus_k_values:
    plt.axvline(x=k, color='gray', linestyle='--', alpha=0.3)

# e. coli - log scale
plt.subplot(2, 2, 2)
plt.semilogy(data['k'], data['ecoli_cpu'], '-', marker='.', label='CPU', color='red')
plt.semilogy(data['k'], data['ecoli_gpu_naive'], '-', marker='.', label='GPU Naive', color='blue')
plt.semilogy(data['k'], data['ecoli_gpu_optimized'], '-', marker='.', label='GPU Optimized', color='green')

plt.xlabel('k-mer size (k)')
plt.ylabel('Execution Time (ms) - log scale')
plt.title('E. coli Performance (Log Scale)')
plt.legend()
plt.grid(True)

# yeast - regular scale
plt.subplot(2, 2, 3)
plt.plot(data['k'], data['yeast_cpu'], '-', marker='.', label='CPU', color='red')
plt.plot(data['k'], data['yeast_gpu_naive'], '-', marker='.', label='GPU Naive', color='blue')
plt.plot(data['k'], data['yeast_gpu_optimized'], '-', marker='.', label='GPU Optimized', color='green')

plt.xlabel('k-mer size (k)')
plt.ylabel('Execution Time (ms)')
plt.title('Yeast Performance')
plt.legend()
plt.grid(True)

# yeast - log scale
plt.subplot(2, 2, 4)
plt.semilogy(data['k'], data['yeast_cpu'], '-', marker='.', label='CPU', color='red')
plt.semilogy(data['k'], data['yeast_gpu_naive'], '-', marker='.', label='GPU Naive', color='blue')
plt.semilogy(data['k'], data['yeast_gpu_optimized'], '-', marker='.', label='GPU Optimized', color='green')

plt.xlabel('k-mer size (k)')
plt.ylabel('Execution Time (ms) - log scale')
plt.title('Yeast Performance (Log Scale)')
plt.legend()
plt.grid(True)

plt.tight_layout()
# save with lower dpi to reduce file size
plt.savefig('../results/plots/performance_vs_k.png', dpi=150)
plt.close() 