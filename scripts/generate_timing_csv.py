#!/usr/bin/env python3

import os
import re
import pandas as pd

def parse_timing(file_path):
    """parse timing from output file"""
    with open(file_path, 'r') as f:
        content = f.read()
        # extract timing using regex
        match = re.search(r'CPU time: ([\d.]+) seconds|GPU time: ([\d.]+) seconds', content)
        if match:
            return float(match.group(1) or match.group(2))
    return None

def main():
    # initialize data structure
    data = {
        'k': [],
        'ecoli_cpu': [],
        'ecoli_gpu_naive': [],
        'ecoli_gpu_optimized': [],
        'yeast_cpu': [],
        'yeast_gpu_naive': [],
        'yeast_gpu_optimized': []
    }
    
    # process e. coli results
    ecoli_dir = '../results/ecoli'
    for k in [4, 6, 8, 10, 12]:
        data['k'].append(k)
        
        # cpu timing
        cpu_file = os.path.join(ecoli_dir, f'ecoli_cpu_k{k}.txt')
        data['ecoli_cpu'].append(parse_timing(cpu_file))
        
        # gpu naive timing
        gpu_naive_file = os.path.join(ecoli_dir, f'ecoli_gpu_naive_k{k}.txt')
        data['ecoli_gpu_naive'].append(parse_timing(gpu_naive_file))
        
        # gpu optimized timing
        gpu_opt_file = os.path.join(ecoli_dir, f'ecoli_gpu_optimized_k{k}.txt')
        data['ecoli_gpu_optimized'].append(parse_timing(gpu_opt_file))
    
    # process yeast results
    yeast_dir = '../results/yeast'
    for k in [4, 6, 8, 10, 12]:
        # cpu timing
        cpu_file = os.path.join(yeast_dir, f'yeast_cpu_k{k}.txt')
        data['yeast_cpu'].append(parse_timing(cpu_file))
        
        # gpu naive timing
        gpu_naive_file = os.path.join(yeast_dir, f'yeast_gpu_naive_k{k}.txt')
        data['yeast_gpu_naive'].append(parse_timing(gpu_naive_file))
        
        # gpu optimized timing
        gpu_opt_file = os.path.join(yeast_dir, f'yeast_gpu_optimized_k{k}.txt')
        data['yeast_gpu_optimized'].append(parse_timing(gpu_opt_file))
    
    # create dataframe and save to csv
    df = pd.DataFrame(data)
    os.makedirs('../results', exist_ok=True)
    df.to_csv('../results/timing_results.csv', index=False)
    print("timing results saved to ../results/timing_results.csv")

if __name__ == '__main__':
    main() 