#!/bin/bash

# create results directories if they don't exist
mkdir -p ../results/ecoli
mkdir -p ../results/yeast
mkdir -p ../results/ncu

# define paths to input datasets
ECOLI_DATASET="../data/ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
YEAST_DATASET="../data/ncbi_dataset_saccharomyces_cerevisiae/ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna"

# function to run experiments for a dataset
run_experiments() {
    local dataset=$1
    local dataset_name=$2
    local output_dir="../results/$dataset_name"
    
    echo "Running experiments for $dataset_name..."
    
    # run CPU implementation
    for k in 4 6 8 10 12; do
        echo "Running CPU k-mer counting for $dataset_name with k=$k..."
        ../src/cpu_kmer "$dataset" "$k" > "$output_dir/${dataset_name}_cpu_k${k}.txt"
    done
    
    # run naive GPU implementation
    for k in 4 6 8 10 12; do
        echo "Running naive GPU k-mer counting for $dataset_name with k=$k..."
        ../src/gpu_kmer "$dataset" "$k" > "$output_dir/${dataset_name}_gpu_naive_k${k}.txt"
    done
    
    # run optimized GPU implementation
    for k in 4 6 8 10 12; do
        echo "Running optimized GPU k-mer counting for $dataset_name with k=$k..."
        ../src/gpu_kmer_optimized "$dataset" "$k" > "$output_dir/${dataset_name}_gpu_optimized_k${k}.txt"
    done
}

# run experiments for both datasets
run_experiments "$ECOLI_DATASET" "ecoli"
run_experiments "$YEAST_DATASET" "yeast"

echo "Experiments completed. Results are in the results directory." 