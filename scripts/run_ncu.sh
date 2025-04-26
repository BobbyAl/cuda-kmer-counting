#!/bin/bash

# create results directories if they don't exist
mkdir -p ../results/ecoli
mkdir -p ../results/yeast
mkdir -p ../results/ncu

# define paths to input datasets
ECOLI_DATASET="../data/ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
YEAST_DATASET="../data/ncbi_dataset_saccharomyces_cerevisiae/ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna"

# function to run ncu on a kernel
run_ncu() {
    local dataset=$1
    local dataset_name=$2
    local k=$3
    local kernel=$4
    
    echo "Running ncu on $kernel for $dataset_name with k=$k..."
    
    # run ncu with default metrics
    ncu --set full ../src/$kernel "$dataset" "$k" > "../results/ncu/${dataset_name}_${kernel}_k${k}.txt"
}

# run ncu on both datasets for k=4,8,12
for k in 4 8 12; do
    # run on naive GPU implementation
    run_ncu "$ECOLI_DATASET" "ecoli" "$k" "gpu_kmer"
    run_ncu "$YEAST_DATASET" "yeast" "$k" "gpu_kmer"
    
    # run on optimized GPU implementation
    run_ncu "$ECOLI_DATASET" "ecoli" "$k" "gpu_kmer_optimized"
    run_ncu "$YEAST_DATASET" "yeast" "$k" "gpu_kmer_optimized"
done

echo "NCU analysis completed. Results are in the results/ncu directory." 