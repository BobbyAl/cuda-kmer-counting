#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>

// max kmer size
#define MAX_K 32

// calculate 4^k (max possible kmers for dna)
#define MAX_KMERS(k) (1ULL << (2 * (k)))

// max numnber of unique kmers to store in shared memory 
#define SHARED_MEM_SIZE 1024

// invalid kmer index
#define INVALID_KMER 0xFFFFFFFFFFFFFFFFULL

// convert kmer string to an int index (A=0, C=1, G=2, T=3)
__device__ __host__ unsigned long long kmer_to_index(const char* kmer, int k) 
{
    unsigned long long index = 0;

    for (int i = 0; i < k; i++) 
    {
        unsigned int nucleotide;
        switch (kmer[i]) 
        {
            case 'A': nucleotide = 0; 
                break;
            case 'C': nucleotide = 1; 
                break;
            case 'G': nucleotide = 2; 
                break;
            case 'T': nucleotide = 3; 
                break;
            default: 
                return INVALID_KMER;
        }
        index = (index << 2) | nucleotide;
    }
    return index;
}

// convert index back to kmer string
void index_to_kmer(unsigned long long index, int k, char* kmer) 
{
    char nucleotides[] = {'A', 'C', 'G', 'T'};
    
    for (int i = k - 1; i >= 0; i--) 
    {
        kmer[i] = nucleotides[index & 3];
        index >>= 2;
    }
    kmer[k] = '\0';
}

// function to read fatsa file and extract dna sequence
std::string read_fatsa_file(const std::string& filename) 
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return "";
    }

    std::string line;
    std::string sequence;
    
    // skipping header line
    std::getline(file, line);
    
    // reading sequence lines
    while (std::getline(file, line)) 
            sequence += line;
    
    file.close();
    return sequence;
}

// function to print top N most frequent kmers
void print_top_kmers(const unsigned int* kmer_counts, int k, int n) 
{
    // create vector of pairs for sorting
    std::vector<std::pair<unsigned long long, unsigned int>> kmer_vector;
    
    // getting only non-zero counts
    unsigned long long max_kmers = MAX_KMERS(k);

    for (unsigned long long i = 0; i < max_kmers; i++) 
    {
        if (kmer_counts[i] > 0)
            kmer_vector.push_back({i, kmer_counts[i]});
    }
    
    // sorting by frequency (descending)
    auto compare = [](const std::pair<unsigned long long, unsigned int>& kmer1, 
                     const std::pair<unsigned long long, unsigned int>& kmer2) {
        return kmer1.second > kmer2.second;
    };

    std::sort(kmer_vector.begin(), kmer_vector.end(), compare);
    
    // printing top N kmers
    std::cout << "Top " << n << " most frequent k-mers:" << std::endl;
    char kmerStr[MAX_K + 1];
    int count = 0;
    for (const auto& pair : kmer_vector) 
    {
        if (count >= n) break;
        index_to_kmer(pair.first, k, kmerStr);
        std::cout << kmerStr << ": " << pair.second << std::endl;
        count++;
    }
}

// kernel to count kmers with coalesced memory access
__global__ void count_kmers_coalesced_kernel(const char* sequence, int sequence_len, int k, unsigned int* global_kmer_counts, int coalesce_size)
{
    int tid = threadIdx.x;
    int bid = blockIdx.x;
    int block_size = blockDim.x;
    
    // each block processes a chunk of the sequence that is coalesce_size in length
    int start_pos = bid * coalesce_size;
    
    // create a local cache in registers for the first part of the sequence this thread will process
    char local_seq[MAX_K];

    if (start_pos + tid < sequence_len) 
    {
        // coalesced memory access
        for (int i = 0; i < k && start_pos + tid + i < sequence_len; i++)
            local_seq[i] = sequence[start_pos + tid + i];
    }
    
    // shared memory for local k-mer counts
    extern __shared__ unsigned int shmem[];
    
    // clearing shared memory
    for (int i = tid; i < SHARED_MEM_SIZE; i += block_size) 
        shmem[i] = 0;

    __syncthreads();
    
    // each thread processes multiple positions based on its lane id
    for (int i = 0; i < coalesce_size && start_pos + i + k - 1 < sequence_len; i++) 
    {
        if (tid + i < coalesce_size) 
        {
            // building kmer from the sequence
            char kmer[MAX_K];
            
            // using cached data for the first part if possible
            int cached_chars = min(k, MAX_K);
            for (int j = 0; j < cached_chars && tid + j < MAX_K; j++) 
                kmer[j] = local_seq[j];
            
            // reading any remaining characters needed
            for (int j = cached_chars; j < k; j++) 
            {
                if (start_pos + tid + i + j < sequence_len)
                    kmer[j] = sequence[start_pos + tid + i + j];
            }
            
            // only process valid kmers
            if (start_pos + tid + i + k - 1 < sequence_len) 
            {
                // kmer index
                unsigned long long kmer_idx = kmer_to_index(kmer, k);
                
                // mapping to shared memory location
                unsigned int sharedIdx = kmer_idx % SHARED_MEM_SIZE;
                
                // updating count in shared memory
                atomicAdd(&shmem[sharedIdx], 1);
            }
        }
    }
    
    __syncthreads();
    
    // updating global memory with shared memory counts
    for (int i = tid; i < SHARED_MEM_SIZE; i += block_size) 
    {
        if (shmem[i] > 0)
            atomicAdd(&global_kmer_counts[i], shmem[i]);
    }
}

int main(int argc, char* argv[]) 
{
    if (argc < 3) 
    {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <k>" << std::endl;
        return 1;
    }
    
    std::string fasta_file = argv[1];
    int k = std::stoi(argv[2]);
    
    if (k > MAX_K) {
        std::cerr << "k value too large. Maximum supported is " << MAX_K << std::endl;
        return 1;
    }
    
    // reading sequence
    std::cout << "Reading sequence from " << fasta_file << "..." << std::endl;
    std::string sequence = read_fatsa_file(fasta_file);
    
    if (sequence.empty()) 
    {
        std::cerr << "Couldn't read sequence or sequence is empty" << std::endl;
        return 1;
    }
    
    std::cout << "Sequence length: " << sequence.length() << " bases" << std::endl;
    
    // allocating memory for sequence on device
    char* sequence_d;
    size_t sequence_size = sequence.length() * sizeof(char);
    cudaMalloc(&sequence_d, sequence_size);
    cudaMemcpy(sequence_d, sequence.c_str(), sequence_size, cudaMemcpyHostToDevice);
    
    // allocating memory for kmer counts on device
    unsigned int* kmer_counts_d;
    size_t kmer_counts_size = MAX_KMERS(k) * sizeof(unsigned int);
    
    // check if we have enough memory for the counts array
    if (k > 13) 
    {  // 4^13 = 67,108,864 possible kmers
        std::cerr << "Warning: Large k value may exceed available GPU memory" << std::endl;
    }
    
    cudaMalloc(&kmer_counts_d, kmer_counts_size);
    cudaMemset(kmer_counts_d, 0, kmer_counts_size);
    
    // setting up cuda kernel params
    int block_size = 256;
    int num_blocks = min(1024, (int)((sequence.length() - k + 1 + block_size - 1) / block_size));
    
    // size of coalesced memory chunk
    int coalesce_size = 32 * 4; // 4 warps
    
    std::cout << "Launching optimized GPU kernel with " << num_blocks << " blocks of " 
              << block_size << " threads each" << std::endl;
    
    // measuring time for k-mer counting
    auto start = std::chrono::steady_clock::now();
    
    // launching optimized kernel with shared memory
    size_t shared_mem_size = SHARED_MEM_SIZE * sizeof(unsigned int);

    count_kmers_coalesced_kernel<<<num_blocks, block_size, shared_mem_size>>>(sequence_d, sequence.length(), k, kmer_counts_d, coalesce_size);
    
    cudaDeviceSynchronize();
    
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // checking for kernel launch errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) 
    {
        std::cerr << "CUDA kernel error: " << cudaGetErrorString(err) << std::endl;
        return 1;
    }
    
    // copying results back to host
    unsigned int* kmer_counts_h = new unsigned int[MAX_KMERS(k)];
    cudaMemcpy(kmer_counts_h, kmer_counts_d, kmer_counts_size, cudaMemcpyDeviceToHost);
    
    // counting unique kmers
    unsigned long long unique_kmers = 0;
    for (unsigned long long i = 0; i < MAX_KMERS(k); i++) 
    {
        if (kmer_counts_h[i] > 0)
            unique_kmers++;
    }
    
    std::cout << "Found " << unique_kmers << " unique " << k << "-mers" << std::endl;
    std::cout << "Optimized GPU time: " << elapsed.count() << " seconds" << std::endl;
    
    // printing top 10 most frequent kmers
    print_top_kmers(kmer_counts_h, k, 10);
    
    // freeing memory
    cudaFree(sequence_d);
    cudaFree(kmer_counts_d);
    delete[] kmer_counts_h;
    
    return 0;
} 