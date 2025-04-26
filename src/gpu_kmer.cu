#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>

// max kmers supported
#define MAX_K 32

// calculate 4^k (max possible kmers for dna)
#define MAX_KMERS(k) (1ULL << (2 * (k)))

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

// naive gpu kernel for kmer counting 
__global__ void count_kmers_kernel(const char* sequence, int sequence_len, int k, unsigned int* kmer_counts, int chunk_size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    // Calculate chunk boundaries
    int chunk_start = blockIdx.y * chunk_size;
    int chunk_end = min(chunk_start + chunk_size, sequence_len - k + 1);
    
    // each thread processes multiple positions in its chunk
    for (int i = chunk_start + idx; i < chunk_end; i += stride) 
    {
        // calculating kmer index
        unsigned long long kmer_idx = kmer_to_index(&sequence[i], k);
        
        // update count using atomic operation
        if (kmer_idx != INVALID_KMER) {
            atomicAdd(&kmer_counts[kmer_idx], 1);
        }
    }
}

// function to read fatsa file and extract dna sequence 
std::string read_fasta_file(const std::string& filename) 
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Couldn't open file " << filename << std::endl;
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
    // vector of pairs for sorting 
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
    std::cout << "Top " << n << " most frequent k-mers: " << std::endl;

    char kmer_str[MAX_K + 1];
    int count = 0;

    for (const auto& pair : kmer_vector) 
    {
        if (count >= n)
            break;
        index_to_kmer(pair.first, k, kmer_str);
        std::cout << kmer_str << ": " << pair.second << std::endl;
        count++;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) 
    {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <k>" << std::endl;
        return 1;
    }
    
    std::string fasta_file = argv[1];
    int k = std::stoi(argv[2]);
    
    if (k > MAX_K) 
    {
        std::cerr << "k value too large. Maximum supported is " << MAX_K << std::endl;
        return 1;
    }
    
    // read the sequence 
    std::cout << "Reading sequence from " << fasta_file << "..." << std::endl;
    std::string sequence = read_fasta_file(fasta_file);
    
    if (sequence.empty()) 
    {
        std::cerr << "Could not read sequence or sequence is empty" << std::endl;
        return 1;
    }
    
    std::cout << "Sequence length: " << sequence.length() << " bases" << std::endl;
    
    // allocating memory for sequences on device
    char* sequence_d;
    size_t sequence_size = sequence.length() * sizeof(char);
    
    cudaMalloc(&sequence_d, sequence_size);
    cudaMemcpy(sequence_d, sequence.c_str(), sequence_size, cudaMemcpyHostToDevice);
    
    // allocating memory for kmer counts on device
    unsigned int* kmer_counts_d;
    size_t kmer_counts_size = MAX_KMERS(k) * sizeof(unsigned int);
    
    // check if we have enough memory for the counts array 
    if (k > 13) 
    {
        std::cerr << "Warning: Large k value may exceed available GPU memory" << std::endl;
        std::cerr << "Required memory: " << kmer_counts_size / (1024 * 1024) << " MB" << std::endl;
    }
    
    cudaError_t err = cudaMalloc(&kmer_counts_d, kmer_counts_size);
    if (err != cudaSuccess) 
    {
        std::cerr << "CUDA memory allocation error: " << cudaGetErrorString(err) << std::endl;
        std::cerr << "Failed to allocate " << kmer_counts_size / (1024 * 1024) << " MB for k-mer counts" << std::endl;
        return 1;
    }
    
    err = cudaMemset(kmer_counts_d, 0, kmer_counts_size);
    if (err != cudaSuccess) 
    {
        std::cerr << "CUDA memset error: " << cudaGetErrorString(err) << std::endl;
        cudaFree(kmer_counts_d);
        return 1;
    }
    
    // setting up cuda kernel params
    int block_size = 256;
    int num_blocks = 1024;  // fixed number of blocks per chunk
    int chunk_size = 1000000;  // processing 1 million bases at a time - this is because of the limited shared memory on the GPU
    
    std::cout << "Processing sequence in chunks of " << chunk_size << " bases" << std::endl;
    
    // measuring time
    auto start = std::chrono::steady_clock::now();
    
    // processing each chunk
    for (int chunk_start = 0; chunk_start < sequence.length() - k + 1; chunk_start += chunk_size) 
    {
        int num_chunks = (sequence.length() - k + 1 - chunk_start + chunk_size - 1) / chunk_size;
        dim3 grid(num_blocks, num_chunks);
        
        count_kmers_kernel<<<grid, block_size>>>(sequence_d, sequence.length(), k, kmer_counts_d, chunk_size);
        
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) 
        {
            std::cerr << "CUDA kernel error: " << cudaGetErrorString(err) << std::endl;
            cudaFree(sequence_d);
            cudaFree(kmer_counts_d);
            return 1;
        }
    }
    
    cudaDeviceSynchronize();
    
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // copying resutls back to host
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
    std::cout << "GPU time: " << elapsed.count() << " seconds" << std::endl;
    
    // print top 10 most frequent kmers
    print_top_kmers(kmer_counts_h, k, 10);
    
    // free memory
    cudaFree(sequence_d);
    cudaFree(kmer_counts_d);
    delete[] kmer_counts_h;
    
    return 0;
}
