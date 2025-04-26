#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <vector>
#include <algorithm>

// function to read the fatsa file and extract the dna sequence
std::string read_fatsa_file(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Error: Couldn't open file:" << filename << std::endl;
        return "";
    }

    std::string line;
    std::string sequence;

    // skipping the header line 
    std::getline(file, line);

    // reading the sequence lines
    while (std::getline(file, line)) 
        sequence += line;

    file.close();
    return sequence;
}

// function to count k-mers using CPU
std::unordered_map<std::string, int> count_kmers(const std::string& sequence, int k)
{
    std::unordered_map<std::string, int> kmer_counts;

    // for this project, int would be fine, but if to test with larger sequences, we would probably need size_t 
    // slide window of size k across the sequence
    // kmer = sequence: 0 to k 
    //      kmer count++
    for (size_t i = 0; i < sequence.length() - k; i++)
    {
        std::string kmer = sequence.substr(i, k);
        kmer_counts[kmer]++;
    }

    return kmer_counts;
}

// function to print top N most frequent kmers 
void print_top_kmers(const std::unordered_map<std::string, int>& kmer_counts, int n)
{
    // converting map to vectors of pairs for sorting 
    std::vector<std::pair<std::string, int>> kmer_vector(kmer_counts.begin(), kmer_counts.end());

    // sorting by frequency (descending) - compare kmer counts
    auto compare = [](const std::pair<std::string, int>& kmer1, const std::pair<std::string, int>& kmer2) {
        return kmer1.second > kmer2.second;
    };

    std::sort(kmer_vector.begin(), kmer_vector.end(), compare);

    std::cout << "Top " << n << " most frequent k-mers: " << std::endl;
    int count = 0;

    // printing top N kmers 
    for (const auto& pair : kmer_vector) 
    {
        if (count >= n)
            break;
        std::cout << pair.first << ": " << pair.second << std::endl;
        count++;
    }
}

int main(int argc, char* argv[]) 
{
    if (argc < 3) 
    {
        std::cerr << "Usage: " << argv[0] << " <fatsa_file> <k>" << std::endl;
        return 1;
    }

    std::string fatsa_file = argv[1];
    int k = std::stoi(argv[2]);

    // reading sequence
    std::cout << "Read sequence from " << fatsa_file << "..." << std::endl;
    std::string sequence = read_fatsa_file(fatsa_file);

    if (sequence.empty()) 
    {
        std::cerr << "Couldn't read sequence" << std::endl;
        return 1;
    }

    std::cout << "Sequence length: " << sequence.length() << " bases" << std::endl; 

    // measuring time for kmer counting - https://www.geeksforgeeks.org/chrono-in-c/
    auto start = std::chrono::steady_clock::now();

    // counting kmers
    std::unordered_map<std::string, int> kmer_counts = count_kmers(sequence, k);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Found " << kmer_counts.size() << " unique " << k << "-mers" << std::endl;
    std::cout << "CPU time: " << elapsed.count() << " seconds" << std::endl;

    // print the top most freq kmers
    print_top_kmers(kmer_counts, 10);

    return 0;
}

