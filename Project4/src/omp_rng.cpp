// src/omp_rng.cpp
// copied from lecture notes
#include <random>
#include <chrono>
#include "omp_rng.hpp"
#include <vector>

using namespace std;

namespace omp_rng
{

  // Initialize the omp_rng_container according to the number of threads
  std::vector<uint32_t> initialize_omp_rng_container(uint32_t base_seed, int walkers) {
    std::seed_seq seq{base_seed, base_seed*2}; // need a range of numbers
    std::vector<std::uint32_t> seeds(walkers);

    seq.generate(seeds.begin(), seeds.end());
    
    return seeds;
  }

    

    } // End parallel region

  

