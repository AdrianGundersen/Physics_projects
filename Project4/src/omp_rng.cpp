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
    std::vector<std::uint32_t> seeds;
    seeds.reserve(walkers);
    for (int i = 0; i < walkers; i++) {
        std::seed_seq ss{base_seed, static_cast<std::uint32_t>(i)}; // unique seed for each walker
        std::uint32_t s = 0; // temp variable
        ss.generate(&s, &s + 1);   // generate a new seed
        seeds.push_back(s);
    }
    return seeds;
    } // End parallel region
}
  

