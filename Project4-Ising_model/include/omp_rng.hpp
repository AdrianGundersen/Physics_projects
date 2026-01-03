// include/omp_rng.hpp
// copied from lecture notes

#ifndef __omp_rng_hpp__  
#define __omp_rng_hpp__

#include "omp.h"
#include <vector>
#include <random>
#include <chrono>

// Put this in its own namespace, just to keep things tidy and neat
namespace omp_rng
{

  

  // Initialize the omp_rng_container according to the number of threads
  std::vector<uint32_t> initialize_omp_rng_container(uint32_t base_seed=-1, int walkers = 1);

 

}

#endif