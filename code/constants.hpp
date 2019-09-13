#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cstddef>

namespace cnst {

constexpr unsigned int blockSize = 512;

// for PageRank in general
constexpr std::size_t normCheckFrequency{1};

// for CSR:
constexpr int numberThreadsPerRow = 1;

}


#endif // CONSTANTS_HPP
