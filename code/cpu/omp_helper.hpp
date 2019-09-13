#ifndef OMP_HELPER_HPP
#define OMP_HELPER_HPP

#define OMP_ENABLE_


namespace ompHelper {

constexpr int chunkSize = 64;   //Must be a power of 2
constexpr int blas1ChunkSize = 8 * 1024 / 4;

}

#endif
