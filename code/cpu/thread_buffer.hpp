#ifndef THREAD_BUFFER_HPP
#define THREAD_BUFFER_HPP

#include <omp.h>
#include <cstdlib>
#include <cinttypes>

class ThreadBuffer {

public:
static ThreadBuffer &getInstance() {
    static ThreadBuffer instance;
    return instance;
}

void reserveBufferSizePerThread(std::size_t sizeInBytes) {
    if (bufferSize >= sizeInBytes) {
        return;
    }
    if (byteBuffer != nullptr) {
        free(byteBuffer);
    }
    #pragma omp parallel
    #pragma omp single
    {
        stride = sizeInBytes + alignment - (sizeInBytes % alignment);
        posix_memalign(reinterpret_cast<void **>(&byteBuffer), alignment, omp_get_num_threads() * stride);
        bufferSize = omp_get_num_threads() * stride;
    }
}

std::int8_t *getBuffer() {
    if (byteBuffer == nullptr)
        return byteBuffer;
    
    return byteBuffer + omp_get_thread_num() * stride;
}


~ThreadBuffer() {
    free(byteBuffer); byteBuffer = nullptr;
    bufferSize = 0;
}

private:
    static constexpr std::size_t alignment = 4096;
    
    ThreadBuffer() : byteBuffer(nullptr), bufferSize{0}, stride{0} {}
    
    std::int8_t *byteBuffer;
    std::size_t bufferSize;
    std::size_t stride;
};






#endif // THREAD_BUFFER_HPP
