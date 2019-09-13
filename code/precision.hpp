#ifndef PRECISION_HPP
#define PRECISION_HPP

#include <limits>
#include <cmath>

template<typename FloatType>
FloatType currentPrecision(MySegmentsNumberType currentNumberSegments, MySegmentsNumberType totalNumberSegments) {
    throw "Unsupported Type";
}

template<>
double currentPrecision<double>(MySegmentsNumberType currentNumberSegments, MySegmentsNumberType totalNumberSegments) {
    constexpr int totalBits = 64;
    constexpr int signBits = 1;
    constexpr int exponentBits = 11;
    //constexpr int mantissaBits = 52;
    
    int currentBitsAvailable = totalBits * currentNumberSegments / totalNumberSegments;
    if (currentBitsAvailable <= signBits) {
        return std::numeric_limits<double>::max();
    }
    else if (currentBitsAvailable <= signBits + exponentBits) {
        int usedExponentBits = currentBitsAvailable - signBits;
        return static_cast<double>(1 << (exponentBits - usedExponentBits));
    }
    else {
        int usedMantissaBits = currentBitsAvailable - signBits - exponentBits;
        return exp2(static_cast<double>(-usedMantissaBits));
    }
}

template<>
float currentPrecision<float>(MySegmentsNumberType currentNumberSegments, MySegmentsNumberType totalNumberSegments) {
    constexpr int totalBits = 32;
    constexpr int signBits = 1;
    constexpr int exponentBits = 8;
    //constexpr int mantissaBits = 23;
    
    int currentBitsAvailable = totalBits * currentNumberSegments / totalNumberSegments;
    if (currentBitsAvailable <= signBits) {
        return std::numeric_limits<float>::max();
    }
    else if (currentBitsAvailable <= signBits + exponentBits) {
        int usedExponentBits = currentBitsAvailable - signBits;
        return static_cast<float>(1 << (exponentBits - usedExponentBits));
    }
    else {
        int usedMantissaBits = currentBitsAvailable - signBits - exponentBits;
        return exp2(static_cast<float>(-usedMantissaBits));
    }
}

#endif // PRECISION_HPP
