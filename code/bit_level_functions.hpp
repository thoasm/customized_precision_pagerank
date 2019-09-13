#ifndef BIT_LEVEL_FUNCTIONS_HPP
#define BIT_LEVEL_FUNCTIONS_HPP

#include "macros.hpp"

#include<cinttypes>

namespace bitLevelFunctions {


namespace detail {

template<typename T, int lowBitsToSet, int currentShift = lowBitsToSet - 1>
struct MaskCreationStruct {
    static_assert(currentShift < lowBitsToSet);
    static_assert(currentShift >= 0);
    static_assert(lowBitsToSet > 0);
    static constexpr T value = (T{1} << currentShift) | MaskCreationStruct<T, lowBitsToSet, currentShift-1>::value;
};


template<typename T>
struct MaskCreationStruct<T, 0, 0> {
    static constexpr T value = 0;
};
template<typename T>
struct MaskCreationStruct<T, 0, -1> {
    static constexpr T value = 0;
};

template<typename T, int lowBitsToSet>
struct MaskCreationStruct<T, lowBitsToSet, 0> {
    static constexpr T value = T{1};
};


template<typename T, std::int64_t number>
struct getTwoExponentStruct {
    static_assert(number > 0);
    static constexpr T value = T{1} + getTwoExponentStruct<T, (number >> 1)>::value;
};

template<typename T>
struct getTwoExponentStruct<T, 1> {
    static constexpr T value = T{0};
};

}   // namespace detail


template<typename T, T lowBitsToSet, T currentShift = 0>
MY_ATTRIBUTES constexpr T createMask() {
    return detail::MaskCreationStruct<T, lowBitsToSet>::value;
}


template<typename T, T number>
MY_ATTRIBUTES constexpr T getTwoExponent() {
    return detail::getTwoExponentStruct<T, number>::value;
}

}   // namespace bitLevelFunctions







#endif //BIT_LEVEL_FUNCTIONS_HPP
