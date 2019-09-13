#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstdint>



namespace detail {

template<std::size_t SizeInBytes>
struct default_type_of_size {
    std::int8_t bytes[SizeInBytes];
};

}


template<std::size_t SizeInBytes>
struct type_of_size {
    using type = detail::default_type_of_size<SizeInBytes>;
};

template<>
struct type_of_size<1> {
    using type = std::int8_t;
};

template<>
struct type_of_size<2> {
    using type = std::int16_t;
};

template<>
struct type_of_size<4> {
    using type = std::int32_t;
};

template<>
struct type_of_size<8> {
    using type = std::int64_t;
};


#endif
