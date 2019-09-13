#ifndef MACROS_HPP
#define MACROS_HPP

#include <utility>

#ifdef __CUDACC__
#define MY_ATTRIBUTES __host__ __device__
#define MY_INLINE __forceinline__
#define MY_RESTRICT __restrict__
#else
#define MY_ATTRIBUTES
#define MY_INLINE inline
#define MY_RESTRICT __restrict__
#endif  // __CUDACC__


namespace macros {

template <int... Values>
struct compile_int_list {
};


/*
//TODO
template <int Start, int End>
struct compile_int_range {
    static_assert(Start > End, "Start needs to be smaller than End!");
    //using list = ;
};


template<typename List, int... AdditionalValues>
struct merge_compile_int_list;

template<int... CurrentValues, int... AdditionalValues>
struct merge_compile_int_list<compile_int_list<CurrentValues...>, AdditionalValues...> {
    using type = compile_int_list<CurrentValues..., AdditionalValues...>;
};
*/


template <typename... Values>
struct compile_type_list {
};



#define GKO_ENABLE_IMPLEMENTATION_SELECTION(_name, _callable)                \
    template <typename Predicate, int... IntArgs, typename... TArgs,         \
              typename... InferredArgs>                                      \
    inline void _name(macros::compile_int_list<>, Predicate,                 \
                      macros::compile_int_list<IntArgs...>,                  \
                      macros::compile_type_list<TArgs...>,                   \
                      InferredArgs...) {throw "whoops"; };                   \
                                                                             \
    template <int K, int... Rest, typename Predicate, int... IntArgs,        \
              typename... TArgs, typename... InferredArgs>                   \
    inline void _name(macros::compile_int_list<K, Rest...>,                  \
                      Predicate is_eligible,                                 \
                      macros::compile_int_list<IntArgs...> int_args,         \
                      macros::compile_type_list<TArgs...> type_args,         \
                      InferredArgs... args)                                  \
    {                                                                        \
        if (is_eligible(K)) {                                                \
            _callable<IntArgs..., TArgs...>(                                 \
                macros::compile_int_list<K>(),                               \
                std::forward<InferredArgs>(args)...);                        \
        } else {                                                             \
            _name(macros::compile_int_list<Rest...>(), is_eligible,          \
                  int_args, type_args, std::forward<InferredArgs>(args)...); \
        }                                                                    \
    }


#define GKO_ENABLE_IMPLEMENTATION_ITERATION(_name, _callable)                \
    template <typename Predicate, int... IntArgs, typename... TArgs,         \
              typename... InferredArgs>                                      \
    inline void _name(macros::compile_int_list<>, Predicate,                 \
                      macros::compile_int_list<IntArgs...>,                  \
                      macros::compile_type_list<TArgs...>,                   \
                      InferredArgs...) {throw "whoops"; };                   \
                                                                             \
    template <int K, int... Rest, typename Predicate, int... IntArgs,        \
              typename... TArgs, typename... InferredArgs>                   \
    inline void _name(compile_int_list<K, Rest...>,                          \
                      Predicate is_satisfied,                                \
                      macros::compile_int_list<IntArgs...> int_args,         \
                      macros::compile_type_list<TArgs...> type_args,         \
                      InferredArgs... args)                                  \
    {                                                                        \
        if (is_satisfied(K)) {                                               \
            _callable<IntArgs..., TArgs...>(                                 \
                macros::compile_int_list<K>(),                               \
                std::forward<InferredArgs>(args)...);                        \
            _name(macros::compile_int_list<Rest...>(), is_satisfied,         \
                  int_args, type_args, std::forward<InferredArgs>(args)...); \
        }                                                                    \
    }

}   // namespace macros

#endif //MACROS_HPP
