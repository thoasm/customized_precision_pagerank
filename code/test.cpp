#include "macros.hpp"
#include "types.hpp"
#include "split_pointer.hpp"

#include <iostream>
#include <iomanip>
#include <ios>
#include <type_traits>
#include <vector>


template<int First>
void readList(macros::compile_int_list<First>) {
    std::cout << First;
}


template<int First, int... RestArgs>
void readList(macros::compile_int_list<First, RestArgs...>) {
    std::cout << First << ' ';
    readList(macros::compile_int_list<RestArgs...>{});
}


template<int First>
void testTypes(macros::compile_int_list<First>) {
    std::cout << First << " Bytes: " << sizeof(typename type_of_size<First>::type) << std::endl;
}


template<int First, int... RestArgs>
void testTypes(macros::compile_int_list<First, RestArgs...>) {
    std::cout << First << " Bytes: " << sizeof(typename type_of_size<First>::type) << std::endl;
    testTypes(macros::compile_int_list<RestArgs...>{});
}


template<int First, MySegmentsNumberType SplitInt, typename SplitType, typename IndexType>
void printSplitPtr(macros::compile_int_list<First>,
    SplitPointer<SplitInt, SplitType> &splitVec,
    IndexType idx)
{
    std::cout << (splitVec.template read<First>(idx));
}


template<int First, int... RestArgs, MySegmentsNumberType SplitInt, typename SplitType, typename IndexType>
void printSplitPtr(macros::compile_int_list<First, RestArgs...>,
    SplitPointer<SplitInt, SplitType> &splitVec,
    IndexType idx)
{
    std::cout << (splitVec.template read<First>(idx)) << ' ';
    printSplitPtr(macros::compile_int_list<RestArgs...>{}, splitVec, idx);
}


int main()
{
    using curRange = macros::compile_int_list<1,2,3,4,5,6,7,8,9,10>;
    readList(curRange{});
    std::cout << std::endl;
    std::cout << "Type for size 1 is std::int8_t: " << (std::is_same<type_of_size<1>::type, std::int8_t>::value ? "true" : "false") << std::endl;
    std::cout << "Size 1 is default type: " << (std::is_same<type_of_size<1>::type, detail::default_type_of_size<1>>::value ? "true" : "false") << std::endl;
    std::cout << "Testing int16: " << (std::is_same<type_of_size<2>::type, std::int16_t>::value ? "true" : "false") << std::endl;
    std::cout << "Testing int32: " << (std::is_same<type_of_size<4>::type, std::int32_t>::value ? "true" : "false") << std::endl;
    std::cout << "Testing int64: " << (std::is_same<type_of_size<8>::type, std::int64_t>::value ? "true" : "false") << std::endl;
    testTypes(curRange{});
    
    constexpr std::size_t numberElements = 129;
    using rangeSegments = macros::compile_int_list<0,1,2,3,4,5,6,7,8>;
    using CurSplit = SplitPointer<8, double>;
    //using CurSplit = SplitPointer<1, double>;
    
    std::vector<std::uint8_t> memory(CurSplit::getMinimumBytesFor(numberElements));
    
    CurSplit splitVector;
    std::cout << "\n\nminSize: " << splitVector.minimumSize << "; Size for " << numberElements << " elements: " << CurSplit::getMinimumBytesFor(numberElements) << "; sizeof(SegmentType) = " << sizeof(decltype(splitVector)::SegmentType) << std::endl;
    splitVector.segmentPtr = reinterpret_cast<CurSplit::SegmentType *>(memory.data());
    
    auto generateValue = [](std::size_t i) {return CurSplit::OriginalType{1.0} / static_cast<CurSplit::OriginalType>(i + 1); };
    
    for (std::size_t i = 0; i < numberElements; ++i) {
        splitVector.writeAll(i, generateValue(i));
    }
    
    std::cout << "Reading vector:\n" << std::setprecision(16) << std::scientific;
    for (std::size_t i = 0; i < numberElements; ++i) {
        std::cout << std::setw(2) << i << ": Expected: " << generateValue(i) << '\n';
        printSplitPtr(rangeSegments{}, splitVector, i);
        std::cout << std::endl << std::endl;
    }
    
    //using mergedRange = macros::merge_compile_int_list<curRange, 5,6,7,8>::type;
    //readList(mergedRange{});
    //std::cout << std::endl;
    return 0;
}
