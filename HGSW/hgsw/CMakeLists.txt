add_library(
        hgsw_interface
        STATIC

        hgsw.hpp
        hgsw_algos.hpp
        hgsw_funcs.hpp
        
        rng/aes256ctr.h
        rng/rng.h
        rng/sampler.h      
)

set_target_properties(hgsw_interface PROPERTIES LINKER_LANGUAGE CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_executable(HGSW_test test_hgsw.cpp)
target_link_libraries(HGSW_test hgsw_interface)

