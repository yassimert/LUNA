cmake_minimum_required(VERSION 3.5.1)
project (hgsw_interface)

include_directories(.)

find_package(Palisade)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PALISADE_CXX_FLAGS} -std=c++17 -O3 -march=native -mtune=native")

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/third-party)
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )

link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )

set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )
link_libraries( ${PALISADE_SHARED_LIBRARIES}  )

add_subdirectory(hgsw)







