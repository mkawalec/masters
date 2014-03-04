cmake_minimum_required(VERSION 2.8)
project(turb)
include_directories(src include alglib)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
#set(CXX mpicxx)

set(alglib_src
    alglib/alglibinternal.cpp
    alglib/alglibmisc.cpp
    alglib/ap.cpp
    alglib/dataanalysis.cpp
    alglib/diffequations.cpp
    alglib/fasttransforms.cpp
    alglib/integration.cpp
    alglib/interpolation.cpp
    alglib/linalg.cpp
    alglib/optimization.cpp
    alglib/solvers.cpp
    alglib/specialfunctions.cpp
    alglib/statistics.cpp
    )

set(alglib_headers
    alglib/statistics.h
    alglib/stdafx.h
    alglib/specialfunctions.h
    alglib/solvers.h
    alglib/optimization.h
    alglib/linalg.h
    alglib/interpolation.h
    alglib/integration.h
    alglib/fasttransforms.h
    alglib/diffequations.h
    alglib/dataanalysis.h
    alglib/ap.h
    alglib/alglibmisc.h
    alglib/alglibinternal.h
    )

add_library(
    Alglib
    STATIC
    ${alglib_src}
    ${alglib_headers}
    )

set(HEADER_FILES 
    include/helpers.hpp 
    include/exceptions.hpp
    include/JacobianElement.hpp
    include/Jacobian.hpp
    include/Base.hpp
    include/WriteCheck.hpp
    include/Integrator.hpp 
    include/integrators/PaperIntegrator.hpp
    include/integrators/PolymerIntegrator.hpp
    include/Serializer.hpp
    include/serializers/NormSerializer.hpp
    include/serializers/ComponentSerializer.hpp
    include/serializers/GenericSerializer.hpp
    include/searchers/NoTauSimple.hpp
    include/searchers/SimpleSearcher.hpp
    include/searchers/CMASimple.hpp
    include/searchers/FireflySimple.hpp
    include/Searcher.hpp
    include/Computer.hpp
    include/computers/SimpleComputer.hpp
    include/computers/MultirunComputer.hpp
    include/computers/DecayMultirunComputer.hpp
    include/computers/DecayPathComputer.hpp
    )

set(CPP_FILES
    src/helpers.cpp
    src/integrator.cpp
    src/Integrator.cpp
    src/integrators/PaperIntegrator.cpp
    src/integrators/PolymerIntegrator.cpp
    src/Computer.cpp
    src/Searcher.cpp
    src/searchers/SimpleSearcher.cpp
    src/searchers/NoTauSimple.cpp
    src/searchers/CMASimple.cpp
    src/searchers/FireflySimple.cpp
    src/serializers/NormSerializer.cpp
    src/serializers/ComponentSerializer.cpp
    src/serializers/GenericSerializer.cpp
    src/computers/SimpleComputer.cpp
    src/computers/DecayMultirunComputer.cpp
    src/computers/DecayPathComputer.cpp
    ${HEADER_FILES}
    )

set(TEST_FILES ${CPP_FILES} src/test-suite.cpp include/test_integrators.hpp)

message(status "${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(integrator ${CPP_FILES})
#add_executable(test-suite ${TEST_FILES})

set(MARCH native)
if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${MARCH} -ffast-math -ftree-vectorize -O3 -Wall -Wextra -Wno-ignored-qualifiers -pedantic -W -std=c++0x -g")
else()
    message(status "gcc not found!")
    return()
endif()

unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIRS CACHE)
unset(FFTW_INCLUDES CACHE)
unset(FFTW_LIBRARIES CACHE)
#unset(ALGLIB_INCLUDES CACHE)
#unset(ALGLIB_LIBRARIES CACHE)
unset(MPI_LIB CACHE)

find_package(FFTW REQUIRED)
find_package(MPI REQUIRED)
find_package(Armadillo)

message(status ${ARMADILLO_INCLUDE_DIRS})
if (test strequal ${ARMADILLO_INCLUDE_DIRS} "")
    set(ARMADILLO_INCLUDE_DIRS
        ${CMAKE_CURRENT_BINARY_DIR}/arma_include
        )
    set(ARMADILLO_LIBRARIES
        ${CMAKE_CURRENT_BINARY_DIR}/libarmadillo.so
        )
endif()

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      OFF)
set(Boost_USE_STATIC_RUNTIME     ON)
set(BOOST_ALL_DYN_LINK           OFF)   # force dynamic linking for all libraries



find_package(Boost 1.52 COMPONENTS unit_test_framework program_options date_time REQUIRED)
#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

message(status "** FFTW Include: ${FFTW_INCLUDES}")
message(status "** FFTW Libraries: ${FFTW_LIBRARIES}")

message(status "** MPI Include: ${MPI_CXX_INCLUDE_PATH}")
message(status "** MPI Libraries: ${MPI_CXX_LIBRARIES}")

message(status "** Armadillo include: ${ARMADILLO_INCLUDE_DIRS}")
message(status "** Armadillo libraries: ${ARMADILLO_LIBRARIES}")


INCLUDE_DIRECTORIES(
                 ${INCLUDE_DIRECTORIES}
                 ${Boost_INCLUDE_DIR} 
                 ${FFTW_INCLUDES} 
                 ${MPI_CXX_INCLUDE_PATH}
                 ${ARMADILLO_INCLUDE_DIRS}
                 )
LINK_DIRECTORIES(
                 ${Boost_LIBRARY_DIRS} 
                 ${FFTW_LIBRARIES} 
                 ${MPI_CXX_LIBRARIES} 
                 ${CMAKE_CURRENT_BINARY_DIR}
                 ${ARMADILLO_LIBRARIES}
                 )

target_link_libraries(integrator -lm ${FFTW_LIBRARIES} 
                                     ${Boost_LIBRARIES} 
                                     Alglib
                                     ${MPI_CXX_LIBRARIES}
                                     ${ARMADILLO_LIBRARIES}
                     )
#target_link_libraries(test-suite -lm ${FFTW_LIBRARIES} ${Boost_LIBRARIES} ${ALGLIB_LIBRARIES})

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)
