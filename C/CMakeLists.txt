cmake_minimum_required(VERSION 2.6)
project(turb)
include_directories(src include)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(HEADER_FILES include/helpers.hpp include/Integrator.hpp)
set(TEST_HEADERS ${HEADER_FILES} include/test_integrators.hpp)

message(status "${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(integrator src/integrator.cpp src/Integrator.cpp src/helpers.cpp ${HEADER_FILES})
add_executable(test-suite src/test-suite.cpp src/test_integrators.cpp src/Integrator.cpp src/helpers.cpp ${TEST_HEADERS})

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)

if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static -O2 -Wall -Wextra -W -pedantic -std=c++0x")
else()
    message(status "gcc not found!")
    return()
endif()

#add_custom_target(debug ALL
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG -g")
#    COMMENT "Compiling in debug mode" VERBATIM
#    )

find_package(FFTW)

find_package(Boost 1.53 COMPONENTS unit_test_framework REQUIRED)
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

message(status "** FFTW Include: ${FFTW_INCLUDE_DIR}")
message(status "** FFTW Libraries: ${FFTW_LIBRARY_D}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} ${FFTW_INCLUDES})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

target_link_libraries(integrator -lm ${FFTW_LIBRARIES})
target_link_libraries(test-suite -lm ${FFTW_LIBRARIES} ${Boost_LIBRARIES})

set(Boost_USE_STATIC_LIBS        OFF)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
set(BOOST_ALL_DYN_LINK           ON)   # force dynamic linking for all libraries
