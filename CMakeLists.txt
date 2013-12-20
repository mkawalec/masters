cmake_minimum_required(VERSION 2.8)
project(turb)
include_directories(src include)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(HEADER_FILES 
    include/helpers.hpp 
    include/exceptions.hpp
    include/JacobianElement.hpp
    include/Jacobian.hpp
    include/Base.hpp
    include/WriteCheck.hpp
    include/Integrator.hpp 
    include/Serializer.hpp
    include/serializers/NormSerializer.hpp
    include/serializers/ComponentSerializer.hpp
    include/serializers/GenericSerializer.hpp
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
    src/Computer.cpp
    src/Searcher.cpp
    src/serializers/NormSerializer.cpp
    src/serializers/ComponentSerializer.cpp
    src/serializers/GenericSerializer.cpp
    src/computers/SimpleComputer.cpp
    src/computers/DecayMultirunComputer.cpp
    src/computers/DecayPathComputer.cpp
    ${HEADER_FILES}
    #    ${alglib}
    )

set(TEST_FILES ${CPP_FILES} src/test-suite.cpp include/test_integrators.hpp)

message(status "${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(integrator ${CPP_FILES})
#add_executable(test-suite ${TEST_FILES})

if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -pthread -lpthread -static -Wall -Wextra -Wno-ignored-qualifiers -pedantic -W -std=c++0x -march=native -static")
else()
    message(status "gcc not found!")
    return()
endif()

#add_custom_target(debug ALL
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG -g")
#    COMMENT "Compiling in debug mode" VERBATIM
#    )
unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIRS CACHE)
unset(FFTW_INCLUDES CACHE)
unset(FFTW_LIBRARIES CACHE)

find_package(FFTW)
find_package(ALGLIB)

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     ON)
set(BOOST_ALL_DYN_LINK           OFF)   # force dynamic linking for all libraries

find_package(Boost 1.53 COMPONENTS unit_test_framework program_options random REQUIRED)
#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 
message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
message(status "** Boost Libraries: ${Boost_LIBRARIES}")

message(status "** FFTW Include: ${FFTW_INCLUDES}")
message(status "** FFTW Libraries: ${FFTW_LIBRARIES}")


message(status "** FFTW Include: ${ALGLIB_INCLUDES}")
message(status "** FFTW Libraries: ${ALGLIB_LIBRARIES}")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} ${FFTW_INCLUDES} ${ALGLIB_INCLUDES})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ${FFTW_LIBRARIES} ${ALGLIB_INCLUDES})


target_link_libraries(integrator -lm -pthread ${FFTW_LIBRARIES} ${Boost_LIBRARIES} ${ALGLIB_LIBRARIES})
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
