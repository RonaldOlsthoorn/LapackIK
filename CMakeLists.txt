cmake_minimum_required(VERSION 3.15)
project(LapackIK)

set(CMAKE_CXX_STANDARD 14)

add_executable(LapackIK main.cpp model/Link.cpp model/Link.h model/RotationJoint.cpp model/RotationJoint.h model/PrismaticLink.cpp model/PrismaticLink.h model/Schunk.cpp model/Schunk.h)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(MKL REQUIRED)

if (MKL_FOUND)

    include_directories(${MKL_INCLUDE_DIR})
    # and for each of your dependent executable/library targets:
    set(MKL_LIBRARIES ${MKL_LIBRARIES} m dl pthread)
    target_link_libraries(LapackIK ${MKL_LIBRARIES})
endif()