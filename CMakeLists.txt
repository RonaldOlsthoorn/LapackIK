cmake_minimum_required(VERSION 3.15)
project(LapackIK)

set(CMAKE_CXX_STANDARD 14)

add_executable(LapackIK main.cpp model/Link.cpp model/Link.h model/RotationJoint.cpp model/RotationJoint.h model/PrismaticLink.cpp model/PrismaticLink.h model/Schunk.cpp model/Schunk.h model/Trajectory.cpp model/Trajectory.h model/operations.h model/operations.cpp model/InverseOperations.h model/InverseOperations.cpp)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(MKL REQUIRED)

if (MKL_FOUND)

    include_directories(${MKL_INCLUDE_DIR})
    # and for each of your dependent executable/library targets:
    set(MKL_LIBRARIES ${MKL_LIBRARIES} m dl pthread)
    target_link_libraries(LapackIK ${MKL_LIBRARIES})
endif()

find_package(PythonLibs 2.7)
target_include_directories(LapackIK PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(LapackIK ${PYTHON_LIBRARIES})

target_include_directories(LapackIK PRIVATE ../matplotlib-cpp)