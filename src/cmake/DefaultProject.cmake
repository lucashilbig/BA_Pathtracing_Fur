string(REPLACE "/" ";" p2list "${CMAKE_SOURCE_DIR}")
string(REPLACE "\\" ";" p2list "${p2list}")
list(REVERSE p2list)
list(GET p2list 0 first)
list(GET p2list 1 ProjectId)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})

include(${CMAKE_MODULE_PATH}/doxygen.cmake)
include(${CMAKE_MODULE_PATH}/macros.cmake)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_CONFIGURATION_TYPES Debug;Release)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif ()
if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8 /W1 /wd4267 /wd4312")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:MSVCRT")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /NODEFAULTLIB:LIBCMT /IGNORE:4099")
endif ()

find_package(OpenGL3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(GLFW3 REQUIRED)
find_package(ASSIMP REQUIRED)
find_package(GLM REQUIRED)

IF (APPLE)
find_package(LibSTD REQUIRED)
endif()

find_package(GTest)

if ("${CMAKE_SYSTEM}" MATCHES "Linux")
    find_package(X11)
    set(ALL_LIBRARIES ${ALL_LIBRARIES} ${X11_LIBRARIES} Xrandr Xxf86vm Xi pthread)
endif ()

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
GENERATE_SUBDIRS(ALL_LIBRARIES ${LIBRARIES_PATH} ${PROJECT_BINARY_DIR}/libraries)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
GENERATE_SUBDIRS(ALL_EXECUTABLES ${EXECUTABLES_PATH} ${PROJECT_BINARY_DIR}/executables)

if (GTEST_FOUND)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
    GENERATE_SUBDIRS(ALL_UNITTESTS ${UNITTESTS_PATH} ${PROJECT_BINARY_DIR}/unittests)
endif ()

if (EXISTS ${SHADERS_PATH})
    add_subdirectory(${SHADERS_PATH})
endif ()

file(MAKE_DIRECTORY ${LOGS_PATH})
