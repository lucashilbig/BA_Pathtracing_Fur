cmake_minimum_required(VERSION 2.8)
IF (${CMAKE_MAJOR_VERSION} GREATER "2")
    cmake_policy(SET CMP0022 OLD)
    cmake_policy(SET CMP0038 OLD)
ENDIF ()

get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})

include_directories(
        ${OpenGL3_INCLUDE_PATH}
        ${GTEST_INCLUDE_DIRS}
        ${GLEW_INCLUDE_PATH}
        ${GLFW3_INCLUDE_PATH}
        ${ASSIMP_INCLUDE_PATH}
        ${GLM_INCLUDE_PATH}
        ${EXTERNAL_LIBRARY_PATHS}
        ${LIBRARIES_PATH}
)

file(GLOB_RECURSE SOURCES *.cpp)
file(GLOB_RECURSE HEADER *.h)

add_definitions(-DSHADERS_PATH="${SHADERS_PATH}")
add_definitions(-DRESOURCES_PATH="${RESOURCES_PATH}")
add_definitions(-DGUI_RESOURCES_PATH="${GUI_RESOURCES_PATH}")
add_definitions(-DSCREENSHOTS_PATH="${SCREENSHOTS_PATH}")
add_definitions(-DLOGS_PATH="${LOGS_PATH}")
add_definitions(-DGLFW_INCLUDE_GLCOREARB)
add_definitions(-DGLEW_STATIC)
# glew.h sets __gl_h_ which makes gl3.h think /gl.h (OpenGL 2) is included. Calm that warning:
add_definitions(-DGL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED)
# ... and really don't include GLU and GL(2)
add_definitions(-DGLEW_NO_GLU)
add_definitions(-DGLM_FORCE_RADIANS)
# some functions (like printf) result in warnings in VS. Calm these warnings:
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

add_executable(${ProjectId} ${SOURCES} ${HEADER})

set_target_properties(${ProjectId} PROPERTIES FOLDER UnitTests)

if(APPLE)
target_link_libraries(
${ProjectId}
${ALL_LIBRARIES}
${ALL_LIBRARIES}
${ALL_LIBRARIES}
${GLFW3_LIBRARIES}
${ASSIMP_LIBRARIES}
${GLEW_LIBRARIES}
${OpenGL3_LIBRARIES}
${LibSTD_LIBRARIES}
${GTEST_BOTH_LIBRARIES}
)
else()
target_link_libraries(
        ${ProjectId}
        ${ALL_LIBRARIES}
        ${ALL_LIBRARIES}
        ${ALL_LIBRARIES}
        ${GLFW3_LIBRARIES}
        ${ASSIMP_LIBRARIES}
        ${GLEW_LIBRARIES}
        ${OpenGL3_LIBRARIES}
        ${GTEST_BOTH_LIBRARIES}
)
endif()
