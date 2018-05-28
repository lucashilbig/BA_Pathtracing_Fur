cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake CACHE PATH "Project specific path. Set manually if it was not found.")
set(DEPENDENCIES_OSX_PATH ${CMAKE_SOURCE_DIR}/../dependencies_OSX CACHE PATH "Project specific path. Set manually if it was not found.")
set(GUI_RESOURCES_PATH ${CMAKE_SOURCE_DIR}/../gui CACHE PATH "Project specific path. Set manually if it was not found.") 
set(EXECUTABLES_PATH ${CMAKE_SOURCE_DIR}/executables CACHE PATH "Project specific path. Set manually if it was not found.")
set(LIBRARIES_PATH ${CMAKE_SOURCE_DIR}/libraries CACHE PATH "Project specific path. Set manually if it was not found.")
set(SHADERS_PATH ${CMAKE_SOURCE_DIR}/shaders CACHE PATH "Project specific path. Set manually if it was not found.")
set(KERNELS_PATH ${CMAKE_SOURCE_DIR}/kernels CACHE PATH "Project specific path. Set manually if it was not found.")
set(UNITTESTS_PATH ${CMAKE_SOURCE_DIR}/unittests CACHE PATH "Project specific path. Set manually if it was not found.")
set(SCREENSHOTS_PATH ${CMAKE_SOURCE_DIR}/../screenshots CACHE PATH "Project specific path. Set manually if it was not found.")
set(LOGS_PATH ${CMAKE_SOURCE_DIR}/../logs CACHE PATH "Project specific path. Set manually if it was not found.")
set(SCENES_PATH ${CMAKE_SOURCE_DIR}/../scenes CACHE PATH "Project specific path. Set manually if it was not found.")
set(RESOURCES_PATH ${SCENES_PATH}/shared_resources CACHE PATH "Project specific path. Set manually if it was not found.")


include(${CMAKE_MODULE_PATH}/DefaultProject.cmake)
