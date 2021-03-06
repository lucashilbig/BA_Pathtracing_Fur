#
# Try to find GLEW library and include path.
# Once done this will define
#
# GLEW_FOUND
# GLEW_INCLUDE_PATH
# GLEW_LIBRARY
#
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    MESSAGE("64-Bit Compiler erkannt")
    set(OpenGL_ROOT_ENV $ENV{OpenGLx64_ROOT})
else (CMAKE_SIZEOF_VOID_P EQUAL 8)
    MESSAGE("32-Bit Compiler erkannt")
    set(OpenGL_ROOT_ENV $ENV{OpenGL_ROOT})
endif (CMAKE_SIZEOF_VOID_P EQUAL 8)

IF (APPLE)
    set(OpenGL_ROOT_ENV $ENV{CVK_DEPENDENCIES_OSX}/OpenGL)
ENDIF (APPLE)

IF (WIN32)
    FIND_PATH(GLEW_INCLUDE_PATH GL/glew.h
            ${OpenGL_ROOT_ENV}/include/
            )

    FIND_LIBRARY(GLEW_LIBRARY
            NAMES glew32s.lib
            PATHS ${OpenGL_ROOT_ENV}/lib
            )


ELSEIF (APPLE)
    FIND_PATH(GLEW_INCLUDE_PATH GL/glew.h
            PATHS ${OpenGL_ROOT_ENV}/include/
            )

    FIND_LIBRARY(GLEW_LIBRARY
            NAMES libGLEW.a
            PATHS ${OpenGL_ROOT_ENV}/lib
            )

ELSE ()
    FIND_PATH(GLEW_INCLUDE_PATH GL/glew.h)
    FIND_LIBRARY(GLEW_LIBRARY
            NAMES GLEW glew32 glew glew32s PATH_SUFFIXES lib64)
ENDIF ()


SET(GLEW_FOUND "NO")
IF (GLEW_INCLUDE_PATH AND GLEW_LIBRARY)
    SET(GLEW_LIBRARIES ${GLEW_LIBRARY})
    SET(GLEW_FOUND "YES")
    message("EXTERNAL LIBRARY 'GLEW' FOUND")
ELSE ()
    message("ERROR: EXTERNAL LIBRARY 'GLEW' NOT FOUND")
ENDIF (GLEW_INCLUDE_PATH AND GLEW_LIBRARY)
