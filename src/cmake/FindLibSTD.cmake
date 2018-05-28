#
# Try to find GLEW library and include path.
# Once done this will define
#
# LibSTD_FOUND
# LibSTD_INCLUDE_PATH
# LibSTD_LIBRARY
#

IF (APPLE)
    set(LibSTD_ROOT_ENV $ENV{CVK_DEPENDENCIES_OSX}/libstd)

    FIND_PATH(LibSTD_INCLUDE_PATH vector
        PATHS ${LibSTD_ROOT_ENV}/include/c++/v1/
        )

    FIND_LIBRARY(LibSTD_LIBRARY
        NAMES libc++experimental.a libc++.a
        PATHS ${LibSTD_ROOT_ENV}/lib/
        )

message(${LibSTD_LIBRARY})


ELSE()
    message("You shall not use this FindLibStdMAC on anything other than MAC.")
ENDIF ()



SET(LibSTD_FOUND "NO")
IF (LibSTD_INCLUDE_PATH AND LibSTD_LIBRARY)
    SET(LibSTD_LIBRARIES ${LibSTD_LIBRARY})
    SET(LibSTD_FOUND "YES")
    message("EXTERNAL LIBRARY 'LibSTD' FOUND")
ELSE ()
    message("ERROR: EXTERNAL LIBRARY 'LibSTD' NOT FOUND")
ENDIF (LibSTD_INCLUDE_PATH AND LibSTD_LIBRARY)
