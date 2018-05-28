MACRO(GENERATE_SUBDIRS result curdir bindir)
    FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
    SET(names "")
    FOREACH (child ${children})
        IF (IS_DIRECTORY ${curdir}/${child})
            IF (NOT ${child} MATCHES "\\..*")
                if (EXISTS ${curdir}/${child}/CMakeLists.txt)
                    string(REPLACE " " "_" child ${child})
                    SET(names ${names} ${child})
                    message("BUILD FOR '${child}' GENERATED")
                endif ()
            ENDIF ()
        ENDIF ()
    ENDFOREACH ()
    SET(${result} ${${result}} ${names})
    FOREACH (n ${names})
        add_subdirectory(${curdir}/${n} ${bindir}/${n})
    ENDFOREACH ()
ENDMACRO()

MACRO(VALIADATE_PATH result)
    set(pathlist "${ARGN}")
    FOREACH (currentpath ${pathlist})
        IF (EXISTS ${currentpath})
            SET(${result} ${currentpath} CACHE PATH "Project specific path. Set manually if it was not found.")
            message("${result} FOUND AT ${currentpath}")
        ENDIF ()
    ENDFOREACH ()
    IF (NOT ${result})
        SET(${result} "PATH-NOTFOUND" CACHE PATH "Project specific path. Set manually if it was not found.")
        message("ERROR: ${result} NOT FOUND")
    ENDIF ()
ENDMACRO()