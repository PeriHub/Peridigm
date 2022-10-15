# add_subdirectory (unit_test)
add_subdirectory(unit_test)

# Collection of libs that we need to link with
set(PD_ADDITIVE_SOURCES
    additive_utilities.cxx
)

# PdMaterialUtilities library
add_library(PdAdditiveUtilities ${PD_ADDITIVE_SOURCES})

IF(INSTALL_PERIDIGM)
    install(TARGETS PdAdditiveUtilities EXPORT peridigm-export
        LIBRARY DESTINATION ${LIB_INSTALL_DIR}
        ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF()
