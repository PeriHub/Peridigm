# add_subdirectory (unit_test)

# Collection of libs that we need to link with
set(PD_DAMAGE_SOURCES
    damage_utilities.cxx    
)

# PdMaterialUtilities library
add_library(PdDamageUtilities ${PD_DAMAGE_SOURCES})

IF (INSTALL_PERIDIGM)
   install(TARGETS PdDamageUtilities EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()
