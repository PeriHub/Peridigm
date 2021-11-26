#
# Add subdirectories
#
add_subdirectory (unit_test/)

# Collection of libs that we need to link with
set(PD_MATH_SOURCES
    matrices.cxx
)

# PdMathUtilities library
add_library(PdMathUtilities ${PD_MATH_SOURCES})

IF (INSTALL_PERIDIGM)
   install(TARGETS PdMathUtilities EXPORT peridigm-export
   LIBRARY DESTINATION ${LIB_INSTALL_DIR}
   ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
ENDIF ()