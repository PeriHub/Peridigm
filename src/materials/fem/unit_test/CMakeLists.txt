add_executable(ut_FEM_routines ./ut_FEM_routines.cpp)
target_link_libraries(ut_FEM_routines
  ${FEMUtilitiesLib}
  ${Trilinos_LIBRARIES}
  ${PdMathUtilitiesLib}
)
add_test (ut_FEM_routines python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_FEM_routines)