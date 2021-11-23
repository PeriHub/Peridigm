add_executable(ut_matrices ./ut_matrices.cpp)
target_link_libraries(ut_matrices
  ${PdMathUtilitiesLib}
  ${Trilinos_LIBRARIES}
)
add_test (ut_matrices python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_matrices)