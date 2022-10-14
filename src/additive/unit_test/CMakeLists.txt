add_executable(ut_additive_utilities ./ut_additive_utilities.cpp)
target_link_libraries(ut_additive_utilities
    ${AdditiveUtilitiesLib}
    ${Trilinos_LIBRARIES}
)
add_test(ut_additive_utilities python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_additive_utilities)