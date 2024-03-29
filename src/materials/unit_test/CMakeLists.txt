add_subdirectory(twoPoint_SLS_Relaxation/)
add_subdirectory(utPeridigm_ElasticPlastic/)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string(COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)

if(NOT ${buildInSource})
  set(UT_DIRS ${UT_DIRS} -d src/materials/unit_test/twoPoint_SLS_Relaxation)
  set(UT_DIRS ${UT_DIRS} -d src/materials/unit_test/utPeridigm_ElasticPlastic)
  execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS})
endif(NOT ${buildInSource})

add_executable(utPeridigm_ElasticMaterial ./utPeridigm_ElasticMaterial.cpp

  # is the ) correct here?
)
target_link_libraries(utPeridigm_ElasticMaterial
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  ${PdMathUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(utPeridigm_ElasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_ElasticMaterial)

add_executable(utPeridigm_Material ./utPeridigm_Material.cpp)
target_link_libraries(utPeridigm_Material
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  ${PdMathUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(utPeridigm_Material python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_Material)

add_executable(utPeridigm_MultiphysicsElasticMaterial ./utPeridigm_MultiphysicsElasticMaterial.cpp)
target_link_libraries(utPeridigm_MultiphysicsElasticMaterial
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  ${PdMathUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(utPeridigm_MultiphysicsElasticMaterial python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./utPeridigm_MultiphysicsElasticMaterial)

add_executable(ut_correspondence ./ut_correspondence.cpp)
target_link_libraries(ut_correspondence
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  ${PdMathUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(ut_correspondence python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_correspondence)

add_executable(ut_elastic_correspondence ./ut_elastic_correspondence.cpp)
target_link_libraries(ut_elastic_correspondence
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  ${PdMathUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(ut_elastic_correspondence python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_elastic_correspondence)

add_executable(ut_material_utilities ./ut_material_utilities.cpp)
target_link_libraries(ut_material_utilities
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  ${PdMathUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(ut_material_utilities python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_material_utilities)

add_executable(ut_user_material_interface ./ut_user_material_interface.cpp)
target_link_libraries(ut_user_material_interface
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  ${PdMathUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(ut_user_material_interface python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_user_material_interface)

add_executable(ut_umats ./ut_umats.cpp)
target_link_libraries(ut_umats
  ${Peridigm_LIBRARY}
  ${PdMaterialUtilitiesLib}
  ${PdMathUtilitiesLib}
  PdField
  QuickGrid
  ${REQUIRED_LIBS}
  ${Trilinos_LIBRARIES}
)
add_test(ut_umats python ${CMAKE_BINARY_DIR}/scripts/run_unit_test.py ./ut_umats)
