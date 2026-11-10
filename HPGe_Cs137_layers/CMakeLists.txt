cmake_minimum_required(VERSION 3.16...3.21)
project(HPGe_Cs137)

find_package(Geant4 11.3 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

add_executable(HPGe_Cs137
  main.cc
  DetectorConstruction.cc
  PrimaryGeneratorAction.cc
  EventAction.cc
  SteppingAction.cc
  ActionInitialization.cc
)

target_link_libraries(HPGe_Cs137 PRIVATE ${Geant4_LIBRARIES})
target_compile_features(HPGe_Cs137 PRIVATE cxx_std_17)

file(COPY macros DESTINATION ${CMAKE_BINARY_DIR})
