# Executable to generate samples with dummy data
add_executable(GenerateSample generate_sample.cxx)
target_link_libraries(GenerateSample ROOT::ROOTVecOps ROOT::ROOTDataFrame logging)
target_include_directories(GenerateSample PRIVATE ${CMAKE_SOURCE_DIR})
install(TARGETS GenerateSample DESTINATION bin)
