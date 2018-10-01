cmake_minimum_required(VERSION 3.9)

enable_testing()

# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})


project(test_bspline_surface)
set(
    TESTS 
        "test_bspline_surface.cpp"
)
add_executable(${PROJECT_NAME} ${TESTS})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME} Catch BSplineSurface)