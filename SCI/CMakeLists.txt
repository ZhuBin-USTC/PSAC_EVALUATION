cmake_minimum_required (VERSION 3.13)

project (SCI)
set(NAME "SCI")


option(NO_REVEAL_OUTPUT "Only output secret shares after 2PC" OFF)
message(STATUS "Option: NO_REVEAL_OUTPUT = ${NO_REVEAL_OUTPUT}")

option(BUILD_TESTS "Build tests" ON)
message(STATUS "Option: BUILD_TESTS = ${BUILD_TESTS}")


option(BUILD_PSAC "Build PSACEvaluation" ON)
message(STATUS "Option: BUILD_PSAC = ${BUILD_PSAC}")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NO_REVEAL_OUTPUT)
    add_compile_definitions(NO_REVEAL_OUTPUT=1)
endif()

add_subdirectory(src)

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()


if (BUILD_PSAC)
    add_subdirectory(PSACEvaluation)
endif()

set(INSTALL_DIR "${CMAKE_INSTALL_PREFIX}")
set(BUILD_DIR "${PROJECT_SOURCE_DIR}/build")
include(CMakePackageConfigHelpers)
configure_package_config_file(
	"${PROJECT_SOURCE_DIR}/cmake/SCIConfig.cmake.in"
	"${PROJECT_SOURCE_DIR}/cmake/SCIConfig.cmake"
	INSTALL_DESTINATION lib/cmake/SCI
    PATH_VARS INSTALL_DIR BUILD_DIR)
install(FILES
        "${PROJECT_SOURCE_DIR}/cmake/SCIConfig.cmake"
        DESTINATION lib/cmake/SCI)
