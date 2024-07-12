cmake_minimum_required(VERSION 3.13)

project(PSACSEAL VERSION 4.1.0 LANGUAGES CXX)

if(NOT DEFINED BUILD_PSACSEAL)
    set(BUILD_PSACSEAL ON)

    # Import Microsoft SEAL
    # find_package(SEAL 4.1.1 EXACT REQUIRED)
    find_package(SEAL REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()

if(BUILD_PSACSEAL)



    add_library(ndsslib SHARED ${CMAKE_CURRENT_LIST_DIR}/Lu-ndss-HE.cpp ${CMAKE_CURRENT_LIST_DIR}/Lu-ndss-HE.h)
    if(TARGET SEAL::seal)
        target_link_libraries(ndsslib PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(ndsslib PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()

    add_executable(Lu-ndss-HE)
    target_sources(Lu-ndss-HE
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/Lu-ndss-HE.cpp
    )
    if(TARGET SEAL::seal)
    target_link_libraries(Lu-ndss-HE PRIVATE SEAL::seal fmt)
    elseif(TARGET SEAL::seal_shared)
    target_link_libraries(Lu-ndss-HE PRIVATE SEAL::seal_shared fmt)
    else()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()



    add_executable(U-test-HE)
    target_sources(U-test-HE
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/U-test-HE.cpp
    )
    if(TARGET SEAL::seal)
        target_link_libraries(U-test-HE PRIVATE SEAL::seal fmt)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(U-test-HE PRIVATE SEAL::seal_shared fmt)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()

endif()