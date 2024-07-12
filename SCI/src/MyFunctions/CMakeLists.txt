# include(FetchContent)

# FetchContent_Declare(fmt
#   GIT_REPOSITORY https://github.com/fmtlib/fmt.git
#   GIT_TAG 10.0.0
# )
# FetchContent_MakeAvailable(fmt)

add_library(My-Functions my-functions.cpp my-functions-class.cpp)
target_link_libraries(My-Functions
    PUBLIC
    SCI-LinearOT
    SCI-BuildingBlocks
    # fmt::fmt
)
