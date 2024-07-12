# include(FetchContent)

# FetchContent_Declare(fmt
#   GIT_REPOSITORY https://github.com/fmtlib/fmt.git
#   GIT_TAG 10.0.0
# )
# FetchContent_MakeAvailable(fmt)

macro (PSAC_add_test_OT _name)
	add_executable(PSAC_OT_${_name} "PSAC_OT_${_name}.cpp")
    target_link_libraries(PSAC_OT_${_name} 
        SCI-BuildingBlocks 
        SCI-OT 
        SCI-HE
        My-Functions 
        fmt
        spdlog
    ) 
endmacro()

macro (PSAC_add_test_GC _name)
	add_executable(PSAC_GC_${_name} "PSAC_GC_${_name}.cpp")
    target_link_libraries(PSAC_GC_${_name} 
        SCI-BuildingBlocks 
        SCI-OT
        My-Functions 
        fmt
    ) 
endmacro()

# -------------------Evaluations for PSAC----------------------
PSAC_add_test_OT(min)
PSAC_add_test_OT(quantile)
PSAC_add_test_OT(quantile_hideK)
PSAC_add_test_OT(outlier)
PSAC_add_test_OT(chi2)
PSAC_add_test_OT(U)

PSAC_add_test_OT(min_rmind)
PSAC_add_test_OT(quantile_rmind)

PSAC_add_test_GC(min)
PSAC_add_test_GC(quantile_hideK)
PSAC_add_test_GC(sort_rmind)
PSAC_add_test_GC(chi2)

# -----------------------------------------

# macro(my_add_test_func _name)
#     add_executable(my-${_name}-func "my_test_func_${_name}.cpp")
#     target_link_libraries(my-${_name}-func
#         SCI-BuildingBlocks
#         SCI-OT
#         SCI-HE
#         SCI-FloatingPoint
#         My-Functions
#         fmt
#     ) 
# endmacro(my_add_test_func _name)
