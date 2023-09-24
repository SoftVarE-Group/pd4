cmake_minimum_required(VERSION 3.0)
project(gpmc CXX)
option(GPMC_STATIC_BUILD "build gpmc as a static execution file" off)

link_directories(/usr/local/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/flow-cutter-pace17/src /usr/local/include)
add_compile_options("$<$<CONFIG:Debug>:-O0>")
add_compile_options(-Wall -Wno-parentheses -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-result)
add_compile_options("$<$<CXX_COMPILER_ID:GNU>:-Wno-literal-suffix;-Wno-class-memaccess>")
add_compile_options("$<$<CXX_COMPILER_ID:Clang>:-Wno-reserved-user-defined-literal;-Wno-gnu-array-member-paren-init>")
add_compile_options("$<$<CXX_COMPILER_ID:AppleClang>:-Wno-redundant-move;-Wno-undefined-var-template;-Wno-reserved-user-defined-literal;-Wno-gnu-array-member-paren-init;-std=c++11>")


set(CMAKE_BUILD_TYPE Release)
add_subdirectory(flow-cutter-pace17)

# create an executable

#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)

add_library(gpmc SHARED 
	core/ComponentCache.cc core/ComponentManager.cc core/Counter.cc
	core/Config.cc core/Instance.cc ddnnf/DecisionTree.cc
	core/Solver.cc 
	utils/Options.cc utils/System.cc
	preprocessor/Preprocessor.cc preprocessor/TestSolver.cc 
	preprocessor/lib_sharpsat_td/subsumer.cpp
	preprocessor/TreeDecomposition.cc preprocessor/IFlowCutter.cc
    core/gpmc.cpp
    core/Solver.h
	)
target_compile_definitions(gpmc 
  PRIVATE $<$<CONFIG:Debug>:DEBUG>
           __STDC_FORMAT_MACROS __STDC_LIMIT_MACROS
)

target_include_directories(gpmc PUBLIC include)

target_include_directories(gpmc PRIVATE ${CMAKE_CURRENT_LIST_DIR} ../..)

target_link_libraries(gpmc PRIVATE flowcutter z mpfr gmp gmpxx $<$<BOOL:${GPMC_STATIC_BUILD}>:-static>)
target_link_libraries(gpmc PRIVATE arjun)
