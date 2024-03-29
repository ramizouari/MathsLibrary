# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
cmake_minimum_required(VERSION 3.12)
include_directories(.)
set(GCC_COVERAGE_COMPILE_FLAGS "-pthread --ffast-math -O3")
SET(CMAKE_CXX_FLAGS -pthread)
set(CMAKE_CXX_STANDARD 20)
project (MathsLibrary)
file(GLOB math_files CONFIGURE_DEPENDS "*.h" "integer.cpp" "complex.cpp" "real_field.cpp"  "guassian_integer.cpp" "absalg/*" "primitive_types/*" "nt/*" "nt/special/*.cpp")
add_executable(MathsLibrary ${math_files} "main.cpp")
add_executable(FactorisationExample ${math_files} "examples/factorisation.cpp")
