cmake_minimum_required(VERSION 3.10)
project(MatlabToCpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找Eigen库
find_package(Eigen3 REQUIRED)

# 设置OpenMP
find_package(OpenMP REQUIRED)

# 添加可执行文件
add_executable(main_program
    main.cpp
    math_utils.cpp
    coordinate_conversion.cpp
    grid_utils.cpp
    ray_tracing_optimized.cpp
    sparse_matrix.cpp
    main_functions_optimized.cpp
)

# 链接库
target_link_libraries(main_program OpenMP::OpenMP_CXX)
target_link_libraries(main_program Eigen3::Eigen)
target_link_libraries(main_program m)

# 设置编译选项（简化版本，确保兼容性）
target_compile_options(main_program PRIVATE 
    -O3 
    -march=native 
    -ffast-math
    -DNDEBUG
)