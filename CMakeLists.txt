cmake_minimum_required(VERSION 3.17)
project(koval_stft)

set(ENABLE_PRECOMPILED_HEADERS OFF)

add_subdirectory(source)

add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)
