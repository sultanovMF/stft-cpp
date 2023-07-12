find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)


find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)

pybind11_add_module(koval_stft stft.cpp stft.hpp)

target_compile_definitions(${PROJECT_NAME} PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_include_directories(${PROJECT_NAME} PRIVATE ${PYBIND11_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)