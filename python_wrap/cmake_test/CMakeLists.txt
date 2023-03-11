cmake_minimum_required(VERSION 3.10)
#set(CMAKE_INSTALL_PREFIX /home/i/innes/)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules)

project(fortran_lib)

enable_language(Fortran)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)


add_subdirectory(src)

add_cython_target(_add_cython _add_cython.pyx)
add_library(_add_cython MODULE ${_add_cython})

message(ADD_CYTHNO = ${_add_cython})
python_extension_module(_add_cython)

target_link_libraries(_add_cython add_cython_fortran)

target_include_directories(_add_cython PUBLIC ${CMAKE_SOURCE_DIR}/src)
#add_executable(c.exe c.f90)
install(TARGETS _add_cython LIBRARY DESTINATION add_cython)

#target_link_libraries(c.exe ab)

