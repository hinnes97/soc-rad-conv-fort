# Attempt to make CMake work with my project. Let's see how it goes

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(soc-rad-conv-fort Fortran)

#Find HDF5
find_package(HDF5 COMPONENTS Fortran HL)
message(HDF5=${HDF5_INCLUDE_DIRS})
message(HDF5=${HDF5_LIBRARIES})
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")


# Set fortran flags, change this if you want
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# Find netcdf
INCLUDE(${CMAKE_MODULE_PATH}/FindNetCDF.cmake)
set (NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

include_directories(${NETCDF_INCLUDES} ${HDF5_INCLUDE_DIRS})
add_definitions(${HDF5_DEFINITIONS})

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

# Define executable name
set(EXE main.exe)  
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(LIB ${CMAKE_SOURCE_DIR}/lib)
set(OBJ ${CMAKE_SOURCE_DIR}/obj)

# Put module files in the obj directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${OBJ})

# Put executable in main directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# Add separate socrates subdirectory
add_subdirectory(${SRC}/socrates ${OBJ})

# Find all source files
file(GLOB_RECURSE SRC_FILES
  "${SRC}/*.f90"
  "${SRC}/*.F90"
 )

 # Remove socrates to avoid double-counting
set (EXCLUDE_DIR "${CMAKE_SOURCE_DIR}/src/socrates/")
foreach (TMP_PATH ${SRC_FILES})
  string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM SRC_FILES ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)

foreach (TMP_PATH ${SRC_FILES})
  message(${TMP_PATH})
endforeach(TMP_PATH)
# Add executable
add_executable(${EXE} ${SRC_FILES})

# Make sure socrates compiled before executable (is this necessary?)
add_dependencies(${EXE} socrates)

# Link SOCRATES library
target_link_libraries(${EXE}
  socrates
  ${NETCDF_LIBRARIES_F90}
  ${HDF5_LIBRARIES}
  )

