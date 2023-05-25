# For compiling part of the code without SOCRATES 
# Requires a dummy fortran file to act as socrates code
CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
PROJECT(soc-rad-conv-fort Fortran)

# Set directory of source files
set(ROOTDIR ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran)

#Find HDF5
find_package(HDF5 COMPONENTS Fortran HL)
SET(CMAKE_MODULE_PATH "${ROOTDIR}/cmake/Modules/")

# Set fortran flags, change this if you want
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# Find netcdf
INCLUDE(${CMAKE_MODULE_PATH}/FindNetCDF.cmake)
set (NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

include_directories(${NETCDF_INCLUDES} ${HDF5_INCLUDE_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/lib)
add_definitions(${HDF5_DEFINITIONS})

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

# Define executable name
set(EXE main.exe)  
set(SRC ${ROOTDIR}/src)
set(LIB ${CMAKE_SOURCE_DIR}/lib)
set(OBJ ${CMAKE_SOURCE_DIR}/obj)
set(BIN ${CMAKE_SOURCE_DIR}/bin)

# Put module files in the obj directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${OBJ})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})

# Dummy socrates
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB})
add_library(socrates ${CMAKE_SOURCE_DIR}/src/dummy_soc.f90)
target_link_libraries(socrates ${NETCDF_LIBRARIES_F90} ${HDF5_LIBRARIES})

# Find all source files
file(GLOB_RECURSE SRC_FILES
  "${SRC}/*.f90"
  "${SRC}/*.F90"
 )

 # Remove socrates to avoid double-counting
set (EXCLUDE_DIR "${ROOTDIR}/src/socrates/")
foreach (TMP_PATH ${SRC_FILES})
  string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM SRC_FILES ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)

# Remove archive directory from compilation
set (EXCLUDE_DIR "${ROOTDIR}/src/archive/")
foreach (TMP_PATH ${SRC_FILES})
  string (FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND)
    if (NOT ${EXCLUDE_DIR_FOUND} EQUAL -1)
        list (REMOVE_ITEM SRC_FILES ${TMP_PATH})
    endif ()
endforeach(TMP_PATH)

# Add executable
add_executable(${EXE} ${SRC_FILES})
target_compile_definitions(${EXE} PUBLIC SOC)

# Make sure socrates compiled before executable (is this necessary?)
add_dependencies(${EXE} socrates)

# Link SOCRATES library
target_link_libraries(${EXE}
  socrates
  ${NETCDF_LIBRARIES_F90}
  ${HDF5_LIBRARIES}
  )


