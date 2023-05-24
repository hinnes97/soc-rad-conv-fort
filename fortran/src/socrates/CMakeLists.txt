file(GLOB_RECURSE SOC_SRC
  *.f90
  *.F90
  *.f
  *.finc)
list(FILTER SOC_SRC EXCLUDE REGEX "/src/aux/*" )
list(FILTER SOC_SRC EXCLUDE REGEX "/src/correlated_k/*")
list(FILTER SOC_SRC EXCLUDE REGEX "/src/general/*")
list(FILTER SOC_SRC EXCLUDE REGEX "src/scatter/*")

add_library(socrates SHARED ${SOC_SRC})

target_link_libraries(socrates ${NETCDF_LIBRARIES_F90}
  ${HDF5_LIBRARIES})
  

 
