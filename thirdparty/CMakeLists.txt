include(FetchContent)

## Add malloc_count
FetchContent_Declare(
  malloc_count
  GIT_REPOSITORY https://github.com/bingmann/malloc_count
  )
  
FetchContent_GetProperties(malloc_count)
if(NOT malloc_count_POPULATED)
  FetchContent_Populate(malloc_count)

  add_library(malloc_count OBJECT ${malloc_count_SOURCE_DIR}/malloc_count.c ${malloc_count_SOURCE_DIR}/malloc_count.h)
  target_link_libraries(malloc_count dl)
  target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")

  add_library(memprofile OBJECT ${malloc_count_SOURCE_DIR}/memprofile.h)
  target_include_directories(memprofile PUBLIC "${malloc_count_SOURCE_DIR}")
endif()

  ## Add klib
FetchContent_Declare(
  klib
  GIT_REPOSITORY https://github.com/attractivechaos/klib
)

FetchContent_GetProperties(klib)
if(NOT klib_POPULATED)
  FetchContent_Populate(klib)
  
  # add_subdirectory(${klib_SOURCE_DIR} ${klib_BINARY_DIR} EXCLUDE_FROM_ALL)
  add_library(klib INTERFACE)

  target_include_directories(klib INTERFACE ${klib_SOURCE_DIR})
endif()

## Add Big-BWT
FetchContent_Declare(
  bigbwt
  GIT_REPOSITORY https://github.com/alshai/Big-BWT.git
  )
  
FetchContent_GetProperties(bigbwt)
if(NOT bigbwt_POPULATED)
  FetchContent_Populate(bigbwt)
  add_subdirectory(${bigbwt_SOURCE_DIR} ${bigbwt_BINARY_DIR})

  endif()


  
## Add gsacak
FetchContent_Declare(
  gsacak
  GIT_REPOSITORY https://github.com/felipelouza/gsa-is.git
)

FetchContent_GetProperties(gsacak)
if(NOT gsacak_POPULATED)
    FetchContent_Populate(gsacak)
    add_library(gsacak OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak PUBLIC "${gsacak_SOURCE_DIR}")
    
    add_library(gsacak64 OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak64 PUBLIC "${gsacak_SOURCE_DIR}")
    target_compile_options(gsacak64 PUBLIC -DM64)
endif()


## Add sdsl
FetchContent_Declare(
  sdsl
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
  FetchContent_Populate(sdsl)

  set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")
  
  add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

## Add divsuffsort
FetchContent_Declare(
  divsufsort
  GIT_REPOSITORY https://github.com/simongog/libdivsufsort.git
  GIT_TAG        2.0.1
)

FetchContent_GetProperties(divsufsort)
if(NOT divsufsort_POPULATED)
  FetchContent_Populate(divsufsort)

  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
  set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
  set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")

  add_subdirectory(${divsufsort_SOURCE_DIR} ${divsufsort_BINARY_DIR} EXCLUDE_FROM_ALL)

  target_include_directories(divsufsort PUBLIC "${divsufsort_BINARY_DIR}/include")
  target_include_directories(divsufsort64 PUBLIC "${divsufsort_BINARY_DIR}/include")
endif()


## Add r-index
FetchContent_Declare(
  r-index
  GIT_REPOSITORY https://github.com/maxrossi91/r-index.git
)

FetchContent_GetProperties(r-index)
if(NOT r-index_POPULATED)
  FetchContent_Populate(r-index)

  add_subdirectory(${r-index_SOURCE_DIR} ${r-index_BINARY_DIR} EXCLUDE_FROM_ALL)
  add_library(ri INTERFACE)
  target_link_libraries(ri INTERFACE klib z)
  target_include_directories(ri INTERFACE ${r-index_SOURCE_DIR}/internal)
endif()

## Add pfp-thresholds
FetchContent_Declare(
  pfp_thresholds
  GIT_REPOSITORY https://github.com/maxrossi91/pfp-thresholds.git
  GIT_TAG develop
  )
  
FetchContent_GetProperties(pfp_thresholds)
if(NOT pfp_thresholds_POPULATED)
  FetchContent_Populate(pfp_thresholds)
  add_subdirectory(${pfp_thresholds_SOURCE_DIR} ${pfp_thresholds_BINARY_DIR})
  
endif()

## Add bigrepair
FetchContent_Declare(
  bigrepair
  GIT_REPOSITORY https://gitlab.com/maxrossi91/bigrepair.git
  # GIT_REPOSITORY https://gitlab.com/manzai/bigrepair.git
  )
  
FetchContent_GetProperties(bigrepair)
if(NOT bigrepair_POPULATED)
  set(DISABLE_PFP ON CACHE BOOL "Build bigrepair without the PFP")
  FetchContent_Populate(bigrepair)
  add_subdirectory(${bigrepair_SOURCE_DIR} ${bigrepair_BINARY_DIR})

  # execute_process(COMMAND make
  #       RESULT_VARIABLE result
  #       WORKING_DIRECTORY ${bigrepair_SOURCE_DIR} )
  # if(result)
  #     message(FATAL_ERROR "CMake step for bigrepair failed: ${result}")
  # endif()
  
endif()

## Add ShapedSlp
FetchContent_Declare(
  shaped_slp
  GIT_REPOSITORY https://github.com/koeppl/ShapedSlp.git
  GIT_TAG master
  )
  
FetchContent_GetProperties(shaped_slp)
if(NOT shaped_slp_POPULATED)
  FetchContent_Populate(shaped_slp)
  add_subdirectory(${shaped_slp_SOURCE_DIR} ${shaped_slp_BINARY_DIR})
  set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
  set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)
endif()

## Add SSW
FetchContent_Declare(
    ssw
    GIT_REPOSITORY https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
    GIT_TAG master
  )
  
  FetchContent_GetProperties(ssw)
  if(NOT ssw_POPULATED)
    FetchContent_Populate(ssw)

    add_library(ssw OBJECT ${ssw_SOURCE_DIR}/src/ssw_cpp.cpp ${ssw_SOURCE_DIR}/src/ssw.c)
    target_include_directories(ssw PUBLIC ${ssw_SOURCE_DIR}/src)
  endif()


## Add Ksw2
FetchContent_Declare(
    ksw2
    GIT_REPOSITORY https://github.com/lh3/ksw2
    GIT_TAG master
  )
  
  FetchContent_GetProperties(ksw2)
  if(NOT ksw2_POPULATED)
    FetchContent_Populate(ksw2)
    
    add_library(ksw2 OBJECT ${ksw2_SOURCE_DIR}/kalloc.c 
                            ${ksw2_SOURCE_DIR}/ksw2_gg.c 
                            ${ksw2_SOURCE_DIR}/ksw2_gg2.c 
                            ${ksw2_SOURCE_DIR}/ksw2_gg2_sse.c 
                            ${ksw2_SOURCE_DIR}/ksw2_extz.c 
                            ${ksw2_SOURCE_DIR}/ksw2_extz2_sse.c
                            ${ksw2_SOURCE_DIR}/ksw2_extd.c 
                            ${ksw2_SOURCE_DIR}/ksw2_extd2_sse.c 
                            ${ksw2_SOURCE_DIR}/ksw2_extf2_sse.c 
                            ${ksw2_SOURCE_DIR}/ksw2_exts2_sse.c)
    target_include_directories(ksw2 PUBLIC ${ksw2_SOURCE_DIR}/src)
  endif()
  