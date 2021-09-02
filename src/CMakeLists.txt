
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary sdsl malloc_count)
target_include_directories(compress_dictionary PUBLIC "../include/common")

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)
FetchContent_GetProperties(ksw2)
FetchContent_GetProperties(klib)
FetchContent_GetProperties(bigbwt)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)


add_executable(ms matching_statistics.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(ms common sdsl divsufsort divsufsort64 malloc_count ri pthread)
target_include_directories(ms PUBLIC    "../include/ms" 
                                        "../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(ms PUBLIC "-std=c++17")

add_executable(mems mems.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(mems common sdsl divsufsort divsufsort64 malloc_count ri pthread)
target_include_directories(mems PUBLIC    "../include/ms" 
                                        "../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(mems PUBLIC "-std=c++17")

add_executable(rlebwt_ms_build rlebwt_ms_build.cpp)
target_link_libraries(rlebwt_ms_build common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(rlebwt_ms_build PUBLIC    "../include/ms" 
                                        "../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(rlebwt_ms_build PUBLIC "-std=c++17")

add_executable(extend_klib extend_klib.cpp ${klib_SOURCE_DIR}/ksw.c ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(extend_klib common malloc_count sdsl divsufsort divsufsort64 ri klib ssw pthread)
target_include_directories(extend_klib PUBLIC    "../include/ms" 
                                            "../include/common" 
                                            "../include/extender" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                            "${klib_SOURCE_DIR}"
                                            "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(extend_klib PUBLIC "-std=c++17")

add_executable(extend_ksw2 extend_ksw2.cpp ${bigbwt_SOURCE_DIR}/xerrors.c)
target_link_libraries(extend_ksw2 common sdsl malloc_count divsufsort divsufsort64 ri ksw2 pthread)
target_include_directories(extend_ksw2 PUBLIC    "../include/ms" 
                                            "../include/common"
                                            "${ksw2_SOURCE_DIR}"
                                            "../include/extender" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${bigbwt_SOURCE_DIR}"
                                        )
target_compile_options(extend_ksw2 PUBLIC "-std=c++17")

add_executable(build_seqidx build_seqidx.cpp)
target_link_libraries(build_seqidx common sdsl divsufsort divsufsort64 malloc_count klib z)
target_include_directories(build_seqidx PUBLIC    "../include/ms" 
                                        "../include/common" 
                                        )
target_compile_options(build_seqidx PUBLIC "-std=c++17")