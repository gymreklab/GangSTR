file(GLOB SOURCES *.cpp *.h *.c)
add_library(gstr ${SOURCES})
target_link_libraries(gstr ${htslib_static} ${zlib_static} ${nlopt_static})
target_include_directories(gstr PUBLIC "${PROJECT_BINARY_DIR}")
