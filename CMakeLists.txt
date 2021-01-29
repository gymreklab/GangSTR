cmake_minimum_required(VERSION 3.16)

include(ExternalProject)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#add_compile_options(-Werror)

enable_testing()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


# set the project name
project(GangSTR VERSION 2.5.0)
configure_file(GangSTRConfig.h.in GangSTRConfig.h)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

######################### htslib and dependencies ############################

ExternalProject_Add(pnlopt
        PREFIX ${CMAKE_BINARY_DIR}/thirdparty/nlopt
        GIT_REPOSITORY "https://github.com/stevengj/nlopt.git"
        GIT_TAG "v2.7.0"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/thirdparty/nlopt -DBUILD_SHARED_LIBS=OFF
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        LOG_DOWNLOAD 1
        )
ExternalProject_Add(zlib
        PREFIX ${CMAKE_BINARY_DIR}/thirdparty/zlib
        GIT_REPOSITORY "https://github.com/madler/zlib.git"
        GIT_TAG "v1.2.8"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/thirdparty/zlib/src/zlib/configure --prefix=${CMAKE_BINARY_DIR}/thirdparty/zlib --static
        INSTALL_DIR ${CMAKE_BINARY_DIR}/thirdparty/zlib
        LOG_DOWNLOAD 1
        LOG_INSTALL 1
        )
ExternalProject_Add(htslib
        PREFIX ${CMAKE_BINARY_DIR}/thirdparty/htslib
        GIT_REPOSITORY "https://github.com/samtools/htslib.git"
        GIT_TAG "1.11"
        UPDATE_COMMAND autoreconf ${CMAKE_BINARY_DIR}/thirdparty/htslib/src/htslib
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/thirdparty/htslib/src/htslib/configure #--disable-bz2 --disable-lzma
        BUILD_COMMAND make
        INSTALL_COMMAND make install prefix=${CMAKE_BINARY_DIR}/thirdparty/htslib
	LOG_DOWNLOAD 1
)
ExternalProject_Add(pgsl
        PREFIX ${CMAKE_BINARY_DIR}/thirdparty/gsl
        GIT_REPOSITORY "https://github.com/ampl/gsl.git"
        GIT_TAG "v2.5.0"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/thirdparty/gsl -DBUILD_SHARED_LIBS=OFF
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        LOG_DOWNLOAD 1
        )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR}/thirdparty/nlopt/include)
include_directories(${CMAKE_BINARY_DIR}/thirdparty/gsl/include)
#set(nlopt_static ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib/lnlopt.a)
set(htslib_static ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib/libhts.a)
include_directories(${CMAKE_BINARY_DIR}/thirdparty/htslib/include)
include_directories(${CMAKE_BINARY_DIR}/thirdparty/zlib/include)
set(zlib_static ${CMAKE_BINARY_DIR}/thirdparty/zlib/lib/libz.a)
add_dependencies(htslib zlib)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)


######################### gangstr code ############################

add_subdirectory(src)

# add the executable
add_executable(${PROJECT_NAME} src/main_gangstr.cpp)


target_compile_features(GangSTR PRIVATE cxx_range_for)
target_include_directories(GangSTR PUBLIC "${PROJECT_BINARY_DIR}")
target_link_directories(GangSTR PUBLIC ${CMAKE_BINARY_DIR}/src PUBLIC ${CMAKE_BINARY_DIR}/thirdparty/zlib/lib PUBLIC ${CMAKE_BINARY_DIR}/thirdparty/gsl/lib PUBLIC ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib64 ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib PUBLIC ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib)

target_link_libraries(GangSTR Threads::Threads)
target_link_libraries(GangSTR hts)
target_link_libraries(GangSTR gstr)
target_link_libraries(GangSTR z)
target_link_libraries(GangSTR m nlopt)
target_link_libraries(GangSTR gsl)

#target_link_libraries(GangSTR gstr pthread)
#target_link_libraries (GangSTR ${CMAKE_THREAD_LIBS_INIT} gstr )
install(TARGETS GangSTR DESTINATION bin)

add_dependencies(htslib zlib)
add_dependencies(gstr htslib pnlopt pgsl)

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "") #required

set(CPACK_PACKAGE_VERSION "${GangSTR_VERSION_MAJOR}.${GangSTR_VERSION_MINOR}.${GangSTR_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_MAJOR "${GangSTR_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${GangSTR_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${GangSTR_VERSION_PATCH}")

set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")

set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/build/;${CMAKE_SOURCE_DIR}/.git/")

include(CPack)
