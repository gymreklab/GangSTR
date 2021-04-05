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

find_package(PkgConfig REQUIRED)
pkg_check_modules(hts REQUIRED htslib)
pkg_check_modules(nlopt REQUIRED nlopt)
pkg_check_modules(gsl REQUIRED gsl)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)


######################### gangstr code ############################

add_subdirectory(src)

# add the executable
add_executable(${PROJECT_NAME} src/main_gangstr.cpp)


target_compile_features(GangSTR PRIVATE cxx_range_for)
target_include_directories(GangSTR PUBLIC "${PROJECT_BINARY_DIR}")
target_link_directories(GangSTR PUBLIC ${CMAKE_BINARY_DIR}/src)

target_link_libraries(GangSTR Threads::Threads)
target_link_libraries(GangSTR hts)
target_link_libraries(GangSTR gstr)
target_link_libraries(GangSTR m nlopt)
target_link_libraries(GangSTR gsl)

#target_link_libraries(GangSTR gstr pthread)
#target_link_libraries (GangSTR ${CMAKE_THREAD_LIBS_INIT} gstr )
install(TARGETS GangSTR DESTINATION bin)

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "") #required

set(CPACK_PACKAGE_VERSION "${GangSTR_VERSION_MAJOR}.${GangSTR_VERSION_MINOR}.${GangSTR_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_MAJOR "${GangSTR_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${GangSTR_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${GangSTR_VERSION_PATCH}")

set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")

set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/build/;${CMAKE_SOURCE_DIR}/.git/")

include(CPack)
