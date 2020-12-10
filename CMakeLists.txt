cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(ZernikeMoments LANGUAGES CXX VERSION 0.0.6)

add_subdirectory(lib)

add_library(picosha2 INTERFACE)
target_sources(picosha2 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/picosha2/picosha2.h)
target_include_directories(picosha2 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/picosha2)

add_library(sqlmoderncpp INTERFACE)

file(GLOB_RECURSE sqlmodern_headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sqlmoderncpp/hdr/*.h)
target_sources(sqlmoderncpp INTERFACE ${sqlmodern_headers})
target_include_directories(sqlmoderncpp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/sqlmoderncpp/hdr)

add_subdirectory(main)
