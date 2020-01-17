
find_package(Boost 1.71 REQUIRED COMPONENTS filesystem program_options)
             

add_executable(zernike3d main.cpp BinvoxReader.cpp)
target_compile_features(zernike3d PUBLIC cxx_std_14)
target_include_directories(zernike3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(zernike3d PUBLIC 3DZM PUBLIC Boost::boost PUBLIC Boost::filesystem PUBLIC Boost::program_options)
