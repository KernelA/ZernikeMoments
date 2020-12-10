find_package(Boost 1.72 REQUIRED COMPONENTS filesystem program_options log log_setup)
find_package(SQLite3 REQUIRED)


add_executable(zernike3d 
	${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/include/binvox_reader.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/include/stdafx.h 
	${CMAKE_CURRENT_SOURCE_DIR}/src/compute_descriptors.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/include/compute_descriptors.h 
	${CMAKE_CURRENT_SOURCE_DIR}/include/loggers.h 
	${CMAKE_CURRENT_SOURCE_DIR}/include/binvox_utils.hpp  
	${CMAKE_CURRENT_SOURCE_DIR}/include/compute_sha256.h 
	${CMAKE_CURRENT_SOURCE_DIR}/include/path_tree.hpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/compute_sha256.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/sqlite_row.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/db.h
)
target_compile_features(zernike3d PRIVATE cxx_std_14)
target_include_directories(zernike3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_precompile_headers(zernike3d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/stdafx.h)
target_link_libraries(zernike3d PRIVATE 3DZM PRIVATE SQLite::SQLite3 PRIVATE Boost::log_setup PRIVATE Boost::log PRIVATE Boost::boost PRIVATE Boost::filesystem PRIVATE Boost::program_options PRIVATE Boost::dynamic_linking PRIVATE picosha2 PRIVATE sqlmoderncpp)

add_custom_command(TARGET zernike3d POST_BUILD 
  COMMAND "${CMAKE_COMMAND}" -E copy 
  	${CMAKE_CURRENT_SOURCE_DIR}/logsettings.ini
     "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION>/logsettings.ini" 
  COMMENT "Copying logsettings")
