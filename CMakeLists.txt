cmake_minimum_required(VERSION 3.9.0)

project(Analysis_Simtel VERSION 1.0)
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/compiled")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_INSTALL_PREFIX}/bin")
set(LIBRARY_OUTPUT_PATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_EXPORT_COMPILE_COMMANDS 1 CACHE BOOL "for clang" FORCE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS " -g -O0 -DH_MAX_TEL=39 -DH_MAX_PIX=8000 ")

set(HESS "/data/home/zhipz/hessioxxx/lib/libhessio.so")
find_package(ROOT  CONFIG REQUIRED COMPONENTS Minuit)
include("${ROOT_USE_FILE}")
root_generate_dictionary(Class ${PROJECT_SOURCE_DIR}/include/MonteCarloRunHeader.h  ${PROJECT_SOURCE_DIR}/include/TImage_Parameter.h  ${PROJECT_SOURCE_DIR}/include/LACTEvent.h ${PROJECT_SOURCE_DIR}/include/LACTTeldata.h ${PROJECT_SOURCE_DIR}/include/HitPix.h
		LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)


add_library(class SHARED) 
target_sources(class PRIVATE ${PROJECT_SOURCE_DIR}/src/MonteCarloRunHeader.cpp ${PROJECT_SOURCE_DIR}/src/TImage_Parameter.cpp ${PROJECT_SOURCE_DIR}/src/rec_tools.c   ${PROJECT_SOURCE_DIR}/src/LACTTeldata.cpp Class.cxx ${PROJECT_SOURCE_DIR}/src/LACTEvent.cpp ${PROJECT_SOURCE_DIR}/src/HitPix.cpp)
target_include_directories(class PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(class PRIVATE  ${ROOT_LIBRARIES} )

add_executable(Convert_Root)
target_sources(Convert_Root PUBLIC ${PROJECT_SOURCE_DIR}/src/convert_root.cpp  ${PROJECT_SOURCE_DIR}/src/LACTree.cpp ${PROJECT_SOURCE_DIR}/src/Detect_config.cpp ${PROJECT_SOURCE_DIR}/src/TCalibData.cpp)
target_include_directories(Convert_Root PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(Convert_Root PRIVATE class ${HESS} ${ROOT_LIBRARIES})

add_executable(Analysis)
target_sources(Analysis PUBLIC ${PROJECT_SOURCE_DIR}/src/Analysis.cpp ${PROJECT_SOURCE_DIR}/src/Detect_config.cpp ${PROJECT_SOURCE_DIR}/src/LACTree.cpp ${PROJECT_SOURCE_DIR}/src/TRecData.cpp ${PROJECT_SOURCE_DIR}/src/rec_tools.c ${PROJECT_SOURCE_DIR}/src/moments.cpp ${PROJECT_SOURCE_DIR}/src/TMcData.cpp ${PROJECT_SOURCE_DIR}/src/TCuts.cpp)
target_include_directories(Analysis PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(Analysis PRIVATE class ${ROOT_LIBRARIES})

add_executable(Angular)
target_sources(Angular PUBLIC ${PROJECT_SOURCE_DIR}/src/Angular.cpp ${PROJECT_SOURCE_DIR}/src/TRecData.cpp ${PROJECT_SOURCE_DIR}/src/rec_tools.c ${PROJECT_SOURCE_DIR}/src/LACTree.cpp ${PROJECT_SOURCE_DIR}/src/TMcData.cpp)
target_include_directories(Angular PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(Angular PRIVATE  ${ROOT_LIBRARIES})
