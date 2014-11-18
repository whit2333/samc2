#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(SAMC CXX Fortran)

# The version number.
set (SAMC_MAJOR_VERSION 0)
set (SAMC_MINOR_VERSION 1)
set (SAMC_PATCH_VERSION 0)

set (SAMC_CXX_CFLAGS)
set (SAMC_CXX_LDFLAGS)
set (SAMC_CXX_LIBS)
set (SAMC_CXX_INC_DIR)
set (SAMC_FC_CFLAGS)

set (SAMC_GRID_DATA_DIR)
set (SAMC_DATA_DIR "${PROJECT_SOURCE_DIR}/main")

#-------------------------------------------------------
set(CMAKE_MODULE_PATH
   ${PROJECT_SOURCE_DIR}/cmake/Modules 
   ${PROJECT_SOURCE_DIR}/cmake
   ${CMAKE_MODULE_PATH}
   )


#-------------------------------------------------------
# checkout out of source directory build
macro(SAMC_CHECK_OUT_OF_SOURCE_BUILD)
   get_filename_component(bindir_parent ${CMAKE_BINARY_DIR} PATH)
   if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
      file(REMOVE_RECURSE ${CMAKE_SOURCE_DIR}/Testing)
      file(REMOVE ${CMAKE_SOURCE_DIR}/DartConfiguration.tcl)
      message(FATAL_ERROR "SAMC should be built as an out of source build, to keep the source directory clean. Please create a extra build directory and run the command 'cmake <path_to_source_dir>' in this newly created directory. You have also to delete the directory CMakeFiles and the file CMakeCache.txt in the source directory. Otherwise cmake will complain even if you run it from an out-of-source directory.")
   elseif(IS_SYMLINK ${CMAKE_BINARY_DIR} AND CMAKE_SOURCE_DIR STREQUAL bindir_parent)
      message(FATAL_ERROR "SAMC cannot be built from a sub-directory of the source tree that is a symlink. This is a current limitation of CMake. Please create a real build directory and run the command 'cmake <path_to_source_dir>' in this newly created directory.")
   endif()
endmacro()
#Doesn't work with older cmake versions:
#SAMC_CHECK_OUT_OF_SOURCE_BUILD()

#-------------------------------------------------------
# Check the compiler and set the compile and link flags
#set(CMAKE_BUILD_TYPE Debug)

#-------------------------------------------------------
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

#-------------------------------------------------------
# not sure this is the best way to do this....
set(SAMC_FC_CFLAGS "${SAMC_FC_CFLAGS} ${CMAKE_Fortran_FLAGS}")
#set(CMAKE_Fortran_FLAGS " -D\"SAMC_PDF_GRID_DIR='${SAMC_PDF_GRID_DIR}'\" -g -O2 -Wall -fPIC -Iinclude -ffixed-line-length-none -x f77-cpp-input ")
set(CMAKE_Fortran_FLAGS " -D\"SAMC_PDF_GRID_DIR='${SAMC_PDF_GRID_DIR}'\" -g -O2 -Wall  -ffixed-line-length-none -x f77-cpp-input ")

#set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} ")

#----------------------------------------------------------------------------
# Find ROOT (required package)
find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
   message(STATUS "ROOT package not found. ")  
   return()
endif()

#----------------------------------------------------------------------------
# Find yaml-cpp (required package)
set(YAMLCPP_LIBRARIES "-lyaml-cpp")

#-------------------------------------------------------
set(LINK_DIRECTORIES
   ${ROOT_LIBRARY_DIR}
   )
link_directories( ${LINK_DIRECTORIES})

#-------------------------------------------------------
set(INCLUDE_DIRECTORIES
   ${ROOT_INCLUDE_DIR} )

include_directories(${PROJECT_SOURCE_DIR}/include 
   "${CMAKE_BINARY_DIR}"
   ${INCLUDE_DIRECTORIES}
   )

message(STATUS " ${INCLUDE_DIRECTORIES} ")  

#---------------------------------------------------------------------------
# ROOT compile flags
EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

EXECUTE_PROCESS(COMMAND insane-config --inc OUTPUT_VARIABLE INSANE_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INSANE_CXX_FLAGS}")

set(SAMC_CXX_CFLAGS "${SAMC_CXX_CFLAGS} ${CMAKE_CXX_FLAGS}")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.C)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the libraries 
#add_library(myObjects ${sources})


#---------------------------------------------------------------------------
# ROOT linking flags 
EXECUTE_PROCESS(COMMAND root-config --evelibs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ROOT_LD_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}  -lGeomPainter -lGeom -lEve -lTreeViewer -lSpectrum -lSpectrumPainter -lFoam ")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}   -Wl,-rpath-link,/home/whit/lib ")

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}  -lyaml-cpp ")
 

set(SAMC_CXX_LIBS "-L${CMAKE_INSTALL_PREFIX}/lib -lSAMC ${SAMC_CXX_LIBS} ")
set(SAMC_CXX_LDFLAGS "${SAMC_CXX_LDFLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
set(SAMC_CXX_INC_DIR "-I${CMAKE_INSTALL_PREFIX}/include/SAMC ")

# configure a header file to pass some of the CMake settings to the source code
configure_file (
   "${PROJECT_SOURCE_DIR}/SAMCConfig.h.in"
   "${PROJECT_BINARY_DIR}/SAMCConfig.h"
   )

#----------------------------------------------------------------------------
# Add the executable
add_executable(samc-config samc-config.cxx )

#target_link_libraries(SAMC ${myObjects} 
#   ${MYSQL_LIBRARIES} 
#   ${CLHEP_LIBRARIES} 
#   ${SQLITE3_LIBRARIES} 
#   ${ROOT_LIBRARIES} )

set(SAMC_LIB_NAMES
   General
   )

set(SAMC_SHARED_LIBS)
set(SAMC_ROOT_DICTS)

# ---------------------------------------
foreach(aSharedLib ${SAMC_LIB_NAMES})
   set(SAMC_SHARED_LIBS ${SAMC_SHARED_LIBS} "${CMAKE_BINARY_DIR}/lib/SAMC${aSharedLib}.so")
   set(SAMC_ROOT_DICTS ${SAMC_ROOT_DICTS} "SAMC${aSharedLib}Dict.h")
   set(SAMC_ROOT_DICT_SRCS ${SAMC_ROOT_DICT_SRCS} "SAMC${aSharedLib}Dict.cxx")
   if(NOT SAMC_USE_ROOT6)
       set(SAMC_ROOT_DICT_HEADERS ${SAMC_ROOT_DICT_HEADERS} "SAMC${aSharedLib}Dict.h")
   endif(NOT SAMC_USE_ROOT6)
endforeach(aSharedLib)

ADD_CUSTOM_TARGET(ROOTDICTS DEPENDS ${SAMC_ROOT_DICT_SRCS} ${SAMC_ROOT_DICT_HEADERS} )

foreach(aSharedLib ${SAMC_LIB_NAMES})
   include(${aSharedLib})
   add_dependencies(SAMC${aSharedLib} ROOTDICTS)
endforeach(aSharedLib)


FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${PROJECT_SOURCE_DIR}/include/*.h" "${PROJECT_BINARY_DIR}/*.h")

# -----------------------------------------------------------
# Set the library version in the main CMakeLists.txt
SET(SAMC_VERSION "${SAMC_MAJOR_VERSION}.${SAMC_MINOR_VERSION}.${SAMC_PATCH_VERSION}")
SET(SAMC_LIBRARY_PROPERTIES ${SAMC_LIBRARY_PROPERTIES}
    VERSION "${SAMC_VERSION}"
    SOVERSION "${SAMC_MAJOR_VERSION}"
    SUFFIX ".so"
)
set(SAMC_SRCS ${General_SRCS})

add_library(SAMC SHARED ${SAMC_SRCS})
target_link_libraries(SAMC ${ROOT_LIBRARIES})
add_dependencies(SAMC ROOTDICTS)
set_target_properties(SAMC PROPERTIES ${SAMC_LIBRARY_PROPERTIES})


# ------------------------------------------------------------------------------
# Build the binary and ins 
#add_executable(samc2 src/SAMC.C)
#add_dependencies(samc2 SAMC)
#target_link_libraries(samc2 ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES} SAMC)

add_executable(samc2 src/main.C)
add_dependencies(samc2 SAMC)
target_link_libraries(samc2 ${YAMLCPP_LIBRARIES} ${ROOT_LIBRARIES} SAMC)

# ------------------------------------------------------------------------------
# Install the shared libraries 

foreach(aLib ${SAMC_LIB_NAMES})
   install(TARGETS "SAMC${aLib}" DESTINATION   lib)
endforeach(aLib)

#install(TARGETS Cluster DESTINATION lib) #${CMAKE_BINARY_DIR}/lib)
#install(TARGETS Cut DESTINATION     lib) #${CMAKE_BINARY_DIR}/lib)
#install(TARGETS Event DESTINATION   lib) #${CMAKE_BINARY_DIR}/lib)
#install(TARGETS General DESTINATION lib) #${CMAKE_BINARY_DIR}/lib)
#install(TARGETS Geo DESTINATION     lib) #${CMAKE_BINARY_DIR}/lib)
#install(TARGETS Math DESTINATION    lib) #${CMAKE_BINARY_DIR}/lib)
#install(TARGETS Physics DESTINATION lib) #${CMAKE_BINARY_DIR}/lib)

# ------------------------------------------------------------------------------
# install the binaries and libs
install(TARGETS SAMC DESTINATION lib)
install(FILES ${files} DESTINATION include/SAMC)
install(TARGETS samc2 DESTINATION bin)
install(TARGETS samc-config DESTINATION bin)


# ------------------------------------------------------------------------------
# Build the binary 
MESSAGE("         CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
MESSAGE("                SAMC_SRCS: ${SAMC_SRCS}")
MESSAGE("            SAMC_CXX_LIBS: ${SAMC_CXX_LIBS}")
MESSAGE("         SAMC_CXX_LDFLAGS: ${SAMC_CXX_LDFLAGS}")
MESSAGE("         SAMC_CXX_INC_DIR: ${SAMC_CXX_INC_DIR}")
MESSAGE("             CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
MESSAGE("              FC_COMPILER: ${CMAKE_FC_COMPILER}")
MESSAGE("      CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
MESSAGE("                CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
MESSAGE("         CMAKE_LINK_FLAGS: ${CMAKE_LINK_FLAGS}")
MESSAGE("CMAKE_SHARED_LINKER_FLAGS: ${CMAKE_SHARED_LINKER_FLAGS}")
MESSAGE(" ROOT_LIBRARIES: ${ROOT_LIBRARIES}")


