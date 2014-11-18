#---------------------------------------------------------------------------
set(General_files 
   SAMCMaterial
   SAMCEvent
   SAMCConstants
   SAMCManager
   SAMCEventGenerator
   SAMCPropagator
   #withxycollimator
   )
set(General_SRCS)
set(General_HEADERS)
foreach(infileName ${General_files})
   SET(General_SRCS ${General_SRCS} ${PROJECT_SOURCE_DIR}/src/${infileName}.C)
   SET(General_HEADERS ${General_HEADERS} ${PROJECT_SOURCE_DIR}/include/${infileName}.h)
endforeach(infileName)

#---------------------------------------------------------------------------
set(GeneralFortran_files
   Left_funcs
   Left_r-function
   Right_funcs
   Right_r-function
   monte_trans_hrs
   )
set(GeneralFortran_SRCS)
set(GeneralFortran_OBJECTS)
foreach(infileName ${GeneralFortran_files})
   SET(GeneralFortran_SRCS ${GeneralFortran_SRCS} "${PROJECT_SOURCE_DIR}/src/${infileName}.f")
   SET(GeneralFortran_OBJECTS ${GeneralFortran_OBJECTS} "${infileName}.o")
endforeach(infileName)

# set everything needed for the root dictonary and create the
# dictionary
set(General_LINKDEF ${PROJECT_SOURCE_DIR}/include/SAMCGeneral_LinkDef.h )
set(General_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/SAMCGeneralDict.cxx) 
ROOT_GENERATE_DICTIONARY("${General_HEADERS}" "${General_LINKDEF}" "${General_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

# add the dictionary to the list of source files
SET(General_SRCS ${GeneralFortran_SRCS} ${General_SRCS} ${General_DICTIONARY}) 

# Set the library version in the main CMakeLists.txt
SET(General_MAJOR_VERSION 0)
SET(General_MINOR_VERSION 0)
SET(General_PATCH_VERSION 0)
SET(General_VERSION "${General_MAJOR_VERSION}.${General_MINOR_VERSION}.${General_PATCH_VERSION}")
SET(General_LIBRARY_PROPERTIES ${General_LIBRARY_PROPERTIES}
    VERSION "${General_VERSION}"
    SOVERSION "${General_MAJOR_VERSION}"
    SUFFIX ".so"
)

add_library(          SAMCGeneral SHARED ${General_SRCS})
target_link_libraries(SAMCGeneral ${ROOT_LIBRARIES})
set_target_properties(SAMCGeneral PROPERTIES ${General_LIBRARY_PROPERTIES})

#install(TARGETS General DESTINATION ${CMAKE_BINARY_DIR}/lib)

