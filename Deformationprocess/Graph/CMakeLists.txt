
cmake_minimum_required( VERSION 3.8 )
PROJECT( Graph VERSION 1.0 DESCRIPTION "The Graph module of the ZuKai package")
set( Target_Name ZuKaiGraph )

#include(GNUInstallDirs)
message( "GNUInstallDirs: " ${GNUInstallDirs} )

set( MAIN_PROJECT_ROOT "${CMAKE_SOURCE_DIR}/.." )
# set( CMAKE_OSX_SYSROOT "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk")
set(CMAKE_CXX_STANDARD 17)

add_library(${Target_Name} SHARED
        src/BaseUndirectedGraph.cpp
        src/Graph.cpp
        src/Metro.cpp
        src/Guide.cpp
        src/BaseVertexProperty.cpp
        src/BaseEdgeProperty.cpp
        src/BaseGraphProperty.cpp
        )

set_target_properties(${Target_Name} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1
        PUBLIC_HEADER "src/Graph.h;src/Metro.h;src/BaseVertexProperty.h;src/BaseEdgeProperty.h;src/BaseGraphProperty.h;src/BaseUndirectedGraph.h"
        )

#----------------------------------------------------------------------
# path setting
#----------------------------------------------------------------------

target_include_directories(${Target_Name}
        PUBLIC
        ${CMAKE_SOURCE_DIR}/Base/src/
        )

link_directories(${Target_Name}
        PUBLIC
        ${CMAKE_BINARY_DIR}/Base/
        )

# set( CUSTOM_LIBRARIES "-lZuKaiBase")
target_link_libraries(${Target_Name} ZuKaiBase)

#----------------------------------------------------------------------
# Boost
#----------------------------------------------------------------------
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS system)

if(Boost_FOUND)
    message(STATUS "-- Linking boost")
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(${Target_Name} ${Boost_LIBRARIES} )
endif()

#----------------------------------------------------------------------
# install setting
#----------------------------------------------------------------------
set( INSTALL_LIBDIR ${CMAKE_SOURCE_DIR}/lib )
set( INSTALL_INCLUDEDIR ${CMAKE_SOURCE_DIR}/lib )

message( "INSTALL_LIBDIR: " ${INSTALL_LIBDIR} )
message( "INSTALL_INCLUDEDIR: " ${INSTALL_INCLUDEDIR} )

#configure_file(${Target_Name}.pc.in ${Target_Name}.pc @ONLY)

target_include_directories(${Target_Name} PRIVATE .)
install(TARGETS ${Target_Name}
        LIBRARY DESTINATION ${INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDEDIR} 
        )

#install(FILES ${CMAKE_BINARY_DIR}/mylib.pc
#        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
