cmake_minimum_required( VERSION 3.8 )
PROJECT( MetroShapes )
set( Target_Name MetroShapes )

#----------------------------------------------------------------------
# cmake setting
#----------------------------------------------------------------------
message( "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR} )
message( "CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR} )
# set output folder
# set(CMAKE_INCLUDE_CURRENT_DIR ON)
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
#set( SEPARATE_BUILD ON)
set(CMAKE_CXX_STANDARD 17)

set( MAIN_PROJECT_ROOT "${CMAKE_SOURCE_DIR}/.." )

#----------------------------------------------------------------------
# compiler setting
#----------------------------------------------------------------------
set(CMAKE_C_COMPILER /usr/bin)

if( APPLE )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -framework Accelerate -DEIGEN_USE_BLAS" )
    #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline -g -Wall -DIL_STD -O3 -DNDEBUG -lboost_thread-mt -lcgal -lgmp" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline -g -Wall -Wno-deprecated -DIL_STD -O3 -DNDEBUG" )
    #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline -g -Wall -Wno-deprecated -DIL_STD -O3 -DNDEBUG -openmp -omp" )
    #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline -g -Wall -Wno-deprecated -DIL_STD" )
elseif( UNIX )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -lboost_graph" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline -g -Wall -DIL_STD -O3 -DNDEBUG" )
    #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline -g -Wall -Wno-deprecated -DIL_STD -O3 -DNDEBUG" )
    #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline -g -Wall -Wno-deprecated -DIL_STD -O3 -DNDEBUG -openmp -omp" )
    #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline -g -Wall -Wno-deprecated -DIL_STD" )
elseif (WIN32)
    set( CMAKE_CXX_STANDARD 11)
endif()

#----------------------------------------------------------------------
# sources
#----------------------------------------------------------------------
set(HEADERS
        core/MetroShapes.h
#        core/RNGraph.h
        optimization/Smooth.h
        optimization/Mixedlayout.h
        matching/Stationmatching.h
        )
set(SOURCES
        core/MetroShapes.cpp
#        core/RNGraph.cpp
        optimization/Smooth.cpp
        optimization/Mixedlayout.cpp
        matching/Stationmatching.cpp
        main.cpp
        )
set(RESOURCES
        ../resources/ui.qrc
        matching/Manuelpath.cpp matching/Manuelpath.h)

add_executable(${Target_Name} ${SOURCES} ${RESOURCES} matching/AutoMatching.cpp matching/AutoMatching.h)


#----------------------------------------------------------------------
# path setting
#----------------------------------------------------------------------
if( APPLE )

    target_include_directories(${Target_Name}
            PUBLIC
            /opt/local/include/
            /usr/local/include/
            /usr/local/opt/qt/include/
            ${CMAKE_SOURCE_DIR}/Base/src/
            ${CMAKE_SOURCE_DIR}/Graph/src/
            ${CMAKE_SOURCE_DIR}/UI/src/
            )

    target_link_directories(${Target_Name}
            PUBLIC
            /opt/local/lib/
            /usr/local/lib/
            /usr/local/opt/qt/lib/
            ${CMAKE_BINARY_DIR}/Base/
            ${CMAKE_BINARY_DIR}/Graph/
            ${CMAKE_BINARY_DIR}/UI/
            )
elseif( UNIX )

    target_include_directories(${Target_Name}
            PUBLIC
            /opt/local/include/
            /usr/include/
            /usr/include/x86_64-linux-gnu/
            ${CMAKE_SOURCE_DIR}/include/
            ${QT5_PATH_DIR}/include/
            ${CMAKE_SOURCE_DIR}/Base/src/
            ${CMAKE_SOURCE_DIR}/Graph/src/
            ${CMAKE_SOURCE_DIR}/UI/src/
#            ${MAIN_PROJECT_ROOT}/Base/lib
#            ${MAIN_PROJECT_ROOT}/Graph/lib
#            ${MAIN_PROJECT_ROOT}/UI/lib
            ${THIRD_PARTY_INCLUDE}
            )

    link_directories( ${Target_Name}
            PUBLIC
            /opt/local/lib/
            /usr/lib/
            /usr/lib/x86_64-linux-gnu/
            ${CMAKE_SOURCE_DIR}/src/
            ${CMAKE_BINARY_DIR}/Base/
            ${CMAKE_BINARY_DIR}/Graph/
            ${CMAKE_BINARY_DIR}/UI/
            ${THIRD_PARTY_LIB}
            )
endif()

# strip the trailing space
string(STRIP "${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS)
SET_PROPERTY(TARGET ${Target_Name} PROPERTY DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

if( APPLE )
    set( CUSTOM_LIBRARIES "-lZuKaiBase -lZuKaiGraph -lZuKaiUI")
    target_link_libraries( ${Target_Name} PUBLIC ${CMAKE_CXX_FLAGS} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${CUSTOM_LIBRARIES})
    set_target_properties( ${Target_Name} PROPERTIES
            OUTPUT_NAME ${Target_Name}
            )
elseif( UNIX )
#    set( CUSTOM_LIBRARIES "ZuKaiBase ZuKaiGraph ZuKaiUI")
    target_link_libraries( ${Target_Name} PUBLIC ${CMAKE_CXX_FLAGS} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ZuKaiBase ZuKaiGraph ZuKaiUI)
    set_target_properties( ${Target_Name} PROPERTIES
            OUTPUT_NAME ${Target_Name}
            )
endif()



#----------------------------------------------------------------------
# QT
#----------------------------------------------------------------------
message(STATUS "-- Linking qt")
set(QT_REQUIRED_COMPONENTS Core Gui Widgets OpenGL)
find_package(Qt5 COMPONENTS ${QT_REQUIRED_COMPONENTS} REQUIRED)

# Add components to the target
foreach(component ${QT_REQUIRED_COMPONENTS})
    target_link_libraries(${Target_Name}
            PUBLIC
            Qt5::${component}
            )
endforeach()

# Enable automoc/uic/rcc
set_target_properties(${Target_Name}
        PROPERTIES
        AUTOMOC ON
        AUTOUIC ON
        AUTORCC ON
        )

#----------------------------------------------------------------------
# EIGEN
#----------------------------------------------------------------------
message(STATUS "-- Linking eigen")

find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(${Target_Name}
        PUBLIC
        Eigen3::Eigen
        )

#----------------------------------------------------------------------
# Boost
#----------------------------------------------------------------------
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED)

if(Boost_FOUND)
    message(STATUS "-- Linking boost")
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(${Target_Name} ${Boost_LIBRARIES})
endif()

#----------------------------------------------------------------------
# CGAL
#----------------------------------------------------------------------
#set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
#find_package(CGAL REQUIRED)

#if(CGAL_FOUND)
#    message(STATUS "-- Linking CGAL")
#    message(CGAL_LIBRARIES ${CGAL_LIBRARIES})
#    include_directories(${CGAL_INCLUDE_DIRS})
#    target_link_libraries(${Target_Name} PUBLIC ${CGAL_LIBRARIES})
#endif()

#----------------------------------------------------------------------
# COPY CONFIG AND DATA FOLDER
#----------------------------------------------------------------------

## Create directories
#file(COPY ${CMAKE_SOURCE_DIR}/${Target_Name}/config/ DESTINATION ${CMAKE_BINARY_DIR}/${Target_Name}/config/)

#add_custom_command(TARGET ${Target_Name}
#        POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E make_directory ${directory} "$<TARGET_FILE_DIR:${Target_Name}>/config/"
#        COMMAND ${CMAKE_COMMAND} -E make_directory ${directory} "$<TARGET_FILE_DIR:${Target_Name}>/data/"
#        COMMAND ${CMAKE_COMMAND} -E make_directory ${directory} "$<TARGET_FILE_DIR:${Target_Name}>/data/metro"
#        )
#if( APPLE OR UNIX )
#    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/data/metro"  ${CMAKE_BINARY_DIR}/bin/data/metro/ )
#elseif (WIN32)
#    #TODO
#endif()
