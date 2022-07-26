cmake_minimum_required( VERSION 3.8 )
PROJECT( ZUKAI )
set( Target_Name ZuKai )

#----------------------------------------------------------------------
# cmake setting
#----------------------------------------------------------------------
message( "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR} )
message( "CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR} )
## set output folder
## set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
#set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
#set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
##set( SEPARATE_BUILD ON)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory( Base )
add_subdirectory( Graph )
add_subdirectory( UI )
#add_subdirectory( FileIO )
#add_subdirectory( FocusContext ) //Commented Out Problems with naming 
#add_subdirectory( FocusContext++ )
#add_subdirectory( Artistic )
add_subdirectory( MetroShapes )

# Copy directories
file(COPY ${CMAKE_SOURCE_DIR}/config/ DESTINATION ${CMAKE_BINARY_DIR}/config/)
file(COPY ${CMAKE_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/data/)
