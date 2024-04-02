cmake_minimum_required (VERSION 3.1)
project (JonkerVolgenant    VERSION 1.0
                            DESCRIPTION "Library for the Jonker-Volgenant algorithm to solve the assignement problem in a bipartite graph"
                            LANGUAGES C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

add_library (JonkerVolgenant SHARED
                src/assign2DCBasic.c
                src/interface.c
            )

set_target_properties(JonkerVolgenant PROPERTIES PUBLIC_HEADER src/assignAlgs2D.h)
set_property(TARGET JonkerVolgenant PROPERTY C_STANDARD 99)
set_property(TARGET JonkerVolgenant PROPERTY C_STANDARD_REQUIRED TRUE)
target_include_directories(JonkerVolgenant PRIVATE .)
install(TARGETS JonkerVolgenant DESTINATION lib)


