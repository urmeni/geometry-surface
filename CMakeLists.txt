cmake_minimum_required(VERSION 3.1...3.22)
project(geometry_surface)

if (APPLE)
    message("-------------- APPLE system found")
    INCLUDE_DIRECTORIES(/opt/homebrew/include)
else ()
    message("-------------- non APPLE system found")

    #    set(CGAL_DIR "/amuhome/bruasse.a/Public/Libs/CGAL-5.6/lib/cmake/CGAL")
    #    list (APPEND CMAKE_PREFIX_PATH "/opt/Qt-6.3/6.3.0/gcc_64")
    #    set (Qt5_DIR "/opt/Qt-6.3/6.3.0/gcc_64/lib/cmake/Qt6")
    #    set(Boost_USE_STATIC_LIBS OFF)
    #    set(Boost_USE_MULTITHREADED ON)
    #    set(Boost_USE_STATIC_RUNTIME OFF)
    #    set(Boost_NO_SYSTEM_PATHS TRUE)

    #    if (Boost_NO_SYSTEM_PATHS)
    #        set(BOOST_ROOT "/amuhome/bruasse.a/Public/Libs/boost_1_78")
    #        set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
    #        set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
    #    endif (Boost_NO_SYSTEM_PATHS)

    # ##### For AMU desktop, uncomment all previous section

endif (APPLE)

include_directories("./libs")

find_package(CGAL REQUIRED COMPONENTS Core)
find_package(CGAL REQUIRED COMPONENTS ImageIO)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)

# create a target per cppfile
set(SRCS
        code_courb.cpp
        mesh_functions.cpp
        courbures.cpp)

set(HEADERS
        mesh_functions.hpp
        starter_defs.h
        distrib.hpp
        courbures.h)

add_executable(code_courb ${SRCS} ${HEADERS})

set_target_properties(code_courb PROPERTIES
        CXX_STANDARD 17
        CXX_EXTENSIONS OFF)