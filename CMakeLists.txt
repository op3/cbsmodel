cmake_minimum_required(VERSION 3.13)

project(cbsmodel
    LANGUAGES CXX
    VERSION 1.2.0
    DESCRIPTION "This program calculates excitation energies E0 and E2 transition strengths within the 'Confined Beta Soft' nuclear structure model as described in [N. Pietralla and O. M. Gorbachenko, Phys. Rev. C 70, 011304(R) (2004)]. Parameter fits are possible.")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_executable(
    cbsmodel
    src/Main.cpp
    src/CBSNucleus.cpp
)

set(READLINE ON CACHE BOOL "Enable GNU Readline")
if(READLINE)
    find_package(Readline)
    set(HAVE_LIBREADLINE ${READLINE_FOUND})
    if(READLINE_FOUND)
        target_link_libraries(
            cbsmodel
            PRIVATE
            readline
        )
        target_compile_definitions(
            cbsmodel
            PRIVATE
            HAVE_LIBREADLINE=1
        )
    else()
        message(FATAL_ERROR "GNU readline was not found.\n"
            "Use -DREADLINE=no to explicitely disable it.\n")
    endif()
endif()

find_package(GSL REQUIRED)

target_link_libraries(
    cbsmodel
    PUBLIC
    m
    GSL::gsl
)

