cmake_minimum_required(VERSION 3.1...3.14)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(basis LANGUAGES CXX)

set (TARGET basis)

set(SOURCES jjjpn_operator.cpp jjjpn_scheme.cpp
            jjjt_operator.cpp jjjt_scheme.cpp
            lsjt_operator.cpp lsjt_scheme.cpp
            m_scheme.cpp
            many_body.cpp
            nlj_operator.cpp nlj_orbital.cpp
            proton_neutron.cpp)

function(add_headers VAR)
  set(headers ${${VAR}})
  foreach (header ${ARGN})
    set(headers ${headers} ${CMAKE_CURRENT_SOURCE_DIR}/${header})
  endforeach()
  set(${VAR} ${headers} PARENT_SCOPE)
endfunction()

add_headers(HEADERS basis.h degenerate.h
                    hypersector.h
                    jjjpn_operator.h jjjpn_scheme.h
                    jjjt_operator.h jjjt_scheme.h
                    jt_operator.h
                    lsjt_operator.h lsjt_scheme.h
                    m_scheme.h
                    many_body.h
                    nlj_operator.h nlj_orbital.h
                    operator.h
                    proton_neutron.h)

add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

target_link_libraries(${TARGET} PUBLIC m)
target_link_libraries(${TARGET} PRIVATE am mcutils fmt::fmt)

find_package(Eigen3 REQUIRED)
target_link_libraries(${TARGET} PUBLIC Eigen3::Eigen)

find_package(GSL REQUIRED)
target_link_libraries(${TARGET} PUBLIC GSL::gsl GSL::gslcblas)

target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
