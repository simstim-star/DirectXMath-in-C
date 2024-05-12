cmake_minimum_required (VERSION 3.20)

project(XMathC DESCRIPTION "DirectXMath SIMD C math library")

set(CMAKE_C_STANDARD 11)
set(LIBRARY_HEADERS
    Inc/DirectXMathC.h
    Inc/DirectXMathCConvert.inl
    Inc/DirectXMathCMatrix.inl)

option(XM_USE_TEST "Enable Tests" OFF)
option(XM_NO_INTRINSICS "Disable intrinsics" OFF)

if(XM_NO_INTRINSICS)
  add_compile_definitions(_XM_NO_INTRINSICS_)
endif()

if(MSVC)
  add_definitions(-D_WINDOWS -D_USRDLL)

  if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DNDEBUG)
    add_compile_options(/W3 /Ox /Gy /Oi /TC)
    foreach(flag_var
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
      string(REGEX REPLACE "/RTC(su|[1su])" "" ${flag_var} "${${flag_var}}")
    endforeach(flag_var)
  endif()
else()
  add_compile_options(-Wall)

  if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(-O3)
  endif()
endif()

add_library(XMathC INTERFACE ${LIBRARY_HEADERS})

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Inc>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/directxmath>)

source_group("Headers" FILES ${LIBRARY_HEADERS})

# Test Configuration
if(XM_USE_TEST)
  include(CTest)
  enable_testing()
  add_subdirectory(test)
endif()
