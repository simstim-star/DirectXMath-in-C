cmake_minimum_required (VERSION 3.20)

project(XMathC DESCRIPTION "DirectXMath SIMD C math library")

set(XMathC_VERSION 0.0.1)

set(CMAKE_C_STANDARD 11)
set(LIBRARY_HEADERS
    Inc/DirectXMathC.h
    Inc/DirectXMathCConvert.inl
    Inc/DirectXMathCMatrix.inl
    Inc/DirectXMathCVector.inl
    Inc/DirectXMathCMisc.inl
    Inc/DirectXCollisionC.h
    Inc/DirectXCollisionC.inl)

option(XM_USE_TEST "Enable Tests" OFF)
option(XM_NO_INTRINSICS "Disable intrinsics" OFF)
option(USE_VECTORCALL "Enable __vectorcall calling convention" ON)


if(XM_NO_INTRINSICS)
  add_compile_definitions(_XM_NO_INTRINSICS_)
endif()

if(USE_VECTORCALL)
    add_compile_definitions(_XM_VECTORCALL_)
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
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/xmathc>)

source_group("Headers" FILES ${LIBRARY_HEADERS})


# https://stackoverflow.com/questions/47718485/install-and-export-interface-only-library-cmake
#--- Package
include(CMakePackageConfigHelpers)

string(TOLOWER ${PROJECT_NAME} PACKAGE_NAME)

write_basic_package_version_file(
  ${PACKAGE_NAME}-config-version.cmake
  VERSION ${XMathC_VERSION}
  COMPATIBILITY AnyNewerVersion
  ARCH_INDEPENDENT)

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PACKAGE_NAME})

install(EXPORT ${PROJECT_NAME}-targets
  FILE ${PROJECT_NAME}-targets.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PACKAGE_NAME})

install(FILES ${LIBRARY_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/xmathc)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PACKAGE_NAME})



# Test Configuration
if(XM_USE_TEST)
  include(CTest)
  enable_testing()
  add_subdirectory(test)
endif()
