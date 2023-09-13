# This file is part of CPCM-X.
# SPDX-Identifier: LGPL-3.0-or-later
#
# CPCM-X is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CPCM-X is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CPCM-X.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.14)

project(
  "cpcmx"
  LANGUAGES "Fortran" "C"
  VERSION "1.0.0"
  DESCRIPTION "Extended conduction-like polarizable continuum model"
)
enable_testing()

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
add_subdirectory("config")

# Dependencies
if(NOT TARGET "numsa::numsa")
  find_package("numsa" REQUIRED)
endif()

if(NOT TARGET "mctc-lib::mctc-lib")
   find_package("mctc-lib" REQUIRED)
endif()

if(NOT TARGET "toml-f::toml-f")
   find_package("toml-f" REQUIRED)
endif()

if(NOT TARGET "test-drive::test-drive")
   find_package("test-drive" REQUIRED)
endif()
find_package("LAPACK" REQUIRED)
find_package("BLAS" REQUIRED)

# Collect source of the project
set(srcs)
add_subdirectory("src")

# CPCM-X library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
target_compile_definitions(
  "${PROJECT_NAME}-lib"
  PRIVATE
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PRIVATE
  "numsa::numsa"
  "mctc-lib::mctc-lib"
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

# Add example application
add_subdirectory("app")

# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${PROJECT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)

# add the testsuite
enable_testing()
add_subdirectory("test")
