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

if(NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

set(
  tests
  "calc"
  "io"
)
set(
  test-srcs
  "main.f90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND test-srcs "test_${t}.f90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
  "test-drive::test-drive"
  "mctc-lib::mctc-lib"
)
set_target_properties(
  "${PROJECT_NAME}-tester"
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/include")
  make_directory("${CMAKE_CURRENT_BINARY_DIR}/include")
endif()

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
  set_tests_properties(
    "${PROJECT_NAME}/${t}"
    PROPERTIES
    ENVIRONMENT CPXHOME=${PROJECT_SOURCE_DIR}
  )
endforeach()


set(CPX_EXE "${PROJECT_BINARY_DIR}/app/cpx")

add_test("EXE_Argparser_print_help" ${CPX_EXE} --help)

set_tests_properties(
  PROPERTIES
  ENVIRONMENT CPXHOME=${PROJECT_SOURCE_DIR}
)
