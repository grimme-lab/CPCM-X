# This file is part of COSMO-X.
# SPDX-Identifier: LGPL-3.0-or-later
#
# COSMO-X is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# COSMO-X is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with COSMO-X.  If not, see <https://www.gnu.org/licenses/>.

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

omp_dep = dependency('openmp')
lib_deps += omp_dep

lapack_dep = dependency('mkl', 'lapack')
lib_deps += lapack_dep

numsa_dep = dependency(
  'numsa',
  fallback: ['numsa', 'numsa_dep'],
  default_options: ['default_library=static'],
)
lib_deps += numsa_dep
