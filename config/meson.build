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
elif fc_id == 'intel' or fc_id == 'intel-llvm'
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

if get_option('openmp')
  omp_dep = dependency('openmp')
  lib_deps += omp_dep
endif

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel' or fc_id == 'intel-llvm'
    lapack_vendor = 'mkl'
  endif
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel' or fc_id == 'intel-llvm'
    mkl_dep += fc.find_library('mkl_intel_lp64')
    if get_option('openmp')
      mkl_dep += fc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += fc.find_library('mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += fc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += fc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += fc.find_library('mkl_core')
  lib_deps += mkl_dep

elif lapack_vendor == 'mkl-rt'
  mkl_dep = fc.find_library('mkl_rt')
  lib_deps += mkl_dep

elif lapack_vendor == 'openblas'
  openblas_dep = dependency('openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = fc.find_library('openblas')
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
    endif
    lib_deps += lapack_dep
  endif
elif lapack_vendor == 'custom'
  foreach lib: get_option('custom_libraries')
    lib_deps += fc.find_library(lib)
  endforeach
else
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack')
  endif
  lib_deps += lapack_dep
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas')
  endif
  lib_deps += blas_dep
endif

numsa_dep = dependency(
  'numsa',
  fallback: ['numsa', 'numsa_dep'],
  default_options: ['default_library=static'],
)
lib_deps += numsa_dep


tomlf_prj = subproject(
  'toml-f',
  version: '>=0.2',
  default_options: [
    'default_library=static',
  ],
)
tomlf_dep = tomlf_prj.get_variable('tomlf_dep')
lib_deps += tomlf_dep
