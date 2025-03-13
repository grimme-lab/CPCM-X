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

project(
  'cpx',
  'fortran', 'c',
  version: '1.1.0',
  meson_version: '>=0.60.0',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=static',
  ],
)

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

lib_deps = []

if get_option('openmp')
  omp_dep = dependency('openmp', required: fc_id != 'intel' and fc_id != 'intel-llvm' and fc_id != 'nvidia_hpc')
  if not omp_dep.found()
    if fc_id == 'intel' or fc_id != 'intel-llvm'
      message('Using -qopenmp to use OpenMP with Intel compilers')
      omp_dep = declare_dependency(
        compile_args: '-qopenmp',
        link_args: '-qopenmp',
      )
    else
      message('Using -mp to use OpenMP with NVHPC compilers')
      omp_dep = declare_dependency(
        compile_args: '-mp',
        link_args: '-mp',
      )
    endif
  endif
  lib_deps += omp_dep
endif

lib_deps += dependency('threads')

subdir('config')

srcs = []
subdir('src')

cpx_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  install: true,
)
cpx_dep = declare_dependency(
  link_with: cpx_lib,
  dependencies: lib_deps,
)

inc_dir=include_directories('.')
subdir('app')

module_id = meson.project_name() / fc_id + '-' + fc.version()
meson.add_install_script(
  find_program(files('config'/'install-mod.py')),
  get_option('includedir') / module_id,
)

pkg = import('pkgconfig')
pkg.generate(
  cpx_lib,
  subdirs: ['', module_id],
)

cpxenv = environment()
cpxenv.set('CPXHOME', meson.current_source_dir())

subdir('test')
