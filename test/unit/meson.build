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

# Find test-drive dependency for testing
testdrive_dep = dependency(
  'test-drive',
  version: '>=0.3',
  fallback: ['test-drive', 'testdrive_dep'],
  required: not meson.is_subproject(),
  default_options: [
    'default_library=static',
    'openmp=false',
  ],
)
# If we do not find test-drive and are a subproject, we just skip testing
if not testdrive_dep.found()
  subdir_done()
endif

tests = [
  'io',
  'calc',
]

test_srcs = files(
  'main.f90',
)

foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [cpx_dep, testdrive_dep],
  link_language: 'fortran',
)

foreach t : tests
  test(t, tester, args: t, suite: 'unit', env: cpxenv)
endforeach

# some very basic checks to see if the executable reacts at all
#test('Argparser print version', cpx_exe, args: '--version', env: cpxenv)
test('Argparser print help', cpx_exe, args: '--help', env: cpxenv)
#test('Argparser print license', cpx_exe, args: '--license', env: cpxenv)
#test('Argparser no arguments', cpx_exe, should_fail: true, env: cpxenv)
