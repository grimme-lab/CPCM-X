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

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
   "${dir}/bonding.f90"
   "${dir}/broyden.f90"
   "${dir}/cpxcalc.f90"
   "${dir}/crs.f90"
   "${dir}/data.f90"
   "${dir}/ElementParameters.f90"
   "${dir}/globals.f90"
   "${dir}/initialize.f90"
   "${dir}/internaldb.f90"
   "${dir}/isodens.f90"
   "${dir}/lapack.f90"
   "${dir}/meson.build"
   "${dir}/pr.f90"
   "${dir}/profile.f90"
   "${dir}/qc_calc.f90"
   "${dir}/sac.f90"
   "${dir}/sigma_av.f90"
   "${dir}/smd.f90"
   "${dir}/sort.f90"
   "${dir}/timer.f90"
   "${dir}/type.f90"
)

set(srcs "${srcs}" PARENT_SCOPE) 
