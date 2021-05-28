module sdm

   contains

      !> Example implementation to calculate surface area for a molecule input
      subroutine get_surface_area(species, symbols, coord, probe, surface, dsdr)
      use mctc_env, only : wp
      use numsa, only : surface_integrator, new_surface_integrator, get_vdw_rad_cosmo, grid_size
      !> Unique chemical species in the input structure, shape: [nat]
      integer, intent(in) :: species(:)
      !> Element symbol for each chemical species, shape: [nsp]
      character(len=*), intent(in) :: symbols(:)
      !> Cartesian coordinates in Bohr, shape: [3, nat]
      real(wp), intent(in) :: coord(:, :)
      !> Probe radius for surface area integration in Bohr
      real(wp), intent(in) :: probe
      !> Accessible surface area in Bohr², shape: [nat]
      real(wp), intent(out) :: surface(:)
      !> Derivative of surface area w.r.t. atomic displacements, shape: [3, nat, nat]
      real(wp), intent(out) :: dsdr(:, :, :)

      type(surface_integrator) :: sasa
      real(wp), allocatable :: rad(:)

      rad = get_vdw_rad_cosmo(symbols)
      call new_surface_integrator(sasa, species, rad, probe, grid_size(8))
      call sasa%get_surface(species, coord, surface, dsdr)

end subroutine get_surface_area
   



end module sdm
