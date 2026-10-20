program project1

implicit none

integer, parameter :: dp = kind(1.0d0) ! double precision
integer, parameter :: jmax=5000 ! need jmax to be a known value at compile time cause we use it to declare array sizes later

! define constants from slides
real(dp), parameter :: msol = 1.989d33 ! g
real(dp), parameter :: kpc = 3.086d21 ! cm
real(dp), parameter :: G = 6.672d-8 ! cgs
real(dp), parameter :: k_B = 1.380649d-16 ! erg/K
real(dp), parameter :: m_p = 1.6726219d-24 ! g
real(dp), parameter :: mu = 0.61_dp

real(dp) :: r(jmax), rr(jmax)
real(dp) :: rmin, rmax, dr
integer  :: j, iu

! set the grid

rmin = 0.0_dp*kpc
rmax = 3000.0_dp*kpc
dr = (rmax-rmin)/real(jmax-1, dp) ! define grid spacing between nodes: the total length divided by the number of intervals
do j=1,jmax
    r(j)  = rmin + real(j-1, dp)*dr ! nodes r_j
    rr(j) = r(j) + 0.5_dp*dr ! half-point grid, shifted by half a step, r_{j+1/2}, defined for every j
enddo

! check that the grid ends at rmax
if (abs(r(jmax) - rmax) > 1.0d-10*rmax) stop 'grid does not end at rmax' ! stop code if grid does not end at rmax

! write to file (in kpc)
open(newunit=iu, file='grid.dat', status='replace')
do j = 1, jmax
   write(iu,'(2(1pe14.6))') r(j)/kpc, rr(j)/kpc
end do
close(iu)

end program project1
