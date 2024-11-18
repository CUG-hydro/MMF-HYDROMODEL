subroutine GW2RIVER(imax, js, je, nzg, slz, deltat, soiltxt, landmask, wtd, maxdepth, riverdepth, width, length, area, fdepth, qrf)
  implicit none
  integer :: i, j, imax, js, je, nsoil, k, iwtd, nzg
  real, dimension(nzg + 1) :: slz
  integer, dimension(2, imax, js:je) :: soiltxt
  integer, dimension(imax, js:je) :: landmask
  real, dimension(imax, js:je) :: wtd, maxdepth, riverdepth, width, length, area, qrf, fdepth
  real :: riversurface, deltat, soilwatercap, rcond, rdepth, hydcon

  soilwatercap = 0.
  qrf = 0.

  do j = js + 1, je - 1
    do i = 2, imax - 1
      if (landmask(i, j) .eq. 0 .or. width(i, j) .eq. 0.) cycle
      rdepth = max(riverdepth(i, j), 0.)

      nsoil = soiltxt(2, i, j)
      riversurface = -(maxdepth(i, j) - rdepth)
      if (riversurface .ge. 0.) cycle      !this just in case...

      hydcon = Ksat(nsoil)*max(min(exp((-maxdepth(i, j) + 1.5)/fdepth(i, j)), 1.), 0.1)
      rcond = width(i, j)*length(i, j)*hydcon

      if (wtd(i, j) .gt. riversurface) then

        qrf(i, j) = rcond*(wtd(i, j) - riversurface)*(deltat/area(i, j))

!limit it to prevent sudden drops , lets say 50mm per day 0.05/86400.
!                  qrf(i,j)=min(qrf(i,j),deltat*0.05/86400.)

      elseif (wtd(i, j) .gt. -maxdepth(i, j)) then   !water table connected to the river, even though below river surface

        soilwatercap = -rcond*(wtd(i, j) - riversurface)*(deltat/area(i, j))
!                  soilwatercap=min(soilwatercap,deltat*0.05/86400.)
        qrf(i, j) = -max(min(soilwatercap, riverdepth(i, j)), 0.)*min(width(i, j)*length(i, j)/area(i, j), 1.)

      else
!water table below river bed, disconnected from the river. No rcond use, just
!infiltration. Assume it occurs at the Ksat rate and water goes directly to the
!water table.
        qrf(i, j) = -max(min(Ksat(nsoil)*deltat, rdepth), 0.)*min(width(i, j)*length(i, j)/area(i, j), 1.)
      end if

    end do
  end do
end subroutine gw2river
