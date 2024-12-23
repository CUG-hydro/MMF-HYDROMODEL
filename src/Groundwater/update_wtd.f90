subroutine update_wtd(nzg, slz, dz, wtd, qspring, totwater, theta, theta_eq, soiltextures, smoiwtd)
  implicit none
  integer :: nzg, iwtd, kwtd, nsoil, nsoil1, k, k1
  real, dimension(nzg + 1) :: slz
  real, dimension(nzg) :: dz
  real :: wtd, qspring, wtdold, totwater, smoiwtd, maxwatup, maxwatdw, wgpmid, syielddw, dzup, tempk, fracliq, smoi_eqwtd
  real, dimension(nzg) :: theta, theta_eq
  real :: sy                                                                    ! specific yield

  integer, dimension(2) :: soiltextures
  integer, dimension(nzg) :: soiltxt

  where (slz .lt. -0.3)
    soiltxt = soiltextures(1)
  elsewhere
    soiltxt = soiltextures(2)
  end where

  iwtd = 1
!case 1: totwater > 0 (water table going up):
  IF (totwater .gt. 0.) then
    if (wtd .ge. slz(1)) then

      do k = 2, nzg
        if (wtd .lt. slz(k)) exit
      end do
      iwtd = k            ! 所在下一层，非饱和层
      kwtd = iwtd - 1     ! 地下水所在层，
      nsoil = soiltxt(kwtd)

!max water that fits in the layer
      maxwatup = dz(kwtd)*(theta_sat(nsoil) - theta(kwtd))  ! 剩余填充能力

      if (totwater .le. maxwatup) then
        theta(kwtd) = theta(kwtd) + totwater/dz(kwtd)
        theta(kwtd) = min(theta(kwtd), theta_sat(nsoil))

        if (theta(kwtd) .gt. theta_eq(kwtd)) then
          ! TODO: 考证这个公式如何来的
          ! dz(kwtd) = slz(iwtd) - slz(kwtd) ! 非饱和 - 饱和，
          sy = theta_sat(nsoil) - theta_eq(kwtd)
          wtd = (theta(kwtd)*dz(kwtd) - theta_eq(kwtd)*slz(iwtd) + theta_sat(nsoil)*slz(kwtd))/sy
          wtd = min(wtd, slz(iwtd))
        end if
        totwater = 0.
      else  !water enough to saturate the layer
        theta(kwtd) = theta_sat(nsoil)
        totwater = totwater - maxwatup
        k1 = iwtd
        do k = k1, nzg + 1
          wtd = slz(k)
          iwtd = k + 1
          if (k .eq. nzg + 1) exit
          nsoil = soiltxt(k)
          maxwatup = dz(k)*(theta_sat(nsoil) - theta(k))
          if (totwater .le. maxwatup) then
            theta(k) = theta(k) + totwater/dz(k)
            theta(k) = min(theta(k), theta_sat(nsoil))
            if (theta(k) .gt. theta_eq(k)) then
              sy = theta_sat(nsoil) - theta_eq(k)
              wtd = (theta(k)*dz(k) - theta_eq(k)*slz(iwtd) + theta_sat(nsoil)*slz(k))/sy
              wtd = min(wtd, slz(iwtd))
            end if

            totwater = 0.
            exit
          else
            theta(k) = theta_sat(nsoil)
            totwater = totwater - maxwatup
          end if

        end do
      end if

    elseif (wtd .ge. slz(1) - dz(1)) then ! wtd below bottom of soil model

      nsoil = soiltxt(1)
      maxwatup = (theta_sat(nsoil) - smoiwtd)*dz(1)

      if (totwater .le. maxwatup) then
        smoi_eqwtd = theta_sat(nsoil)*(slpots(nsoil)/ &
                                       (slpots(nsoil) - dz(1)))**(1./slbs(nsoil))
        smoi_eqwtd = max(smoi_eqwtd, theta_cp(nsoil))

        smoiwtd = smoiwtd + totwater/dz(1)
        smoiwtd = min(smoiwtd, theta_sat(nsoil))
        if (smoiwtd .gt. smoi_eqwtd) wtd = min((smoiwtd*dz(1) &
                                                - smoi_eqwtd*slz(1) + theta_sat(nsoil)*(slz(1) - dz(1)))/ &
                                               (theta_sat(nsoil) - smoi_eqwtd), slz(1))
        totwater = 0.
      else
        smoiwtd = theta_sat(nsoil)
        totwater = totwater - maxwatup
        do k = 1, nzg + 1
          wtd = slz(k)
          iwtd = k + 1
          if (k .eq. nzg + 1) exit
          nsoil = soiltxt(k)
          maxwatup = dz(k)*(theta_sat(nsoil) - theta(k))
          if (totwater .le. maxwatup) then
            theta(k) = min(theta(k) + totwater/dz(k), theta_sat(nsoil))
            if (theta(k) .gt. theta_eq(k)) wtd = min((theta(k)*dz(k) &
                                                      - theta_eq(k)*slz(iwtd) + theta_sat(nsoil)*slz(k))/ &
                                                     (theta_sat(nsoil) - theta_eq(k)), slz(iwtd))
            totwater = 0.
            exit
          else
            theta(k) = theta_sat(nsoil)
            totwater = totwater - maxwatup
          end if
        end do
      end if

!deep water table
    else
      nsoil = soiltxt(1)
      maxwatup = (theta_sat(nsoil) - smoiwtd)*(slz(1) - dz(1) - wtd)
      if (totwater .le. maxwatup) then
        wtd = wtd + totwater/(theta_sat(nsoil) - smoiwtd)
        totwater = 0.
      else
        totwater = totwater - maxwatup
        wtd = slz(1) - dz(1)
        maxwatup = (theta_sat(nsoil) - smoiwtd)*dz(1)
        if (totwater .le. maxwatup) then
          smoi_eqwtd = theta_sat(nsoil)*(slpots(nsoil)/ &
                                         (slpots(nsoil) - dz(1)))**(1./slbs(nsoil))
          smoi_eqwtd = max(smoi_eqwtd, theta_cp(nsoil))

          smoiwtd = smoiwtd + totwater/dz(1)
          smoiwtd = min(smoiwtd, theta_sat(nsoil))
          wtd = (smoiwtd*dz(1) &
                 - smoi_eqwtd*slz(1) + theta_sat(nsoil)*(slz(1) - dz(1)))/ &
                (theta_sat(nsoil) - smoi_eqwtd)
          totwater = 0.
        else
          smoiwtd = theta_sat(nsoil)
          totwater = totwater - maxwatup
          do k = 1, nzg + 1
            wtd = slz(k)
            iwtd = k + 1
            if (k .eq. nzg + 1) exit
            nsoil = soiltxt(k)
            maxwatup = dz(k)*(theta_sat(nsoil) - theta(k))

            if (totwater .le. maxwatup) then
              theta(k) = theta(k) + totwater/dz(k)
              theta(k) = min(theta(k), theta_sat(nsoil))
              if (theta(k) .gt. theta_eq(k)) wtd = (theta(k)*dz(k) &
                                                    - theta_eq(k)*slz(iwtd) + theta_sat(nsoil)*slz(k))/ &
                                                   (theta_sat(nsoil) - theta_eq(k))
              totwater = 0.
              exit
            else
              theta(k) = theta_sat(nsoil)
              totwater = totwater - maxwatup
            end if
          end do
        end if
      end if
    end if

!water springing at the surface
    qspring = totwater

!case 2: totwater < 0 (water table going down):
  ELSEIF (totwater .lt. 0.) then

    if (wtd .ge. slz(1)) then !wtd in the resolved layers

      do k = 2, nzg
        if (wtd .lt. slz(k)) exit
      end do
      iwtd = k

      k1 = iwtd - 1
      do kwtd = k1, 1, -1

        nsoil = soiltxt(kwtd)

!max water that the layer can yield
        maxwatdw = dz(kwtd)*(theta(kwtd) - theta_eq(kwtd))

        if (-totwater .le. maxwatdw) then
          theta(kwtd) = theta(kwtd) + totwater/dz(kwtd)
          if (theta(kwtd) .gt. theta_eq(kwtd)) then
            wtd = (theta(kwtd)*dz(kwtd) &
                   - theta_eq(kwtd)*slz(iwtd) + theta_sat(nsoil)*slz(kwtd))/ &
                  (theta_sat(nsoil) - theta_eq(kwtd))
          else
            wtd = slz(kwtd)
            iwtd = iwtd - 1
          end if
          totwater = 0.
          exit
        else
          wtd = slz(kwtd)
          iwtd = iwtd - 1
          if (maxwatdw .ge. 0.) then
            theta(kwtd) = theta_eq(kwtd)
            totwater = totwater + maxwatdw
          end if
        end if

      end do

      if (iwtd .eq. 1 .and. totwater .lt. 0.) then
        nsoil = soiltxt(1)
        smoi_eqwtd = theta_sat(nsoil)*(slpots(nsoil)/ &
                                       (slpots(nsoil) - dz(1)))**(1./slbs(nsoil))
        smoi_eqwtd = max(smoi_eqwtd, theta_cp(nsoil))

        maxwatdw = dz(1)*(smoiwtd - smoi_eqwtd)

        if (-totwater .le. maxwatdw) then

          smoiwtd = smoiwtd + totwater/dz(1)
          wtd = max((smoiwtd*dz(1) &
                     - smoi_eqwtd*slz(1) + theta_sat(nsoil)*(slz(1) - dz(1)))/ &
                    (theta_sat(nsoil) - smoi_eqwtd), slz(1) - dz(1))

        else
          wtd = slz(1) - dz(1)
          smoiwtd = smoiwtd + totwater/dz(1)
!and now even further down
          dzup = (smoi_eqwtd - smoiwtd)*dz(1)/(theta_sat(nsoil) - smoi_eqwtd)
          wtd = wtd - dzup
          smoiwtd = smoi_eqwtd
        end if
      end if

    elseif (wtd .ge. slz(1) - dz(1)) then

!if wtd was already below the bottom of the resolved soil crust
      nsoil = soiltxt(1)
      smoi_eqwtd = theta_sat(nsoil)*(slpots(nsoil)/ &
                                     (slpots(nsoil) - dz(1)))**(1./slbs(nsoil))
      smoi_eqwtd = max(smoi_eqwtd, theta_cp(nsoil))

      maxwatdw = dz(1)*(smoiwtd - smoi_eqwtd)

      if (-totwater .le. maxwatdw) then

        smoiwtd = smoiwtd + totwater/dz(1)
        wtd = max((smoiwtd*dz(1) &
                   - smoi_eqwtd*slz(1) + theta_sat(nsoil)*(slz(1) - dz(1)))/ &
                  (theta_sat(nsoil) - smoi_eqwtd), slz(1) - dz(1))
      else
        wtd = slz(1) - dz(1)
        smoiwtd = smoiwtd + totwater/dz(1)
!and now even further down
        dzup = (smoi_eqwtd - smoiwtd)*dz(1)/(theta_sat(nsoil) - smoi_eqwtd)
        wtd = wtd - dzup
        smoiwtd = smoi_eqwtd
      end if

    else
!gmmequilibrium soil moisture content
      nsoil = soiltxt(1)
      wgpmid = theta_sat(nsoil)*(slpots(nsoil)/ &
                                 (slpots(nsoil) - (slz(1) - wtd)))**(1./slbs(nsoil))
      wgpmid = max(wgpmid, theta_cp(nsoil))
      syielddw = theta_sat(nsoil) - wgpmid
      wtdold = wtd
      wtd = wtdold + totwater/syielddw
!update wtdwgp
      smoiwtd = (smoiwtd*(slz(1) - wtdold) + wgpmid*(wtdold - wtd))/(slz(1) - wtd)
    end if
    qspring = 0.
  END IF

end subroutine update_wtd
