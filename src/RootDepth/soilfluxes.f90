subroutine SOILFLUXES(i, j, nzg, freedrain, dtll, slz, dz, soiltxt, smoiwtd, transp, transpdeep &
                      , smoi, wtd, rech, deeprech, precip, pet_s, et_s, runoff, flux, fdepth, qlat &
                      , qlatflux, qrf, qrfcorrect, flood, icefactor)!pppendepth)
  implicit none

  integer :: nzg, freedrain, nsoil, nsoil1, k, iwtd, kwtd, i, j
  real, dimension(nzg + 1) :: slz
  real, dimension(nzg) :: dz, vctr2, vctr4, vctr5, vctr6
  real, dimension(nzg) :: transp, smoi, kfmid, diffmid &
                          , aa, bb, cc, rr
  real, dimension(nzg + 1) :: vt3di, flux
  real, dimension(0:nzg + 1) :: qlatflux
  integer, dimension(2) :: soiltxt
  integer*1, dimension(nzg) :: icefactor
  real :: precip, runoff, pet_s, et_s, transpdeep, pppendepth
  real :: wgpmid, kfup, kfdw, hydcon, newwgp, smoiwtd, rech, deeprech, wtd, deeptemp &
          , fracliqwtd, wmid, wtdold, dzup, vt3dbdw, vt3dcdw, dtll, smoibot, dsmoi, icefac, ddw, dup &
          , smoisat, psisat, smoicp, fdepth, qlat, qrf, qrfcorrect, qgw, flood

  do k = 1, nzg
    vctr2(k) = 1./dz(k)
    vctr4(k) = 0.5*(slz(k) + slz(k + 1))
  end do
  do k = 2, nzg
    vctr5(k) = vctr4(k) - vctr4(k - 1)
    vctr6(k) = 1./vctr5(k)
  end do

  kfmid = 0.
  diffmid = 0.

  vt3di = 0.
  rech = 0.
  runoff = 0.

  qgw = qlat - qrf
!top boundary condition, infiltration + potential et from soil

  vt3di(nzg + 1) = (-precip + pet_s)*1.e-3 - flood

  if (freedrain .eq. 0) then
    do k = 1, nzg
      if (wtd .lt. slz(k)) exit
    end do
    iwtd = k
  else
    iwtd = 0
  end if

!if(i.eq.886.and.j.eq.564)write(6,*)'mirar soilfluxes',iwtd,wtd,qgw

  k = max(iwtd - 1, 1)
  qlatflux(k) = qlatflux(k) + qgw
!         do k = 2,nzg
  do k = max(iwtd - 1, 2), nzg
!gmmdiffusivity and conductivity at the interlayers
    wgpmid = smoi(k) + (smoi(k) - smoi(k - 1))*(slz(k) - vctr4(k))*vctr6(k)

    if (slz(k) .lt. -0.30) then
      nsoil = soiltxt(1)
    else
      nsoil = soiltxt(2)
    end if

    hydcon = Ksat(nsoil)*max(min(exp((slz(k) + 1.5)/fdepth), 1.), 0.1)
    smoisat = theta_sat(nsoil)*max(min(exp((slz(k) + 1.5)/fdepth), 1.), 0.1)
    psisat = slpots(nsoil)*min(max(exp(-(slz(k) + 1.5)/fdepth), 1.), 10.)

    wgpmid = min(wgpmid, smoisat)
!            icefac=fracliq(k)** (2. * slbs(nsoil) + 3.)
!            icefac=1.
    if (icefactor(k) .eq. 0) then
      icefac = 1.
    else
      icefac = 0.
    end if

    kfmid(k) = icefac*hydcon &
               *(wgpmid/smoisat)**(2.*slbs(nsoil) + 3.)
    diffmid(k) = -icefac*(hydcon*psisat*slbs(nsoil)/smoisat) &
                 *(wgpmid/smoisat)**(slbs(nsoil) + 2.)

!write(6,*)k,diffmid(k),kfdw,kfup,ddw,dup,smoi(k),smoi(k-1)
  end do
!calculate tridiagonal matrix elements
!       do k=2,nzg-1
  do k = max(iwtd - 2, 2), nzg
    aa(k) = diffmid(k)*vctr6(k)
    cc(k) = diffmid(k + 1)*vctr6(k + 1)
    bb(k) = -(aa(k) + cc(k) + dz(k)/dtll)
    rr(k) = -smoi(k)*dz(k)/dtll - kfmid(k + 1) + kfmid(k) + transp(k)/dtll
    if (k .eq. iwtd - 1) rr(k) = rr(k) - qgw/dtll
  end do
!boundary conditions

!top boundary
  aa(nzg) = diffmid(nzg)*vctr6(nzg)
  bb(nzg) = -aa(nzg) - dz(nzg)/dtll
  rr(nzg) = vt3di(nzg + 1)/dtll - smoi(nzg)*dz(nzg)/dtll + kfmid(nzg) + transp(nzg)/dtll
  if (iwtd - 1 .eq. nzg) rr(nzg) = rr(nzg) - qgw/dtll

!now bottom boundary condition
  IF (freedrain .ne. 1) then

    if (iwtd .le. 3) then
      aa(1) = 0.
      cc(1) = diffmid(2)*vctr6(2)
      bb(1) = -(aa(1) + cc(1) + dz(1)/dtll)
      rr(1) = -smoi(1)*dz(1)/dtll - kfmid(2) + transp(1)/dtll
      if (iwtd .le. 2) rr(1) = rr(1) - qgw/dtll

    else
      do k = 1, max(iwtd - 3, 1)
        aa(k) = 0.
        cc(k) = 0.
        bb(k) = 1.
        rr(k) = smoi(k)
      end do
    end if
  ELSE

!gmmgravitational drainage at the bottom
    nsoil = soiltxt(1)

    hydcon = Ksat(nsoil)*max(min(exp((slz(1) + 1.5)/fdepth), 1.), 0.1)
    smoisat = theta_sat(nsoil)*max(min(exp((slz(1) + 1.5)/fdepth), 1.), 0.1)

    kfmid(1) = hydcon &
               *(smoi(1)/smoisat)**(2.*slbs(nsoil) + 3.)

    aa(1) = 0.
    cc(1) = diffmid(2)*vctr6(2)
    bb(1) = -(cc(1) + dz(1)/dtll)
    rr(1) = -smoi(1)*dz(1)/dtll - kfmid(2) + kfmid(1) + transp(1)/dtll
  END IF
!if(i.eq.25.and.j.eq.30)write(6,*)'mirar wtd',wtd,theta_cp(nsoil),theta_sat(nsoil),nsoil
!if(i.eq.25.and.j.eq.30)write(6,*)'smoi antes',(smoi(k),k=1,nzg)

!solve tridiagonal system and update smoi
  call tridag(aa, bb, cc, rr, smoi, nzg)

!calculate the fluxes
  do k = 2, nzg
    vt3di(k) = (-aa(k)*(smoi(k) - smoi(k - 1)) - kfmid(k))*dtll
  end do
  vt3di(1) = 0.

!if(i.eq.25.and.j.eq.30)write(6,*)'fluxes antes',(vt3di(k),k=1,nzg)
!if(i.eq.25.and.j.eq.30)write(6,*)'capillarity',(-aa(k)*(smoi(k)-smoi(k-1)),k=2,nzg)
!if(i.eq.25.and.j.eq.30)write(6,*)'drainage',(-kfmid(k)*dtll,k=2,nzg)

! now check that soil moisture values are within bounds (theta_sat and theta_cp)
! if not, correct fluxes

!if(i.eq.2017.and.j.eq.751)write(6,*)'mirar',smoi(nzg),vt3di(nzg)
!if(i.eq.25.and.j.eq.30)write(6,*)'transp',(transp(k),k=1,nzg)

!if(i.eq.300.and.j.eq.300)write(6,*)'mirar antes',wtd
!if(i.eq.300.and.j.eq.200)write(6,*)'mirar flujo antes',vt3di(iwtd),vt3di(iwtd-1),vt3di(1),iwtd,wtd,smoi(1),smoi(2),smoiwtd
!if(i.eq.300.and.j.eq.300)write(6,*)smoi
  do k = 1, nzg

    if (slz(k) .lt. -0.30) then
      nsoil = soiltxt(1)
    else
      nsoil = soiltxt(2)
    end if

    smoisat = theta_sat(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)

!if(i.eq.54.and.j.eq.49)write(6,*)'soilflux',k,smoi(k),smoisat
    if (smoi(k) .gt. smoisat) then
      dsmoi = max((smoi(k) - smoisat)*dz(k), 0.)
      if (k .lt. nzg) then
        smoi(k + 1) = smoi(k + 1) + dsmoi*vctr2(k + 1)
        vt3di(k + 1) = vt3di(k + 1) + dsmoi
      else
        runoff = dsmoi
      end if
      smoi(k) = smoisat
    end if
  end do

!if(i.eq.54.and.j.eq.49)write(6,*)'mirar runoff',runoff,flood,precip,qlat,qrf
  do k = 1, nzg - 1

    if (slz(k) .lt. -0.30) then
      nsoil = soiltxt(1)
    else
      nsoil = soiltxt(2)
    end if
    smoicp = theta_cp(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)

    if (smoi(k) .lt. smoicp) then
      dsmoi = max((smoicp - smoi(k))*dz(k), 0.)
      smoi(k + 1) = smoi(k + 1) - dsmoi*vctr2(k + 1)
      vt3di(k + 1) = vt3di(k + 1) - dsmoi
      smoi(k) = smoicp
    end if
  end do

  k = nzg

  if (slz(k) .lt. -0.30) then
    nsoil = soiltxt(1)
  else
    nsoil = soiltxt(2)
  end if

  smoicp = theta_cp(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)

  if (smoi(k) .lt. smoicp) then
!first reduce soil evaporation from PET
    dsmoi = max((smoicp - smoi(k))*dz(k), 0.)
    if (vt3di(k + 1) .gt. dsmoi) then
      et_s = max(0., pet_s - dsmoi*1.e3)
      smoi(k) = smoicp
    else
      et_s = max(0., pet_s - max(vt3di(k + 1), 0.)*1.e3)
      smoi(k) = smoi(k) + max(vt3di(k + 1), 0.)/dz(k)
!take water from below
      dsmoi = max((smoicp - smoi(k))*dz(k), 0.)
      smoi(k - 1) = smoi(k - 1) - dsmoi*vctr2(k - 1)
      vt3di(k) = vt3di(k) + dsmoi
      smoi(k) = smoicp
!then go down all the way to the bottom
      do k = nzg - 1, 1, -1
        if (slz(k) .lt. -0.30) then
          nsoil = soiltxt(1)
        else
          nsoil = soiltxt(2)
        end if

        smoicp = theta_cp(nsoil)*max(min(exp((vctr4(k) + 1.5)/fdepth), 1.), 0.1)
        if (smoi(k) .lt. smoicp) then
          !take water from below
          dsmoi = max((smoicp - smoi(k))*dz(k), 0.)
          if (k .gt. 1) smoi(k - 1) = smoi(k - 1) - dsmoi*vctr2(k - 1)
          vt3di(k) = vt3di(k) + dsmoi
          smoi(k) = smoicp
        end if
      end do
    end if
  else
    et_s = pet_s
  end if

  if (vt3di(1) .gt. 0.) then
    qrfcorrect = -min(vt3di(1), max(qrf, 0.))
!       write(6,*)'too much qrf',i,j,qrf,qrfcorrect,vt3di(1),qlat,wtd
  else
    qrfcorrect = 0.
  end if

!if(i.eq.300.and.j.eq.200)write(6,*)'mirar flujo despues',vt3di(iwtd),vt3di(iwtd-1),vt3di(1),smoi(1),smoi(2),smoiwtd
!if(i.eq.300.and.j.eq.300)write(6,*)vt3di
!if(i.eq.300.and.j.eq.300)write(6,*)smoi

!save rain penetration depth

  flux = flux + vt3di

!     do k=1,nzg
!       if(vt3di(k).lt.-1.e-6)then
!            if(pppendepth.gt.slz(k))pppendepth=slz(k)
!            exit
!       endif
!     enddo
!if(i.eq.25.and.j.eq.30)write(6,*)'smoi despues',(smoi(k),k=1,nzg)
!if(i.eq.25.and.j.eq.30)write(6,*)'fluxes',(vt3di(k),k=1,nzg)

!if(i.eq.2017.and.j.eq.751)write(6,*)'mirar1',smoiwtd
!if(i.eq.251.and.j.eq.369)write(6,*)'mirar',smoi(14),smoi(13),smoi(12),vt3di(15),vt3di(14),vt3di(13),vt3di(12),vt3di(11)
  IF (freedrain .eq. 1) then

!accumulate gravitational drainage
    rech = vt3di(1)
!smoiwtd is now the bucket of water at the bottom, to save the water and put it later into the rivers
    smoiwtd = smoiwtd - vt3di(1)
  END IF

end subroutine soilfluxes
