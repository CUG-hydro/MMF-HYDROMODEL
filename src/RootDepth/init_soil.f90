subroutine init_soil_depth_clm(nzg, slz, dz)
  integer :: nzg, k, kk
  real, dimension(nzg + 1) :: slz, slz2
  real, dimension(nzg) :: dz, dz2, vctr4

  do k = 1, nzg
    vctr4(k) = 0.025*(exp(0.5*(float(k) - 0.5)) - 1.)
  end do

!write(6,*)'soil nodes',(-vctr4(k),k=nzg,1,-1)
  do k = 2, nzg - 1
    dz2(k) = 0.5*(vctr4(k + 1) - vctr4(k - 1))
  end do
  dz2(1) = 0.5*(vctr4(1) + vctr4(2))
  dz2(nzg) = vctr4(nzg) - vctr4(nzg - 1)

  do k = 1, nzg
    slz2(k) = 0.5*(vctr4(k) + vctr4(k + 1))
  end do

  slz2(nzg) = vctr4(nzg) + 0.5*dz2(nzg)

  do k = 1, nzg
    kk = nzg - k + 1
    slz(k) = -slz2(kk)
    dz(k) = dz2(kk)
  end do

  slz(nzg + 1) = 0.
end subroutine init_soil_depth_clm

subroutine init_soil_depth(nzg, slz, dz)
  integer :: nzg, k, kk
  real, dimension(nzg + 1) :: slz
  real, dimension(nzg) :: dz
  real, dimension(40) :: dz2
  data dz2/.1, .1, .1, .1, .1, .2, .2, .2, .2, .2, .3, .3, .3, .3, .4, .4 &
    , .4, .5, .5, .6, .7, .7, .8, .9, 1., 1., 1.2, 1.2, 1.5, 1.5, 2., 2. &
    , 3., 6., 11., 20., 50., 100., 250., 540./

!5.,10.,25.,50.,100.,200.,500.,1000./

  slz(nzg + 1) = 0.
  do k = nzg, 1, -1
    dz(k) = dz2(nzg - k + 1)
    slz(k) = slz(k + 1) - dz(k)
  end do

end subroutine init_soil_depth

FUNCTION khyd(smoi, i)
  integer :: i ! soil type
  real :: khyd, smoi
  khyd = Ksat(i)*(smoi/theta_sat(i))**(2.*slbs(i) + 3.)
END FUNCTION khyd

SUBROUTINE init_soil_param(fieldcp, nzg)
  real, parameter :: potwilt = -153. !matric potential at wilting point
  integer :: nsoil, k, irec, nzg
  real, dimension(nzg, nstyp) :: fieldcp

!define theta_cp, the wilting point in terms of matric potential
  do nsoil = 1, nstyp
    slwilt(nsoil) = theta_sat(nsoil)*(slpots(nsoil)/potwilt)**(1./slbs(nsoil))
  end do
end subroutine init_soil_param
