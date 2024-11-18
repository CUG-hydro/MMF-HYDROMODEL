MODULE module_rootdepth


   implicit none

   integer, parameter ::nvtyp=30,nstyp=13
   real, save, dimension(nstyp) :: slmsts,soilcp,slbs,slcons,slpots,slwilt,klatfactor
   data slmsts/0.395, 0.410, 0.435, 0.485, 0.451, 0.420  &
      ,0.477, 0.476, 0.426, 0.492, 0.482, 0.863, 0.476/
!data fieldcp/.135,.150,.195,.255,.240,.255,.322,.325  &
!            ,.310,.370,.367,.535,.325/
   data soilcp/.050  ,.052  ,.092  ,.170  ,.125  ,.148 &
      ,.195  ,.235  ,.202  ,.257  ,.268  ,.195 ,.235/
   data slbs  /4.05, 4.38, 4.9, 5.3, 5.39, 7.12, 7.75, 8.52  &
      ,10.4, 10.4, 11.4, 7.75, 8.52/
   data slcons /.000176   ,.0001563  ,.00003467  &
      ,.0000072  ,.00000695 ,.0000063  &
      ,.0000017  ,.00000245 ,.000002167  &
      ,.000001033,.000001283,.0000080,.000005787/
!2m/day            ,.0000017  ,.000023148 ,.000002167  &
!2m/day            ,.000001033,.000023148,.0000080/
!1m/day            ,.0000017  ,.000011574 ,.000002167  &
!1m/day            ,.000001033,.000011574,.0000080/
!0.5m/day           ,.000001033,.000005787,.0000080/
   data slpots/-0.121, -0.090, -0.218, -0.786, -0.478, -0.299  &
      ,-0.356, -0.630, -0.153, -0.490, -0.405, -0.356, -0.630/
   data klatfactor /2.,3.,4.,10.,12.,14.,20.,24.,28.,40.,48.,48.,48./
!data klatfactor /2.,3.,4.,10.,12.,14.,20.,100.,28.,40.,100.,48./
!data klatfactor /2.,3.,4.,10.,12.,14.,20.,48.,28.,40.,48.,48./

CONTAINS

   include "tridag.f90"
   include "update_shallow_wtd.f90"
   include "update_wtb_qsat.f90"
   include "extraction.f90"
   include "potevap.f90"
   include "soilfluxes.f90"
   include "interception.f90"
   include "init_soil.f90"


   SUBROUTINE ROOTDEPTH(freedrain,imax,js,je,nzg,slz,dz,deltat,landmask,veg,hveg,soiltxt,wind,temp,qair,press,netrad,rshort &
      ,lai,precip,qsrun,smoi,smoieq,smoiwtd,wtd,waterdeficit,watext,watextdeep,rech,deeprech &
      ,et_s,et_i,et_c,intercepstore,ppacum,pppendepth,pppendepthold &
      ,qlat,qlatsum,qsprings,inactivedays,maxinactivedays,fieldcp,fdepth,steps,floodheight &
      ,qrf,delsfcwat,icefactor,wtdflux,et_s_daily,et_c_daily,transptop,infilk)

!real, parameter :: minpprate=1./3. !pp (ammount in mm per timestep) above which there is no intercpetion loss
      real, parameter :: minpprate=0.01 !pp (ammount in mm per timestep) above which there is no intercpetion loss
      integer :: imax,js,je,nzg,i,j,k,freedrain,itime,maxinactivedays,floodflag
      real :: deltat,steps
      real, dimension(nzg,nstyp) :: fieldcp
      real, dimension(nzg+1) :: slz,flux
      real, dimension(0:nzg+1) :: qlatflux
      integer, dimension(0:nzg+1,imax,js:je) :: inactivedays
      integer*1, dimension(imax,js:je,26:40) :: icefactor
      integer*1, dimension(nzg) :: icefac
      integer*1, dimension(imax,js:je) :: infilk,pppendepthold
      real, dimension(nzg) :: dz,dsmoi
      integer, dimension(imax,js:je) :: landmask
      real, dimension(imax,js:je) :: veg,hveg,wind,temp,qair,press,netrad,rshort,lai,precip,qsrun,ppacum,waterdeficit,intercepstore &
         ,et_s,et_i,et_c,watextdeep,pppendepth,qlat,qlatsum,qsprings,floodheight,qrf,delsfcwat
      real, dimension(imax,js:je) ::  wtdflux,et_s_daily,et_c_daily,transptop
      integer, dimension(2,imax,js:je) :: soiltxt
      real, dimension(nzg,imax,js:je) :: smoi,watext,smoieq

      real :: petstep_s,petstep_c,petstep_w,petstep_i,etstep_s,etstep_c,etstep_i,runoff,rechstep,ppdrip,watdef &
         ,dsmoideep,qlatstep,pppendepthstep,qrfstep,qrfcorrect,floodstep,wtdold
      real :: delta,gamma,lambda,ra_a,ra_c,rs_c,R_a,R_s,petfactor_s,petfactor_c
      real, dimension(imax,js:je) :: smoiwtd,wtd,rech,deeprech,fdepth
      integer*1 :: infilkstep

!maxinactivedays = 30*8*nint(steps)   !30 days * 8 3h periods * steps in the rootactivy calculation

!smoiwtd=0.3
!wtd=0.
!deeprech=0.
      icefac=0

      DO j=js+1,je-1
         DO i=1,imax

            if(landmask(i,j).eq.0)cycle

!       ppacum(i,j) = ppacum(i,j) + precip(i,j)

!gmm calculate PET

!if(i.eq.2761.and.j.eq.571)write(6,*)'temp,rad,press,qair,rshort,wind',temp(i,j),netrad(i,j),press(i,j),qair(i,j),rshort(i,j),wind(i,j)

!     call potevap_priestly_taylor(i,j,temp(i,j),rad(i,j),press(i,j),petstep)
!     call potevap_Penman_Monteith(i,j,temp(i,j),netrad(i,j),rshort(i,j),press(i,j),qair(i,j) &
!                                  ,wind(i,j),lai(i,j),veg(i,j),hveg(i,j),petstep)

            icefac(26:40)=icefactor(i,j,26:40)

            if(floodheight(i,j).gt.0.05)then
               floodflag=1
            else
               floodflag=0
            endif

            call potevap_Shutteworth_Wallace(i,j,deltat,temp(i,j),netrad(i,j),rshort(i,j),press(i,j),qair(i,j) &
               ,wind(i,j),lai(i,j),veg(i,j),hveg(i,j) &
               ,delta,gamma,lambda,ra_a,ra_c,rs_c,R_a,R_s &
!                                  ,petstep_s,petstep_c,petstep_w,petstep_i,floodflag)
               ,petfactor_s,petfactor_c,petstep_w,petstep_i,floodflag)

            et_s(i,j) = et_s(i,j) + petstep_w
            if(floodflag.eq.1.and.nint(veg(i,j)).le.1)delsfcwat(i,j) = delsfcwat(i,j) - petstep_w *1.e-3

            if(nint(veg(i,j)).le.1)cycle

!if(i.eq.2761.and.j.eq.571)write(6,*)'pet',petstep,precip(i,j)

!gmm first interception

            call interception(minpprate,precip(i,j),lai(i,j),intercepstore(i,j),ppdrip,petstep_i,etstep_i)

            et_i(i,j) = et_i(i,j) + etstep_i

!gmm then extraction
!gmm now see where pet has to be transpired

            ppdrip=ppdrip/steps !in mm
            floodstep = floodheight(i,j)/steps !in m
!petstep_s=petstep_s/steps !in mm
!petstep_c=petstep_c/steps !in mm
            qlatstep=qlat(i,j)/steps !in m
            qrfstep=qrf(i,j)/steps !in m

!if(i.eq.3417.and.j.eq.6320)write(6,*),'mirar antes',floodheight(i,j),floodstep,delsfcwat(i,j)

            flux=0.
            qlatflux=0.

            wtdold=wtd(i,j)

            do itime=1,nint(steps)

               call extraction(i,j,nzg,slz,dz,deltat/steps,soiltxt(1,i,j),wtd(i,j),smoi(1,i,j),smoiwtd(i,j) &
                  ,delta,gamma,lambda,lai(i,j),ra_a,ra_c,rs_c,R_a,R_s,petfactor_s,petfactor_c,petstep_s &
                  ,petstep_c,watdef,dsmoi,dsmoideep,inactivedays(0,i,j),maxinactivedays,fieldcp,hveg(i,j),fdepth(i,j) &
                  ,icefac)

               et_c(i,j) = et_c(i,j) + petstep_c - watdef*1.e3
               waterdeficit(i,j) = waterdeficit(i,j) + watdef*1.e3
               watext(:,i,j) = watext(:,i,j) + dsmoi(:)*1.e3
               transptop(i,j) = transptop(i,j) + dsmoi(nzg)*1.e3
               et_c_daily(i,j) = et_c_daily(i,j) + petstep_c - watdef*1.e3
!     watextdeep(i,j) = watextdeep(i,j) + dsmoideep*1.e3

!now update soil moisture from transpiration, evaporation, infiltration and soil
!fluxes

!dsmoi=dsmoi/steps
!dsmoideep=dsmoideep/steps
!ppdrip=ppdrip/steps
!petstep_s=petstep_s/steps

!do itime=1,36

!if(i.eq.46.and.j.eq.61)write(6,*)'mirar antes soilfluxes',qlat(i,j),wtd(i,j)
               call soilfluxes(i,j,nzg,freedrain,deltat/steps,slz,dz,soiltxt(1,i,j),smoiwtd(i,j),dsmoi,dsmoideep  &
                  ,smoi(1,i,j),wtd(i,j),rechstep,deeprech(i,j),ppdrip,petstep_s,etstep_s,runoff,flux &
                  ,fdepth(i,j),qlatstep,qlatflux,qrfstep,qrfcorrect,floodstep,icefac)

               delsfcwat(i,j) = delsfcwat(i,j) - max(floodstep-runoff,0.) !in m
               qsrun(i,j) = qsrun(i,j) + max(runoff-floodstep,0.) !in m
!if(i.eq.46.and.j.eq.61)write(6,*)'mirar after soilflux',delsfcwat(i,j),qsrun(i,j),floodstep,runoff,ppdrip
               rech(i,j) = rech(i,j) + rechstep*1.e3
               et_s(i,j) = et_s(i,j) + etstep_s
               et_s_daily(i,j) = et_s_daily(i,j) + etstep_s
               ppacum(i,j) = ppacum(i,j) + ppdrip
!correct qrf. qrfstep should be zero after soilfluxes if there is no problem
               qrf(i,j) = qrf(i,j) + qrfcorrect

!now adjust wtd

!if(i.eq.46.and.j.eq.61)write(6,*)'mirar antes shallowwtd',wtd(i,j)
               call update_shallow_wtd(i,j,nzg,freedrain,slz,dz,soiltxt(1,i,j),smoieq(1,i,j),smoiwtd(i,j),smoi(1,i,j),wtd(i,j),rechstep,fdepth(i,j))
!if(i.eq.46.and.j.eq.61)write(6,*)'mirar despues shallowwtd',wtd(i,j)

               rech(i,j) = rech(i,j) + rechstep*1.e3

!       wtdflux(i,j) = wtdflux(i,j) + (rechstep - qlatstep + qrfstep + qrfcorrect) * 1e3

!if(i.eq.50.and.j.eq.50)write(6,*)'mirar2',wtd(i,j),rechstep*1e3,qlatstep*1e3,qrfstep*1e3,qrfcorrect*1e3

!               call updatewtdqlat(nzg,slz,dz,wtd(i,j),runoff,qlatstep,smoi(1,i,j) &
!                      ,smoieq(1,i,j),soiltxt(1,i,j),smoiwtd(i,j),qlatflux,fdepth(i,j))
!               qsrun(i,j) = qsrun(i,j) + runoff*1.e3
!               qsprings(i,j) = qsprings(i,j) + runoff*1.e3
               qlatsum(i,j)=qlatsum(i,j)+qlatstep
!if(i.eq.300.and.j.eq.300)write(6,*)'mirar qlat antes',qlat(i,j),qlatstep,qlatsum(i,j),wtd(i,j)

!if(i.eq.886.and.j.eq.564)write(6,*)'mirar rootdepth',wtd(i,j),qrf(i,j)/steps,qrfstep,qlatstep,ppdrip*1.e-3,floodstep,rechstep
            enddo
!     do k=1,nzg
!       if(0.5*(flux(k)+flux(k+1)).lt.-1.e-6)then
!            if(pppendepth(i,j).gt.0.5*(slz(k)+slz(k+1)))pppendepth(i,j)=0.5*(slz(k)+slz(k+1))
!            exit
!       endif
!     enddo

            infilkstep=nzg+1
            pppendepthstep=0.
            flux(nzg+1)=-1.
            do k=nzg,0,-1
               if(k.le.nzg-2)then
!           if(pppendepthold(i,j).gt.slz(k+2)+0.5*dz(k+2))exit
                  if(pppendepthold(i,j).ge.k+3)exit
               endif
!       if(flux(k+1).lt.-1.e-5)then
               if(flux(k+1).lt.-0.333e-5)then  !since the timestep was reduced from 3 to 1 h, change the threshold accordingly
                  if(k.eq.0)then
                     if(-flux(1).gt.-qlatflux(k).and.pppendepthstep.gt.slz(1))then
                        pppendepthstep=slz(1)
                        infilkstep=1
                     endif
                  elseif(-flux(k+1)+flux(k).gt.-qlatflux(k)+dsmoi(k).and.pppendepthstep.gt.slz(k+1))then
                     pppendepthstep=slz(k+1)
                     infilkstep=k+1
                  endif
               endif
            enddo

            pppendepthold(i,j) = infilkstep

            if(pppendepth(i,j).gt.pppendepthstep)pppendepth(i,j)=pppendepthstep

            if(slz(max(infilkstep-1,1)).le.wtdold)wtdflux(i,j)=wtdflux(i,j)-flux(infilkstep)*1.e3
            if(infilk(i,j).gt.infilkstep)infilk(i,j)=infilkstep

!if(i.eq.300.and.j.eq.200)write(6,*)'mirar ppdepth',k,flux(k),flux(k+1),pppendepth(i,j),slz(k)

!if(i.eq.2761.and.j.eq.571)write(6,*)'smoi',(smoi(k,i,j),k=1,nzg)
!if(i.eq.2761.and.j.eq.571)write(6,*)'ncount',(ncount(k,i,j),k=1,nzg)
         ENDDO
      ENDDO
   end subroutine rootdepth


!**********************************************************************************************
   FUNCTION khyd(smoi,nsoil)
      integer :: nsoil
      real :: khyd,smoi

      khyd = slcons(nsoil) * (smoi  / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)
   END FUNCTION khyd

END MODULE MODULE_ROOTDEPTH
