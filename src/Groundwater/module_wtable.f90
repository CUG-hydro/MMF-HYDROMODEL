MODULE module_wtable

   use module_parallel
   use module_rootdepth

   implicit none
   real, parameter :: pi4=3.1415927 * 4.

CONTAINS
   include "flowdir.f90"
   include "flooding.f90"
   include "GW2river.f90"
   include "lateralflow.f90"
   include "moveqrf.f90"
   include "rivers_dw_flood.f90"
   include "rivers_kw_flood.f90"
   include "update_deep_wtb.f90"
   include "update_wtd.f90"


   subroutine WTABLE(imax,jmax,js,je,nzg,slz,dz,area,soiltxt,wtd,bottomflux,rech,qslat,fdepth,topo,landmask,deltat &
      ,smoi,smoieq,smoiwtd,qsprings)

      integer :: imax,jmax,js,je,nzg,i,j,nsoil
      real :: deltat,totwater,qspring,wgpmid,kfup,vt3dbdw,newwgp
      real , dimension(nzg+1) :: slz
      real , dimension(nzg) :: dz
      real,dimension(imax,js:je):: area,fdepth,wtd,rech,bottomflux,qslat,topo &
         ,smoiwtd,klat,qsprings,qlat,deeprech
      real,dimension(nzg,imax,js:je):: smoi,smoieq
      integer, dimension(2,imax,js:je)::soiltxt
      integer, dimension(imax,js:je) :: landmask
      integer :: reqsu,reqsd,reqru,reqrd

      if(numtasks.gt.1)call SENDBORDERS(imax,js,je,wtd,reqsu,reqsd,reqru,reqrd)

!Calculate lateral flow
      qlat=0.
      do j=js,je
         do i=1,imax
            nsoil=soiltxt(1,i,j)
            klat(i,j)=slcons(nsoil)*klatfactor(nsoil)
         enddo
      enddo

!make sure that the borders are received before calculating lateral flow
      if(pid.eq.1)then
         call  MPI_wait(reqru,status,ierr)
      elseif(pid.eq.numtasks-2)then
         call  MPI_wait(reqrd,status,ierr)
      elseif(pid.gt.1.and.pid.lt.numtasks-2)then
         call  MPI_wait(reqru,status,ierr)
         call  MPI_wait(reqrd,status,ierr)
      endif

      call lateralflow(imax,jmax,js,je,wtd,qlat,fdepth,topo,landmask,deltat,area,klat)

      qslat=qslat+qlat*1.e3

!now calculate deep recharge
      deeprech=0.

      DO j=js+1,je-1
         DO i=1,imax
            if(landmask(i,j).eq.1)then

               if(wtd(i,j).lt.slz(1)-dz(1))then

!calculate k for drainage
                  nsoil=soiltxt(1,i,j)
                  wgpmid = 0.5 * (smoiwtd(i,j) + slmsts(nsoil))
                  kfup =  slcons(nsoil) &
                     * (wgpmid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)

!now calculate moisture potential
                  vt3dbdw = slpots(nsoil)  &
                     * (slmsts(nsoil) / smoiwtd(i,j)) ** slbs(nsoil)

!and now flux (=recharge)
                  deeprech(i,j) = deltat * kfup &
                     * ( (slpots(nsoil)-vt3dbdw)/(slz(1)-wtd(i,j))  - 1. )
!now update smoiwtd
                  newwgp=smoiwtd(i,j) + (deeprech(i,j) - bottomflux(i,j)) / (slz(1)-wtd(i,j))
                  if(newwgp.lt.soilcp(nsoil))then
                     deeprech(i,j)=deeprech(i,j)+(soilcp(nsoil)-newwgp)*(slz(1)-wtd(i,j))
                     newwgp=soilcp(nsoil)
                  endif
                  if(newwgp.gt.slmsts(nsoil))then
                     deeprech(i,j)=deeprech(i,j)-(slmsts(nsoil)-newwgp)*(slz(1)-wtd(i,j))
                     newwgp=slmsts(nsoil)
                  endif

                  smoiwtd(i,j)=newwgp
                  rech(i,j) = rech(i,j)  + deeprech(i,j)*1.e3
               endif
            endif
         ENDDO
      ENDDO

      bottomflux=0.
!before changing wtd make sure that the borders have been received

      if(pid.eq.1)then
         call  MPI_wait(reqsu,status,ierr)
      elseif(pid.eq.numtasks-2)then
         call  MPI_wait(reqsd,status,ierr)
      elseif(pid.gt.1.and.pid.lt.numtasks-2)then
         call  MPI_wait(reqsu,status,ierr)
         call  MPI_wait(reqsd,status,ierr)
      endif

!Now update water table and soil moisture
!write(6,*)'now to updatewtd'
      DO j=js+1,je-1
         DO i=1,imax
            if(landmask(i,j).eq.1)then

!Total groundwater balance in the cell
               totwater = qlat(i,j)  - deeprech(i,j)
               if(qlat(i,j).ne.qlat(i,j))write(6,*)'gran problema!',wtd(i,j),qlat(i,j),i,j
               if(i.eq.54.and.j.eq.49)write(6,*)'mirar antes updatewtd',wtd(i,j),qlat(i,j),deeprech(i,j),totwater

               call update_wtd(nzg,slz,dz,wtd(i,j),qspring,totwater,smoi(1,i,j) &
                  ,smoieq(1,i,j),soiltxt(1,i,j),smoiwtd(i,j))

               qsprings(i,j) = qsprings(i,j) + qspring*1.e3

               if(i.eq.54.and.j.eq.49)write(6,*)'mirar despues updatewtd',wtd(i,j)

            endif
         ENDDO
      ENDDO
   end subroutine wtable

END MODULE MODULE_WTABLE
