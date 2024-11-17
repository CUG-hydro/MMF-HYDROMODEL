subroutine update_deep_wtb(imax,jmax,js,je,nzg,slz,dz,soiltxt,wtd,bottomflux,rech &
   ,qslat,qlat,landmask,deltat,smoi,smoieq,smoiwtd,qsprings)

   integer :: imax,jmax,js,je,nzg,i,j,nsoil
   real :: deltat,totwater,qspring,wgpmid,kfup,vt3dbdw,newwgp
   real , dimension(nzg+1) :: slz
   real , dimension(nzg) :: dz
   real,dimension(imax,js:je):: wtd,rech,bottomflux,qslat,qlat &
      ,smoiwtd,qsprings,deeprech
   real,dimension(nzg,imax,js:je):: smoi,smoieq
   integer, dimension(2,imax,js:je)::soiltxt
   integer, dimension(imax,js:je) :: landmask

!calculate deep recharge
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

   DO j=js+1,je-1
      DO i=1,imax
         if(landmask(i,j).eq.1)then
            if(i.eq.300.and.j.eq.300)write(6,*)'mirar qlat',qlat(i,j),qslat(i,j),wtd(i,j)
!Total groundwater balance in the cell
            totwater = qlat(i,j) -qslat(i,j) - deeprech(i,j)

            call update_wtd(nzg,slz,dz,wtd(i,j),qspring,totwater,smoi(1,i,j) &
               ,smoieq(1,i,j),soiltxt(1,i,j),smoiwtd(i,j))

            qsprings(i,j) = qsprings(i,j) + qspring*1.e3
         endif
      ENDDO
   ENDDO
!qlat=qlat*1.e3
end subroutine update_deep_wtb
