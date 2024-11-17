!     ******************************************************************
subroutine update_wtb_qlat(nzg,slz,dz,wtd,qspring,qlat,smoi,smoieq,soiltextures,smoiwtd,qlatflux,fdepth)
   implicit none
   integer :: nzg,iwtd,kwtd,nsoil,nsoil1,k,k1
   real , dimension(nzg+1) :: slz
   real , dimension(0:nzg+1) :: qlatflux
   real , dimension(nzg) :: dz,vctr4
   real :: wtd,qspring,wtdold,qlat,totwater,smoiwtd,maxwatup,maxwatdw,wgpmid,syielddw,dzup,tempk,fracliq,smoieqwtd,fdepth,smoisat
   real, dimension(nzg) :: smoi,smoieq

   integer, dimension(2) :: soiltextures
   integer, dimension(nzg) :: soiltxt

   do k=1,nzg
      vctr4(k) = 0.5 * (slz(k) + slz(k+1))
   enddo

   where(slz.lt.-0.3)
      soiltxt=soiltextures(1)
   elsewhere
      soiltxt=soiltextures(2)
   endwhere

   qspring=0.
   totwater=qlat

   iwtd=1

!case 1: totwater > 0 (water table going up):
   IF(totwater.gt.0.)then

      do k=2,nzg
         if(wtd.lt.slz(k))exit
      enddo
      iwtd=k
      kwtd=iwtd-1
      nsoil=soiltxt(kwtd)
      smoisat = slmsts(nsoil)*max(min(exp((vctr4(kwtd)+1.5)/fdepth),1.),0.1)
!max water that fits in the layer
      maxwatup=dz(kwtd)*(smoisat-smoi(kwtd))

      if(totwater.le.maxwatup)then
         smoi(kwtd) = smoi(kwtd) + totwater / dz(kwtd)
         qlatflux(kwtd) = qlatflux(kwtd) + totwater
         smoi(kwtd) = min(smoi(kwtd),smoisat)
         if(smoi(kwtd).gt.smoieq(kwtd))wtd = min ( ( smoi(kwtd)*dz(kwtd) &
            - smoieq(kwtd)*slz(iwtd) + smoisat*slz(kwtd) ) / &
            ( smoisat-smoieq(kwtd) ) , slz(iwtd) )
         totwater=0.
      else   !water enough to saturate the layer
         smoi(kwtd) = smoisat
         qlatflux(kwtd) = qlatflux(kwtd) + maxwatup
         totwater=totwater-maxwatup
         k1=iwtd
         do k=k1,nzg+1
            wtd = slz(k)
            iwtd=k+1
            if(k.eq.nzg+1)exit
            nsoil=soiltxt(k)
            smoisat =slmsts(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth),1.),0.1)
            maxwatup=dz(k)*(smoisat-smoi(k))
            if(totwater.le.maxwatup)then
               smoi(k) = smoi(k) + totwater / dz(k)
               qlatflux(k) = qlatflux(k) + totwater
               smoi(k) = min(smoi(k),smoisat)
               if(smoi(k).gt.smoieq(k))wtd = min ( ( smoi(k)*dz(k) &
                  - smoieq(k)*slz(iwtd) + smoisat*slz(k) ) / &
                  ( smoisat-smoieq(k) ) , slz(iwtd) )
               totwater=0.
               exit
            else
               smoi(k) = smoisat
               qlatflux(k) = qlatflux(k) + maxwatup
               totwater=totwater-maxwatup
            endif
         enddo
      endif

!water springing at the surface
      qspring=totwater

!case 2: totwater < 0 (water table going down):
   ELSEIF(totwater.lt.0.)then

      do k=2,nzg
         if(wtd.lt.slz(k))exit
      enddo
      iwtd=k

      k1=iwtd-1
      do kwtd=k1,1,-1

         nsoil=soiltxt(kwtd)
         smoisat = slmsts(nsoil)*max(min(exp((vctr4(kwtd)+1.5)/fdepth),1.),0.1)

!max water that the layer can yield
         maxwatdw=dz(kwtd)*(smoi(kwtd)-smoieq(kwtd))

         if(-totwater.le.maxwatdw)then
            smoi(kwtd) = smoi(kwtd) + totwater / dz(kwtd)
            qlatflux(kwtd) = qlatflux(kwtd) + totwater
            if(smoi(kwtd).gt.smoieq(kwtd))then
               wtd = ( smoi(kwtd)*dz(kwtd) &
                  - smoieq(kwtd)*slz(iwtd) + smoisat*slz(kwtd) ) / &
                  ( smoisat-smoieq(kwtd) )
            else
               wtd=slz(kwtd)
               iwtd=iwtd-1
            endif
            totwater=0.
            exit
         else
            wtd = slz(kwtd)
            iwtd=iwtd-1
            if(maxwatdw.ge.0.)then
               smoi(kwtd) = smoieq(kwtd)
               qlatflux(kwtd) = qlatflux(kwtd) + maxwatdw
               totwater = totwater + maxwatdw
            endif
         endif

      enddo

      if(iwtd.eq.1.and.totwater.lt.0.)then
         nsoil=soiltxt(1)
         smoisat = slmsts(nsoil)*max(min(exp((vctr4(1)+1.5)/fdepth),1.),0.1)

         smoi(1) = smoi(1) + totwater / dz(1)
         qlatflux(1) = qlatflux(1) + totwater
         wtd = max( ( smoi(1)*dz(1) &
            - smoieq(1)*slz(2) + smoisat*slz(1) ) / &
            ( smoisat-smoieq(1) ) , slz(1) )
      endif
      qspring=0.
   ENDIF
end subroutine update_wtb_qlat
