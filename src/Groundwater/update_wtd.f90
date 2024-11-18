subroutine update_wtd(nzg,slz,dz,wtd,qspring,totwater,smoi,smoieq,soiltextures,smoiwtd)
   implicit none
   integer :: nzg,iwtd,kwtd,nsoil,nsoil1,k,k1
   real , dimension(nzg+1) :: slz
   real , dimension(nzg) :: dz
   real :: wtd,qspring,wtdold,totwater,smoiwtd,maxwatup,maxwatdw,wgpmid,syielddw,dzup,tempk,fracliq,smoieqwtd
   real, dimension(nzg) :: smoi,smoieq

   integer, dimension(2) :: soiltextures
   integer, dimension(nzg) :: soiltxt

   where(slz.lt.-0.3)
      soiltxt=soiltextures(1)
   elsewhere
      soiltxt=soiltextures(2)
   endwhere

   iwtd=1
!case 1: totwater > 0 (water table going up):
   IF(totwater.gt.0.)then


      if(wtd.ge.slz(1))then

         do k=2,nzg
            if(wtd.lt.slz(k))exit
         enddo
         iwtd=k
         kwtd=iwtd-1
         nsoil=soiltxt(kwtd)

!max water that fits in the layer
         maxwatup=dz(kwtd)*(theta_sat(nsoil)-smoi(kwtd))

         if(totwater.le.maxwatup)then
            smoi(kwtd) = smoi(kwtd) + totwater / dz(kwtd)
            smoi(kwtd) = min(smoi(kwtd),theta_sat(nsoil))
            if(smoi(kwtd).gt.smoieq(kwtd))wtd = min ( ( smoi(kwtd)*dz(kwtd) &
               - smoieq(kwtd)*slz(iwtd) + theta_sat(nsoil)*slz(kwtd) ) / &
               ( theta_sat(nsoil)-smoieq(kwtd) ) , slz(iwtd) )
            totwater=0.
         else   !water enough to saturate the layer
            smoi(kwtd) = theta_sat(nsoil)
            totwater=totwater-maxwatup
            k1=iwtd
            do k=k1,nzg+1
               wtd = slz(k)
               iwtd=k+1
               if(k.eq.nzg+1)exit
               nsoil=soiltxt(k)
               maxwatup=dz(k)*(theta_sat(nsoil)-smoi(k))
               if(totwater.le.maxwatup)then
                  smoi(k) = smoi(k) + totwater / dz(k)
                  smoi(k) = min(smoi(k),theta_sat(nsoil))
                  if(smoi(k).gt.smoieq(k))wtd = min ( ( smoi(k)*dz(k) &
                     - smoieq(k)*slz(iwtd) + theta_sat(nsoil)*slz(k) ) / &
                     ( theta_sat(nsoil)-smoieq(k) ) , slz(iwtd) )
                  totwater=0.
                  exit
               else
                  smoi(k) = theta_sat(nsoil)
                  totwater=totwater-maxwatup
               endif

            enddo

         endif

      elseif(wtd.ge.slz(1)-dz(1))then ! wtd below bottom of soil model

         nsoil=soiltxt(1)
         maxwatup=(theta_sat(nsoil)-smoiwtd)*dz(1)

         if(totwater.le.maxwatup)then
            smoieqwtd = theta_sat(nsoil) * ( slpots(nsoil) / &
               (slpots(nsoil) - dz(1)) ) ** (1./slbs(nsoil))
            smoieqwtd = max(smoieqwtd,theta_cp(nsoil))

            smoiwtd = smoiwtd + totwater / dz(1)
            smoiwtd = min(smoiwtd,theta_sat(nsoil))
            if(smoiwtd.gt.smoieqwtd)wtd = min( ( smoiwtd*dz(1) &
               - smoieqwtd*slz(1) + theta_sat(nsoil)*(slz(1)-dz(1)) ) / &
               ( theta_sat(nsoil)-smoieqwtd ) , slz(1) )
            totwater=0.
         else
            smoiwtd=theta_sat(nsoil)
            totwater=totwater-maxwatup
            do k=1,nzg+1
               wtd=slz(k)
               iwtd=k+1
               if(k.eq.nzg+1)exit
               nsoil=soiltxt(k)
               maxwatup=dz(k)*(theta_sat(nsoil)-smoi(k))
               if(totwater.le.maxwatup)then
                  smoi(k) = min(smoi(k) + totwater / dz(k),theta_sat(nsoil))
                  if(smoi(k).gt.smoieq(k))wtd = min ( ( smoi(k)*dz(k) &
                     - smoieq(k)*slz(iwtd) + theta_sat(nsoil)*slz(k) ) / &
                     ( theta_sat(nsoil)-smoieq(k) ) , slz(iwtd) )
                  totwater=0.
                  exit
               else
                  smoi(k) = theta_sat(nsoil)
                  totwater=totwater-maxwatup
               endif
            enddo
         endif

!deep water table
      else
         nsoil=soiltxt(1)
         maxwatup=(theta_sat(nsoil)-smoiwtd)*(slz(1)-dz(1)-wtd)
         if(totwater.le.maxwatup)then
            wtd = wtd + totwater/(theta_sat(nsoil)-smoiwtd)
            totwater=0.
         else
            totwater=totwater-maxwatup
            wtd=slz(1)-dz(1)
            maxwatup=(theta_sat(nsoil)-smoiwtd)*dz(1)
            if(totwater.le.maxwatup)then
               smoieqwtd = theta_sat(nsoil) * ( slpots(nsoil) / &
                  (slpots(nsoil) - dz(1)) ) ** (1./slbs(nsoil))
               smoieqwtd = max(smoieqwtd,theta_cp(nsoil))

               smoiwtd = smoiwtd + totwater / dz(1)
               smoiwtd = min(smoiwtd,theta_sat(nsoil))
               wtd = ( smoiwtd*dz(1) &
                  - smoieqwtd*slz(1) + theta_sat(nsoil)*(slz(1)-dz(1)) ) / &
                  ( theta_sat(nsoil)-smoieqwtd )
               totwater=0.
            else
               smoiwtd=theta_sat(nsoil)
               totwater=totwater-maxwatup
               do k=1,nzg+1
                  wtd=slz(k)
                  iwtd=k+1
                  if(k.eq.nzg+1)exit
                  nsoil=soiltxt(k)
                  maxwatup=dz(k)*(theta_sat(nsoil)-smoi(k))

                  if(totwater.le.maxwatup)then
                     smoi(k) = smoi(k) + totwater / dz(k)
                     smoi(k) = min(smoi(k),theta_sat(nsoil))
                     if(smoi(k).gt.smoieq(k))wtd = ( smoi(k)*dz(k) &
                        - smoieq(k)*slz(iwtd) + theta_sat(nsoil)*slz(k) ) / &
                        ( theta_sat(nsoil)-smoieq(k) )
                     totwater=0.
                     exit
                  else
                     smoi(k) = theta_sat(nsoil)
                     totwater=totwater-maxwatup
                  endif
               enddo
            endif
         endif
      endif

!water springing at the surface
      qspring=totwater

!case 2: totwater < 0 (water table going down):
   ELSEIF(totwater.lt.0.)then

      if(wtd.ge.slz(1))then !wtd in the resolved layers

         do k=2,nzg
            if(wtd.lt.slz(k))exit
         enddo
         iwtd=k

         k1=iwtd-1
         do kwtd=k1,1,-1

            nsoil=soiltxt(kwtd)

!max water that the layer can yield
            maxwatdw=dz(kwtd)*(smoi(kwtd)-smoieq(kwtd))

            if(-totwater.le.maxwatdw)then
               smoi(kwtd) = smoi(kwtd) + totwater / dz(kwtd)
               if(smoi(kwtd).gt.smoieq(kwtd))then
                  wtd = ( smoi(kwtd)*dz(kwtd) &
                     - smoieq(kwtd)*slz(iwtd) + theta_sat(nsoil)*slz(kwtd) ) / &
                     ( theta_sat(nsoil)-smoieq(kwtd) )
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
                  totwater = totwater + maxwatdw
               endif
            endif

         enddo

         if(iwtd.eq.1.and.totwater.lt.0.)then
            nsoil=soiltxt(1)
            smoieqwtd = theta_sat(nsoil) * ( slpots(nsoil) / &
               (slpots(nsoil) - dz(1)) ) ** (1./slbs(nsoil))
            smoieqwtd = max(smoieqwtd,theta_cp(nsoil))

            maxwatdw=dz(1)*(smoiwtd-smoieqwtd)

            if(-totwater.le.maxwatdw)then

               smoiwtd = smoiwtd + totwater / dz(1)
               wtd = max( ( smoiwtd*dz(1) &
                  - smoieqwtd*slz(1) + theta_sat(nsoil)*(slz(1)-dz(1)) ) / &
                  ( theta_sat(nsoil)-smoieqwtd ) , slz(1)-dz(1) )

            else
               wtd=slz(1)-dz(1)
               smoiwtd = smoiwtd + totwater / dz(1)
!and now even further down
               dzup=(smoieqwtd-smoiwtd)*dz(1)/(theta_sat(nsoil)-smoieqwtd)
               wtd=wtd-dzup
               smoiwtd=smoieqwtd
            endif
         endif

      elseif(wtd.ge.slz(1)-dz(1))then

!if wtd was already below the bottom of the resolved soil crust
         nsoil=soiltxt(1)
         smoieqwtd = theta_sat(nsoil) * ( slpots(nsoil) / &
            (slpots(nsoil) - dz(1)) ) ** (1./slbs(nsoil))
         smoieqwtd = max(smoieqwtd, theta_cp(nsoil))

         maxwatdw=dz(1)*(smoiwtd-smoieqwtd)

         if(-totwater.le.maxwatdw)then

            smoiwtd = smoiwtd + totwater / dz(1)
            wtd = max( ( smoiwtd*dz(1) &
               - smoieqwtd*slz(1) + theta_sat(nsoil)*(slz(1)-dz(1)) ) / &
               ( theta_sat(nsoil)-smoieqwtd ) , slz(1)-dz(1) )
         else
            wtd=slz(1)-dz(1)
            smoiwtd = smoiwtd + totwater / dz(1)
!and now even further down
            dzup=(smoieqwtd-smoiwtd)*dz(1)/(theta_sat(nsoil)-smoieqwtd)
            wtd=wtd-dzup
            smoiwtd=smoieqwtd
         endif

      else
!gmmequilibrium soil moisture content
         nsoil=soiltxt(1)
         wgpmid = theta_sat(nsoil) * ( slpots(nsoil) / &
            (slpots(nsoil) - (slz(1)-wtd)) ) ** (1./slbs(nsoil))
         wgpmid=max(wgpmid,theta_cp(nsoil))
         syielddw=theta_sat(nsoil)-wgpmid
         wtdold=wtd
         wtd = wtdold + totwater/syielddw
!update wtdwgp
         smoiwtd = (smoiwtd*(slz(1)-wtdold)+wgpmid*(wtdold-wtd) ) / (slz(1)-wtd)
      endif
      qspring=0.
   ENDIF

end subroutine update_wtd
