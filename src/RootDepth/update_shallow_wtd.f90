SUBROUTINE update_shallow_wtd(i,j,nzg,freedrain,slz,dz,soiltxt,smoieq,smoiwtd,smoi,wtd,rech,fdepth)

   integer :: nzg,freedrain,nsoil,nsoil1,k,iwtd,kwtd,i,j,flag
   real, dimension(nzg+1) :: slz
   real, dimension(nzg) :: dz,vctr4
   real, dimension(nzg) :: smoieq,smoi
   integer, dimension(2) :: soiltxt
   real :: wtd,wtdold,wgpmid,rech,smoiwtd,dzup,smoieqwtd,fdepth,smoisat

   rech=0.
   flag=0

   do k=1,nzg
      vctr4(k) = 0.5 * (slz(k) + slz(k+1))
   enddo

   do k=1,nzg
!     if(wtd+1.e-6.lt.slz(k))exit
      if(wtd.lt.slz(k))exit
   enddo
   iwtd=k

!if(i.eq.25.and.j.eq.30)write(6,*)'mirar',iwtd,wtd,slz(iwtd),(wtd-slz(iwtd))*1.e3
   kwtd=iwtd-1
   if(kwtd.gt.0)then    !wtd in the resolved layers
      wtdold=wtd

      if(slz(kwtd).lt.-0.30)then
         nsoil=soiltxt(1)
      else
         nsoil=soiltxt(2)
      endif

      if(kwtd.gt.1)then
         smoisat = slmsts(nsoil)*max(min(exp((vctr4(kwtd-1)+1.5)/fdepth),1.),0.1)
         if(wtd.lt.slz(kwtd)+0.01.and.smoi(kwtd-1).lt.smoisat)flag=1
      endif

      smoisat = slmsts(nsoil)*max(min(exp((vctr4(kwtd)+1.5)/fdepth),1.),0.1)

      if(smoi(kwtd).gt.smoieq(kwtd).and.flag.eq.0)then

         if(smoi(kwtd).eq.smoisat)then !wtd went to the layer above
            wtd=slz(iwtd)
            rech=(wtdold-wtd) * (smoisat-smoieq(kwtd))
            iwtd=iwtd+1
            kwtd=kwtd+1
!if(i.eq.300.and.j.eq.200)write(6,*)'mirar 1',rech,wtdold,wtd,kwtd,iwtd,smoi(kwtd),slmsts(nsoil),smoieq(kwtd)
            if(kwtd.le.nzg)then
               if(smoi(kwtd).gt.smoieq(kwtd))then
                  wtdold=wtd

                  if(slz(kwtd).lt.-0.30)then
                     nsoil=soiltxt(1)
                  else
                     nsoil=soiltxt(2)
                  endif

                  smoisat = slmsts(nsoil)*max(min(exp((vctr4(kwtd)+1.5)/fdepth),1.),0.1)
                  wtd = min( ( smoi(kwtd)*dz(kwtd) &
                     - smoieq(kwtd)*slz(iwtd) + smoisat*slz(kwtd) ) / &
                     ( smoisat-smoieq(kwtd) ), slz(iwtd))
                  rech=rech+(wtdold-wtd) * (smoisat-smoieq(kwtd))
!if(i.eq.300.and.j.eq.200)write(6,*)'mirar 2',rech,wtdold,wtd,kwtd,iwtd,smoi(kwtd),slmsts(nsoil),smoieq(kwtd)
               endif
            endif
         else  !wtd stays in the layer
            wtd = min( ( smoi(kwtd)*dz(kwtd) &
               - smoieq(kwtd)*slz(iwtd) + smoisat*slz(kwtd) ) / &
               ( smoisat-smoieq(kwtd) ), slz(iwtd))
            rech=(wtdold-wtd) * (smoisat-smoieq(kwtd))
!if(i.eq.300.and.j.eq.200)write(6,*)'mirar 3',rech,wtdold,wtd,kwtd,iwtd,smoi(kwtd),slmsts(nsoil),smoieq(kwtd),smoi(iwtd),smoieq(iwtd)

         endif

      else    !wtd has gone down to the layer below
         wtd=slz(kwtd)
         rech=(wtdold-wtd) * (smoisat-smoieq(kwtd))
!if(i.eq.300.and.j.eq.200)write(6,*)'mirar 4',rech,wtdold,wtd,kwtd,iwtd,smoi(kwtd),slmsts(nsoil),smoieq(kwtd)
         kwtd=kwtd-1
         iwtd=iwtd-1
!wtd crossed to the layer below. Now adjust it there
         if(kwtd.ge.1)then
            wtdold=wtd

            if(slz(kwtd).lt.-0.30)then
               nsoil=soiltxt(1)
            else
               nsoil=soiltxt(2)
            endif

            smoisat = slmsts(nsoil)*max(min(exp((vctr4(kwtd)+1.5)/fdepth),1.),0.1)

            if(smoi(kwtd).gt.smoieq(kwtd))then
               wtd = min( ( smoi(kwtd)*dz(kwtd) &
                  - smoieq(kwtd)*slz(iwtd) + smoisat*slz(kwtd) ) / &
                  ( smoisat-smoieq(kwtd) ) , slz(iwtd) )
            else
               wtd=slz(kwtd)
            endif
            rech = rech + (wtdold-wtd) * &
               (smoisat-smoieq(kwtd))
!if(i.eq.300.and.j.eq.200)write(6,*)'mirar 5',rech,wtdold,wtd,kwtd,iwtd,smoi(kwtd),slmsts(nsoil),smoieq(kwtd)
         endif

      endif

   endif

   if(wtd.lt.slz(1))write(6,*)'problem with wtd',wtd,i,j

end subroutine update_shallow_wtd
