subroutine GW2RIVER(imax,js,je,nzg,slz,deltat,soiltxt,landmask,wtd,maxdepth,riverdepth,width,length,area,fdepth,qrf)
   implicit none
   integer :: i,j,imax,js,je,nsoil,k,iwtd,nzg
   real, dimension(nzg+1) :: slz
   integer ,dimension(2,imax,js:je) :: soiltxt
   integer, dimension(imax,js:je) :: landmask
   real, dimension(imax,js:je) :: wtd,maxdepth,riverdepth,width,length,area,qrf,fdepth
   real :: riversurface,deltat,soilwatercap,rcond,rdepth,hydcon

   soilwatercap=0.
   qrf=0.

   do j=js+1,je-1
      do i=2,imax-1
         if(landmask(i,j).eq.0.or.width(i,j).eq.0.)cycle
         rdepth=max(riverdepth(i,j),0.)

         nsoil=soiltxt(2,i,j)
         riversurface= -( maxdepth(i,j)-rdepth )
         if(riversurface.ge.0.)cycle      !this just in case...

         if(wtd(i,j).gt.riversurface)then
            hydcon = Ksat(nsoil)*max(min(exp((-maxdepth(i,j)+1.5)/fdepth(i,j)),1.),0.1)
            rcond=width(i,j)*length(i,j)*hydcon
            qrf(i,j)=rcond*(wtd(i,j)-riversurface) * ( deltat / area(i,j) )

!if(i.eq.201.and.j.eq.191)write(6,*)'mirar',rivercond(i,j),pexp(i,j),rcond(i,j),wtd(i,j),eqwtd(i,j),qrf(i,j),riversurface,riverbed(i,j)

!limit it to prevent sudden drops , lets say 50mm per day 0.05/86400.
!                  qrf(i,j)=min(qrf(i,j),deltat*0.05/86400.)

         elseif(wtd(i,j).gt.-maxdepth(i,j))then   !water table connected to the river, even though below river surface

            hydcon = Ksat(nsoil)*max(min(exp((-maxdepth(i,j)+1.5)/fdepth(i,j)),1.),0.1)
            rcond=width(i,j)*length(i,j)*hydcon

            soilwatercap=-rcond*(wtd(i,j)-riversurface) * ( deltat / area(i,j) )
!                  soilwatercap=min(soilwatercap,deltat*0.05/86400.)
            qrf(i,j)=-max(min(soilwatercap,riverdepth(i,j)),0.)*min(width(i,j)*length(i,j)/area(i,j),1.)

         else
!water table below river bed, disconnected from the river. No rcond use, just
!infiltration. Assume it occurs at the Ksat rate and water goes directly to the
!water table.
            qrf(i,j) = -max(min(Ksat(nsoil)*deltat,rdepth),0.)  * min(width(i,j)*length(i,j)/area(i,j),1.)
         endif

!if(i.eq.154.and.j.eq.443)write(6,*)'mirar qrf 1',wtd(i,j),maxdepth(i,j),riverdepth(i,j),width(i,j),length(i,j),nsoil,rcond,soilwatercap,qrf(i,j)
!if(i.eq.886.and.j.eq.564)write(6,*)'mirar qrf 2',wtd(i,j),maxdepth(i,j),riverdepth(i,j),width(i,j),length(i,j),nsoil,rcond,soilwatercap,qrf(i,j)
      enddo
   enddo
end subroutine gw2river
