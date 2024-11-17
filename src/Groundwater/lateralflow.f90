subroutine LATERAL(imax,jmax,js,je,soiltxt,wtd,qlat,fdepth,topo,landmask,deltat,area,lats,dxy)
   integer :: imax,jmax,js,je,i,j,nsoil
   real :: deltat,dxy
   real,dimension(imax,js:je):: area,fdepth,wtd,qlat,topo,klat,lats
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
   call lateralflow4(imax,jmax,js,je,wtd,qlat,fdepth,topo,landmask,deltat,area,klat,lats,dxy)

!before changing wtd make sure that the borders have been received
   if(pid.eq.1)then
      call  MPI_wait(reqsu,status,ierr)
   elseif(pid.eq.numtasks-2)then
      call  MPI_wait(reqsd,status,ierr)
   elseif(pid.gt.1.and.pid.lt.numtasks-2)then
      call  MPI_wait(reqsu,status,ierr)
      call  MPI_wait(reqsd,status,ierr)
   endif
end subroutine lateral



subroutine LATERALFLOW(imax,jmax,js,je,wtd,qlat,fdepth,topo,landmask,deltat,area,klat)
   implicit none
   real :: deltat,fangle,q
   integer :: imax,jmax,js,je,i,j
   integer, dimension(imax,js:je):: landmask
   real,dimension(imax,js:je)::fdepth,wtd,qlat,topo,area,kcell,klat,head

   fangle=sqrt(tan(pi4/32.))/(2.*sqrt(2.))
!gmmlateral flow calculation
!WHERE(fdepth.lt.1.e-6)
!  kcell=0.
!ELSEWHERE(wtd.lt.-1.5)
!   kcell=fdepth*klat*exp((wtd+1.5)/fdepth)
!ELSEWHERE
!   kcell=klat*(wtd+1.5+fdepth)
!END WHERE
   do j=max(js,1),min(je,jmax)
      do i=1,imax
         if(fdepth(i,j).lt.1.e-6)then
            kcell(i,j)=0.
         elseif(wtd(i,j).lt.-1.5)then
            kcell(i,j)=fdepth(i,j)*klat(i,j)*exp((wtd(i,j)+1.5)/fdepth(i,j))
         else
            kcell(i,j)=klat(i,j)*(wtd(i,j)+1.5+fdepth(i,j))
         endif

         head(i,j) = topo(i,j) + wtd(i,j)
      enddo
   enddo

!head=topo+wtd
   do j=js+1,je-1
      do i=2,imax-1
         IF(landmask(i,j).eq.1) then
            q=0.

            q  = q + (kcell(i-1,j+1)+kcell(i,j)) &
               * (head(i-1,j+1)-head(i,j))/sqrt(2.)

            q  = q +  (kcell(i-1,j)+kcell(i,j)) &
               *  (head(i-1,j)-head(i,j))

            q  = q +  (kcell(i-1,j-1)+kcell(i,j)) &
               * (head(i-1,j-1)-head(i,j))/sqrt(2.)

            q  = q +  (kcell(i,j+1)+kcell(i,j)) &
               * (head(i,j+1)-head(i,j))

            q  = q +  (kcell(i,j-1)+kcell(i,j)) &
               * (head(i,j-1)-head(i,j))

            q  = q +  (kcell(i+1,j+1)+kcell(i,j)) &
               * (head(i+1,j+1)-head(i,j))/sqrt(2.)

            q  = q +  (kcell(i+1,j)+kcell(i,j)) &
               * (head(i+1,j)-head(i,j))

            q  = q +  (kcell(i+1,j-1)+kcell(i,j)) &
               * (head(i+1,j-1)-head(i,j))/sqrt(2.)

            qlat(i,j) = fangle* q * deltat / area(i,j)

         ENDIF
      enddo
   enddo
end subroutine lateralflow



subroutine LATERALFLOW4(imax,jmax,js,je,wtd,qlat,fdepth,topo,landmask,deltat,area,klat,xlat,dxy)
   implicit none
   double precision,parameter :: d2r = 0.0174532925199
   real :: deltat,fangle,q,dxy
   integer :: imax,jmax,js,je,i,j
   integer, dimension(imax,js:je):: landmask
   real,dimension(imax,js:je)::fdepth,wtd,qlat,topo,area,kcell,klat,head,xlat

!gmmlateral flow calculation
!WHERE(fdepth.lt.1.e-6)
!  kcell=0.
!ELSEWHERE(wtd.lt.-1.5)
!   kcell=fdepth*klat*exp((wtd+1.5)/fdepth)
!ELSEWHERE
!   kcell=klat*(wtd+1.5+fdepth)
!END WHERE
   do j=max(js,1),min(je,jmax)
      do i=1,imax
         if(fdepth(i,j).lt.1.e-6)then
            kcell(i,j)=0.
         elseif(wtd(i,j).lt.-1.5)then
            kcell(i,j)=fdepth(i,j)*klat(i,j)*exp((wtd(i,j)+1.5)/fdepth(i,j))
         else
            kcell(i,j)=klat(i,j)*(wtd(i,j)+1.5+fdepth(i,j))
         endif

         head(i,j) = topo(i,j) + wtd(i,j)
      enddo
   enddo

!head=topo+wtd

   do j=js+1,je-1
      do i=2,imax-1
         IF(landmask(i,j).eq.1) then
            q=0.
!north
            q  = q + (kcell(i,j+1)+kcell(i,j)) &
               * (head(i,j+1)-head(i,j)) &
               * cos( d2r * (xlat(i,j) + 0.5*dxy) )
!south
            q  = q +  (kcell(i,j-1)+kcell(i,j)) &
               *  (head(i,j-1)-head(i,j)) &
               * cos( d2r * (xlat(i,j) - 0.5*dxy) )
!west
            q  = q +  (kcell(i-1,j)+kcell(i,j)) &
               * (head(i-1,j)-head(i,j)) &
               / cos( d2r * xlat(i,j) )
!east
            q  = q +  (kcell(i+1,j)+kcell(i,j)) &
               * (head(i+1,j)-head(i,j)) &
               / cos( d2r * xlat(i,j) )

            qlat(i,j) =  0.5 * q * deltat / area(i,j)
         ENDIF
      enddo
   enddo
end subroutine lateralflow4
