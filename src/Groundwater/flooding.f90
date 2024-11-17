subroutine flooding(imax,js,je,deltat,fd,bfd,topo,area,riverwidth,riverlength,riverdepth,floodheight,delsfcwat)
!use module_parallel
   implicit none

   integer, parameter :: ntsplit=1
   integer :: imax,js,je,i,j,i1,j1,i2,j2,ilow,jlow,ii,jj,k,ksn,n
   integer, dimension(imax,js:je) :: fd,bfd
   real, dimension(imax,js:je) :: topo,area,floodheight,dflood,dflood2,delsfcwat,riverwidth,riverlength,riverdepth
   real :: deltat,dh,dhmax,dij,dtotal
   integer :: reqsu,reqsd,reqru,reqrd,reqsu2,reqsd2,reqru2,reqrd2


   dflood=0.
   dflood2=0.

   DO n=1,ntsplit
!communicate flood water height to neighboring cells
      if(numtasks .gt. 1)then
         call sendborders(imax,js,je,floodheight,reqsu,reqsd,reqru,reqrd)
      endif

!make sure that the borders are received before calculating anything
      if(pid.eq.1)then
         call  MPI_wait(reqru,status,ierr)
      elseif(pid.eq.numtasks-2)then
         call  MPI_wait(reqrd,status,ierr)
      elseif(pid.gt.1.and.pid.lt.numtasks-2)then
         call  MPI_wait(reqru,status,ierr)
         call  MPI_wait(reqrd,status,ierr)
      endif

      do j=js+1,je-1
         do i=2,imax-1
!if(i.eq.3417.and.j.eq.6320)write(6,*)'mirar en flood',delsfcwat(i,j),floodheight(i,j)

            if(fd(i,j).eq.0)cycle
            if(floodheight(i,j).gt.0.05)then
!                   call flowdir(imax,js,je,fd,i,j,i1,j1)
!                   call flowdir(imax,js,je,bfd,i,j,i2,j2)

!find the lowest elevation neighbour that is not along the main river channel
               dhmax=0.
               dh=0.
               ilow=i
               jlow=j
               do jj=j-1,j+1
                  do ii=i-1,i+1
!                   if(ii.eq.i1.and.jj.eq.j1)cycle
!                   if(ii.eq.i2.and.jj.eq.j2)cycle
                     if(ii.eq.i.and.jj.eq.j)cycle

                     dh=floodheight(i,j)+topo(i,j)-(floodheight(ii,jj)+topo(ii,jj))
                     if(ii.ne.i.and.jj.ne.j)dh=dh/sqrt(2.)
                     if(dh.gt.dhmax)then
                        ilow=ii
                        jlow=jj
                        dhmax=dh
                     endif
                  enddo
               enddo
!now flood the lowest elevation neighbour

               if(dhmax.gt.0.)then
                  call flowdir(imax,js,je,fd,i,j,i1,j1)

                  dtotal=floodheight(i,j)+floodheight(ilow,jlow)
!         dij=max(0.5*(topo(ilow,jlow)-topo(i,j)+dtotal),0.)
!         dflood(i,j)=dflood(i,j)+dij-floodheight(i,j)
!         dflood(ilow,jlow)=dflood(ilow,jlow)+(dtotal-dij)-floodheight(ilow,jlow)
                  dij = max( floodheight(i,j)-max(0.5*(topo(ilow,jlow)-topo(i,j)+dtotal),0.) , 0.)
                  if(ilow.eq.i1.and.jlow.eq.j1) &
                     dij=max(dij-(riverwidth(i,j)*floodheight(i,j)*riverlength(i,j))/area(i,j),0.) !the flow along the river channel is taken care of by the river routine
                  if(delsfcwat(i,j).lt.0.)dij=max(min(dij,floodheight(i,j)+delsfcwat(i,j)),0.)
                  dflood(i,j) = dflood(i,j) - dij
                  dflood(ilow,jlow) = dflood(ilow,jlow) + dij*area(i,j)/area(ilow,jlow)
               endif

!if(ilow.eq.540.and.jlow.eq.641)write(6,*)'mirar flood',floodheight(i,j),dhmax,dij,dflood(i,j),i,j,floodheight(ilow,jlow),dflood(ilow,jlow),floodheight(i,j)+topo(i,j),floodheight(ilow,jlow)+topo(ilow,jlow)
            endif
         enddo
      enddo

      if(numtasks.gt.1)then
         dflood2=dflood
         call sendbordersflood(imax,js,je,dflood2,reqsu2,reqsd2,reqru2,reqrd2)
      endif
!make sure that the borders are received before calculating anything
      if(pid.eq.1)then
         call  MPI_wait(reqru2,status,ierr)
      elseif(pid.eq.numtasks-2)then
         call  MPI_wait(reqrd2,status,ierr)
      elseif(pid.gt.1.and.pid.lt.numtasks-2)then
         call  MPI_wait(reqru2,status,ierr)
         call  MPI_wait(reqrd2,status,ierr)
      endif

      if(pid.eq.1)then
         dflood(1:imax,je-1)  = dflood(1:imax,je-1) +  dflood2(1:imax,je-1)
      elseif(pid.eq.numtasks-2)then
         dflood(1:imax,js+1)  = dflood(1:imax,js+1) +  dflood2(1:imax,js+1)
      elseif(pid.gt.1.and.pid.lt.numtasks-2)then
         dflood(1:imax,js+1)  = dflood(1:imax,js+1) +  dflood2(1:imax,js+1)
         dflood(1:imax,je-1)  = dflood(1:imax,je-1) +  dflood2(1:imax,je-1)
      endif

!before changing floodheight make sure that the borders have been received
      if(pid.eq.1)then
         call  MPI_wait(reqsu,status,ierr)
      elseif(pid.eq.numtasks-2)then
         call  MPI_wait(reqsd,status,ierr)
      elseif(pid.gt.1.and.pid.lt.numtasks-2)then
         call  MPI_wait(reqsu,status,ierr)
         call  MPI_wait(reqsd,status,ierr)
      endif

!update floodheight and riverdepth
      delsfcwat=delsfcwat+dflood
!     floodheight=floodheight+delsfcwat
!     riverdepth=riverdepth+delsfcwat
!     delsfcwat=0.

!before changing dflood make sure that the borders have been received
      if(pid.eq.1)then
         call  MPI_wait(reqsu2,status,ierr)
      elseif(pid.eq.numtasks-2)then
         call  MPI_wait(reqsd2,status,ierr)
      elseif(pid.gt.1.and.pid.lt.numtasks-2)then
         call  MPI_wait(reqsu2,status,ierr)
         call  MPI_wait(reqsd2,status,ierr)
      endif
   ENDDO
end subroutine flooding
