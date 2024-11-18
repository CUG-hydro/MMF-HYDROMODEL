subroutine RIVERS_DW_FLOOD(imax, js, je, deltat, dtlr, fd, bfd, qnew, qs, qrf, delsfcwat &
                           , slope, depth, width, length, maxdepth, area, riverarea, floodarea, riverchannel &
                           , qmean, floodheight, topo)

  implicit none
  real, parameter :: gg = 9.81
!integer, parameter :: ntsplit=15
!integer, parameter :: ntsplit=4
  integer :: i, j, imax, js, je, n, i1, j1, i2, j2
  integer, dimension(imax, js:je) :: fd, bfd
  real, dimension(imax, js:je) :: q, qin, qnew, qs, qrf, qext, delsfcwat &
                                  , slope, depth, width, length, maxdepth, area &
                                  , riverarea, qmean, floodheight, riverchannel, floodarea, topo
  real :: deltat, snew, aa, wi, speed, frwtd, dtlr, dsnew, flowwidth, slopeinst &
          , dtopo, dcommon, qmax, vmax, waterelevij, waterelevi1j1, waterelevi2j2, slopefor, slopeback
  integer :: reqsu, reqsd, reqru, reqrd

  do j = js + 1, je - 1
    do i = 2, imax - 1
      IF (fd(i, j) .ne. 0) then

        qext(i, j) = (qrf(i, j) + qs(i, j) + delsfcwat(i, j))/deltat*area(i, j)
!          riverarea(i,j) = width(i,j)*length(i,j)
!          floodarea(i,j) = max( area(i,j)-riverarea(i,j) , 0. )
!          riverchannel(i,j) = maxdepth(i,j)*riverarea(i,j)

   if (i .eq. 72 .and. j .eq. 26) write (6, *) 'mirar rivers 0', qrf(i, j), qs(i, j), delsfcwat(i, j), qnew(i, j), floodheight(i, j)

!if(i.eq.543.and.j.eq.172)write(6,*)'mirar rivers',width(i,j),riverchannel(i,j)
      END IF
    end do
  end do

!dtlr = deltat/float(ntsplit)
!do n=1,ntsplit

  if (numtasks .gt. 1) call sendborders(imax, js, je, qnew, reqsu, reqsd, reqru, reqrd)

!make sure that the borders are received before calculating anything
  if (pid .eq. 1) then
    call MPI_wait(reqru, status, ierr)
  elseif (pid .eq. numtasks - 2) then
    call MPI_wait(reqrd, status, ierr)
  elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
    call MPI_wait(reqru, status, ierr)
    call MPI_wait(reqrd, status, ierr)
  end if

  q = qnew

  qin = 0.

  do j = js, je
    do i = 1, imax
      if (fd(i, j) .gt. 0) then
        call flowdir(imax, js, je, fd, i, j, i1, j1)
        if (i1 .gt. 1 .and. i1 .lt. imax .and. j1 .gt. js .and. j1 .lt. je) then
          qin(i1, j1) = qin(i1, j1) + q(i, j)
        end if
      end if
    end do
  end do

  do j = js + 1, je - 1
    do i = 2, imax - 1

      IF (fd(i, j) .ne. 0) then

!calculate total inflow into cell i j
!fd (flow direction) tells where the river in cell i j is flowing to

        dsnew = qin(i, j) - q(i, j)

!Taquari
        if (i .eq. 4498 .and. j .eq. 4535) dsnew = dsnew + q(4499, 4534)/4.
        if (i .eq. 4498 .and. j .eq. 4534) dsnew = dsnew - q(4499, 4534)/4.
!Taquari
        if (i .eq. 4464 .and. j .eq. 4536) dsnew = dsnew + q(4465, 4535)/2.
        if (i .eq. 4465 .and. j .eq. 4534) dsnew = dsnew - q(4465, 4535)/2.
!Taquari
        if (i .eq. 4346 .and. j .eq. 4560) dsnew = dsnew + q(4346, 4561)/3.
        if (i .eq. 4345 .and. j .eq. 4561) dsnew = dsnew - q(4346, 4561)/3.
!Taquari
        if (i .eq. 4444 .and. j .eq. 4551) dsnew = dsnew + q(4444, 4552)/3.
        if (i .eq. 4443 .and. j .eq. 4553) dsnew = dsnew - q(4444, 4552)/3.
!Taquari
        if (i .eq. 4350 .and. j .eq. 4497) dsnew = dsnew + q(4352, 4496)/3.
        if (i .eq. 4351 .and. j .eq. 4496) dsnew = dsnew - q(4352, 4496)/3.

!Sao Lourenco
        if (i .eq. 4439 .and. j .eq. 4772) dsnew = dsnew + q(4440, 4773)/2.
        if (i .eq. 4440 .and. j .eq. 4772) dsnew = dsnew - q(4440, 4773)/2.

        if (i .eq. 4400 .and. j .eq. 4685) dsnew = dsnew + q(4401, 4685)/2.
        if (i .eq. 4400 .and. j .eq. 4684) dsnew = dsnew - q(4401, 4685)/2.

        if (i .eq. 4418 .and. j .eq. 4688) dsnew = dsnew + q(4418, 4689)/5.
        if (i .eq. 4417 .and. j .eq. 4689) dsnew = dsnew - q(4418, 4689)/5.

        if (i .eq. 4367 .and. j .eq. 4698) dsnew = dsnew + q(4368, 4699)/6.
        if (i .eq. 4368 .and. j .eq. 4698) dsnew = dsnew - q(4368, 4699)/6.

        if (i .eq. 4363 .and. j .eq. 4667) dsnew = dsnew + q(4364, 4668)/6.
        if (i .eq. 4364 .and. j .eq. 4667) dsnew = dsnew - q(4364, 4668)/6.

        if (i .eq. 4475 .and. j .eq. 4718) dsnew = dsnew + q(4475, 4717)/6.
        if (i .eq. 4474 .and. j .eq. 4717) dsnew = dsnew - q(4475, 4717)/6.
!new river store
        snew = depth(i, j)*riverarea(i, j) + floodheight(i, j)*floodarea(i, j) + (dsnew + qext(i, j))*dtlr

!if(i.eq.1360.and.j.eq.1065)write(6,*)'mirar
!rivers',snew,floodheight(i,j),dsnew,qext(i,j),depth(i,j),q(i,j)

!now redistribute water between river channel and floodplain and calculate new
!riverdepth and floodheight
   if (snew .ne. snew) write (6, *) 'problem with snew', i, j, dsnew, qin(i, j), q(i, j), floodheight(i, j), depth(i, j), qext(i, j)
        if (snew .ge. riverchannel(i, j)) then

          floodheight(i, j) = (snew - riverchannel(i, j))/max(area(i, j), riverarea(i, j))
          depth(i, j) = floodheight(i, j) + maxdepth(i, j)

!if(i.eq.1360.and.j.eq.1065)write(6,*)'mirar rivers
!2',floodheight(i,j),depth(i,j)
        else

          floodheight(i, j) = 0.
          if (riverarea(i, j) .gt. 0.) then
            depth(i, j) = snew/riverarea(i, j)
          else
            depth(i, j) = 0.
          end if
!if(i.eq.1360.and.j.eq.1065)write(6,*)'mirar rivers
!2',floodheight(i,j),depth(i,j)

        end if
        if(depth(i,j).ne.depth(i,j))write(6,*)'problem with depth',i,j,qrf(i,j),qs(i,j),delsfcwat(i,j),qnew(i,j),floodheight(i,j)
!if(i.eq.3417.and.j.eq.6320)write(6,*)'mirar rivers 1',qrf(i,j),qs(i,j),delsfcwat(i,j),qnew(i,j),floodheight(i,j),depth(i,j)
!if(i.eq.158.and.j.eq.441)write(6,*)'mirar rivers 2',qrf(i,j),qs(i,j),delsfcwat(i,j),qnew(i,j),floodheight(i,j),depth(i,j)
      END IF
    end do
  end do

!before changing qnew make sure that the borders have been received
  if (pid .eq. 1) then
    call MPI_wait(reqsu, status, ierr)
  elseif (pid .eq. numtasks - 2) then
    call MPI_wait(reqsd, status, ierr)
  elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
    call MPI_wait(reqsu, status, ierr)
    call MPI_wait(reqsd, status, ierr)
  end if

  if (numtasks .gt. 1) call sendborders(imax, js, je, depth, reqsu, reqsd, reqru, reqrd)

!make sure that the borders are received before calculating anything
  if (pid .eq. 1) then
    call MPI_wait(reqru, status, ierr)
  elseif (pid .eq. numtasks - 2) then
    call MPI_wait(reqrd, status, ierr)
  elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
    call MPI_wait(reqru, status, ierr)
    call MPI_wait(reqrd, status, ierr)
  end if

  do j = js + 1, je - 1
    do i = 2, imax - 1

      if (fd(i, j) .gt. 0) then

        call flowdir(imax, js, je, fd, i, j, i1, j1)
        if (floodheight(i, j) .gt. 0.05 .or. depth(i1, j1) .gt. maxdepth(i1, j1) + 0.05) then

          if (width(i, j) .gt. 0) then
            aa = depth(i, j)*width(i, j)/(2.*depth(i, j) + width(i, j))
          else
            aa = depth(i, j)
          end if

          waterelevij = topo(i, j) - maxdepth(i, j) + depth(i, j)
          waterelevi1j1 = topo(i1, j1) - maxdepth(i1, j1) + max(depth(i1, j1), 0.)
          slopefor = (waterelevi1j1 - waterelevij)/(0.5*(length(i, j) + length(i1, j1)))
          if (bfd(i, j) .gt. 0) then
            call flowdir(imax, js, je, bfd, i, j, i2, j2)
            waterelevi2j2 = topo(i2, j2) - maxdepth(i2, j2) + max(depth(i2, j2), 0.)
            slopeback = (waterelevij - waterelevi2j2)/(0.5*(length(i2, j2) + length(i, j)))
            slopeinst = 0.5*(slopefor + slopeback)
          else
            slopeinst = slopefor
          end if

          qnew(i, j) = (q(i, j) - gg*depth(i, j)*dtlr*slopeinst)/ &
                       (1.+gg*dtlr*0.03**2.*q(i, j)/(aa**(4./3.)*depth(i, j)))

          if (width(i, j) .eq. 0.) then
            flowwidth = sqrt(area(i, j))
          else
            flowwidth = width(i, j)
          end if
          qnew(i, j) = qnew(i, j)*flowwidth

        else
          aa = depth(i, j)*width(i, j)/(2.*depth(i, j) + width(i, j))
          speed = (aa**(2./3.))*sqrt(slope(i, j))/0.03
          speed = max(min(speed, length(i, j)/dtlr), 0.01)

!now calculate the new q
          qnew(i, j) = speed*depth(i, j)*width(i, j)
        end if
      else
        qnew(i, j) = 0.
      end if
    end do
  end do

  qmean = qmean + qnew*dtlr

!before changing depth make sure that the borders have been received
  if (pid .eq. 1) then
    call MPI_wait(reqsu, status, ierr)
  elseif (pid .eq. numtasks - 2) then
    call MPI_wait(reqsd, status, ierr)
  elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
    call MPI_wait(reqsu, status, ierr)
    call MPI_wait(reqsd, status, ierr)
  end if

!enddo
end subroutine rivers_dw_flood
