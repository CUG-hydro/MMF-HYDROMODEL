subroutine MOVEQRF(imax, js, je, fd, qrf, area, width)
  integer :: imax, js, je, i, j, ii, jj, iout, jout
  integer, dimension(imax, js:je) :: fd
  real, dimension(imax, js:je) :: qrf, area, width, qrfextra, qrfextra2
  integer :: reqsu2, reqsd2, reqru2, reqrd2

  qrfextra = 0.

  do j = js + 1, je - 1
    do i = 2, imax - 1
      IF (fd(i, j) .gt. 0) then
        if (width(i, j) .lt. 1.) then
          call flowdir(imax, js, je, fd, i, j, iout, jout)
          qrfextra(iout, jout) = qrfextra(iout, jout) + qrf(i, j)*area(i, j)/area(iout, jout)
          qrf(i, j) = 0.
        end if
      END IF
    end do
  end do

  if (numtasks .gt. 1) then
    qrfextra2 = qrfextra
    call sendbordersflood(imax, js, je, qrfextra2, reqsu2, reqsd2, reqru2, reqrd2)
  end if

!make sure that the borders are received before calculating anything
  if (pid .eq. 1) then
    call MPI_wait(reqru2, status, ierr)
  elseif (pid .eq. numtasks - 2) then
    call MPI_wait(reqrd2, status, ierr)
  elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
    call MPI_wait(reqru2, status, ierr)
    call MPI_wait(reqrd2, status, ierr)
  end if

  if (pid .eq. 1) then
    qrfextra(1:imax, je - 1) = qrfextra(1:imax, je - 1) + qrfextra2(1:imax, je - 1)
  elseif (pid .eq. numtasks - 2) then
    qrfextra(1:imax, js + 1) = qrfextra(1:imax, js + 1) + qrfextra2(1:imax, js + 1)
  elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
    qrfextra(1:imax, js + 1) = qrfextra(1:imax, js + 1) + qrfextra2(1:imax, js + 1)
    qrfextra(1:imax, je - 1) = qrfextra(1:imax, je - 1) + qrfextra2(1:imax, je - 1)
  end if

!change qrf
  qrf = qrf + qrfextra

!before changing qrfextra make sure that the borders have been received
  if (pid .eq. 1) then
    call MPI_wait(reqsu2, status, ierr)
  elseif (pid .eq. numtasks - 2) then
    call MPI_wait(reqsd2, status, ierr)
  elseif (pid .gt. 1 .and. pid .lt. numtasks - 2) then
    call MPI_wait(reqsu2, status, ierr)
    call MPI_wait(reqsd2, status, ierr)
  end if

end subroutine moveqrf
