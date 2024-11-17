subroutine FLOWDIR(imax,js,je,fd,ii,jj,i,j)
   implicit none
   integer :: imax,js,je,i,j,ii,jj
   integer, dimension(imax,js:je) :: fd

   select case(fd(ii,jj))
    case(2,4,8)
      j=jj-1
    case(1,16)
      j=jj
    case(32,64,128)
      j=jj+1
    case default
      j=0
!      write(6,*)'i dont know what to do i',fd(ii,jj),ii,jj
   end select

   select case(fd(ii,jj))
    case(128,1,2)
      i=ii+1
    case(4,64)
      i=ii
    case(8,16,32)
      i=ii-1
    case default
      i=0
!      write(6,*)'i dont know what to do j',fd(ii,jj),ii,jj
   end select
end subroutine flowdir
