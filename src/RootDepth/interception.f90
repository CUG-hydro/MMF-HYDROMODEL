subroutine INTERCEPTION(minpprate, precip, lai, intercepstore, ppdrip, pet_i, et_i)
  real :: minpprate, precip, lai, intercepstore, ppdrip, pet_i, et_i
  real :: intercepmax, deficit

  intercepmax = 0.2*lai
  deficit = intercepmax - intercepstore

  if (precip .gt. deficit) then
    if (precip .lt. minpprate) then
      et_i = min(intercepmax, pet_i)
    else
      et_i = 0.
    end if
    intercepstore = intercepmax - et_i
    ppdrip = precip - deficit
  else
    if (precip .lt. minpprate) then
      et_i = min(intercepstore + precip, pet_i)
    else
      et_i = 0.
    end if
    intercepstore = intercepstore + precip - et_i
    ppdrip = 0.
  end if

end subroutine interception
