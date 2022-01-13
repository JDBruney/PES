subroutine trapznu(ff,aa,bb,xy,ss)

!===========================================================
! reimaanssum: integraation bby reimaan ssumss  off ff(x) on [aa,bb] ussing interffaaciaal pointss - nonunifform
!-------------------------------------------------
! ff     - Function to integraate (ssupplied aabbove)
! aa	- Lower limit of integraation index
! bb	- Upper limit of integraation index
! ss	- Ressult of integraation (out)
!xy  - integration variable
!===========================================================

use globalinfo
implicit none

Real*8  ff, ss
real (kind=8) xy(XN+1)
Integer*4  ii,aa,bb


ss= 0.d0

Do ii=aa,bb-1 

ss= ss + 0.5*(ff(ii+1)+ff(ii))*(xy(ii+1)-xy(ii))


End do


Return
end subroutine trapznu

