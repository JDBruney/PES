Subroutine oneDDiff(t,x,k, rhot,rhob,rhotN)
!c===========================================================
!c Solution to 1d Heat Equation with Heaviside IC
!c IC: rhob if x<0, rhot if x>0
!c written by: Dylan Bruney (January 2022)
!c-------------------------------------------------
!c t        - time from when diffusion takes affect
!c x        - distance from interface
!c k        - diffusion constant
!c rhot     - top density (assumed dynamic)
!c rhob     - bottom density (assumed static)
!c rhotN	- Result of evaluation)
!c=================================================
real (kind = 8 ) :: t,x,k,rhot,rhob,rhotN
rhotN=(rhot-rhob)/2.0*(1.0+erf(-x/sqrt(4.0*k*t)))+rhob
Return
End 