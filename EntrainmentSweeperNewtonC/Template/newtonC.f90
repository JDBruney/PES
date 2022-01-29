Subroutine newtonC(t,D, rhot,rhob,rhotN)
!c===========================================================
!c Solution to Newtons Law of Cooling
!c IC: rhob for one side of interface, rhob for other side of interface
!c written by: Dylan Bruney (January 2022)
!c-------------------------------------------------
!c t        - time from when diffusion takes affect
!c D        - diffusion constant
!c rhot     - top density 
!c rhob     - bottom density 
!c rhotN	- Result of evaluation)
!c=================================================
real (kind = 8 ) :: t,D,rhot,rhob,rhotN
rhotN=(rhob-rhot)exp(-D*t)
Return
End 