!For Non-Uniform Interface !Full2d & FF
subroutine cylindervelinit()
	use globalinfo
	implicit none
	
	integer index
	do index = 1, NZ
		WZ(index) = real(index-1.0, kind=8)*(LZ-epsilon)/NZ
		!myzrangeNZ(index) =  real(index-1.0, kind=8)*2*pi/LZ
		myk(index) = real(index-1.0, kind =8)*2.0*pi/(LZ-epsilon)
	end do
	
	do index = 1, NR
		WR(index) = real(index-1.0, kind=8)*LR/NR
		!myzrangeNR (index) =  real(index-1.0, kind=8)*2*pi/LR
	end do
	
!	print *, "myk", myk(10)
!	
!	print *, "WZ", WZ(1)
!	print *, "WR", WR(1)
	
	call Hfunc(NZ, WZ+epsilon, myHZ)
	call Gfunc(NZ, WZ+epsilon, myGZ)
	call Hfunc(NR, WR, myHR)
	call Gfunc(NR, WR, myGR)
	
!	print *, "myHZ", myHZ(1)
!	print *, "myGZ", myGZ(1)
!	print *, "myHR", myHR(1)
!	print *, "myGR", myGR(1)
!	stop

	firstpartZ = (WZ+epsilon)/2.0*(myHZ+myGZ);
	firstpartR = WR/2.0*(myHR+myGR);

!	print *, "firstpartZ", firstpartZ(10)
!	print *, "firstpartR", firstpartR(10)	
	do index = 1, upperRangeZ
		zcoordinateZ(index) = (index-1.0)*2.0*pi/(LZ-epsilon)
	end do
	
	do index = 1, upperRangeR
		zcoordinateR(index) = (index-1.0)*2.0*pi/LR
	end do
	
	call cylindervelgrid()
	
!	print *, cylindervelZ
!	stop
end subroutine cylindervelinit

subroutine cylindervelgrid()
	use globalinfo
	implicit none
	
	real (kind=8) :: FZ(NZ), FR(NR), BESSI
	!real (kind=8), external :: sign
	integer myi, WRi, WZi
  	real    ( kind = 4 ) wsavez(4*NZ+15), wsaver(4*NR+15)
  	complex ( kind = 4 ) tempFR(NR), tempFZ(NZ)
  	real (kind = 8) :: besseli0R(NR), besseli1R(NR), besseli0Z(NZ), besseli1Z(NZ)
	

	do myi = 1, XN+1
		!r-component of velocity
		do WRi = 1, NR
			besseli0R(WRi) = BESSI(0,WR(WRi)*x(myi))
			besseli1R(WRi) = BESSI(1,WR(WRi)*x(myi))
		end do
		
		tempFR = real((x(myi)*firstpartR*besseli0R-myGR*besseli1R)*0.5, kind=4)
		tempFR(1) = 0.0
!		print *, firstpartR
!		
!		print *, myGR
!    
!		print *, tempFR(16384)
!		stop
		call cffti ( NR, wsaver )
  		call cfftb ( NR, tempFR, wsaver )
  		FR = real(aimag(tempFR)/NR*LR/pi, kind=8)
  		cylindervelR(1:upperRangeR, myi) = FR(1:upperRangeR)
  		
		!z-component of velocity
		do WZi = 1, NZ
			besseli0Z(WZi) = BESSI(0,(WZ(WZi)+epsilon)*x(myi))
			besseli1Z(WZi) = BESSI(1,(WZ(WZi)+epsilon)*x(myi))
		end do
		
		tempFZ = real((x(myi)*firstpartZ*besseli1Z+myHZ*besseli0Z)*0.5, kind=4)
		call cffti ( NZ, wsavez )
		call cfftb ( NZ, tempFZ, wsavez )
		FZ = real(exp(imagi*epsilon*abs(myk))*real(tempFZ/NZ, kind=8)*(LZ-epsilon)/pi+3.0*R/pi*epsilon*((-2.0*R**2/(3.0*R0**2)+1.0)*x(myi)**2/R0**2+log(epsilon*0.5*R0)-1.0), kind=8)
		cylindervelZ(1:upperRangeZ, myi) = FZ(1:upperRangeZ)
  		
    end do
 !   print *, x(10), y(10)
!    
!    print *, "besseli0R", besseli0R(10)
!    print *, "besseli1R", besseli1R(10)
!    print *, "tempFR", tempFR(10)
!    print *, "tempFZ", tempFZ(10)
!
!    print *, "FZ", FZ(10)
!    print *, "FR", FR(10)
!    
!    print *, upperRangeR, upperRangeZ
!    stop
end subroutine cylindervelgrid

subroutine stokes () !(su, sv)
	use globalinfo
	implicit none
	
	real (kind=8) :: k1(XN+1), k2(XN+1),u3r(XN+1), u3z(XN+1),  myr(XN+1), tempx(4), tempy(4),&
				myinterp1, myinterp2, myinterp3, myinterp4 !, su(XN+1), sv(XN+1)
	real (kind=8), external :: sign
	integer (kind=4) :: starti, sizecyl, i, j
	
	k1 = 0
	k2 = 0
	sizecyl = size(cylindervelR, 2)
	myr  = sqrt(sx**2+sy**2)

	!print *, "stokes 1"
!	print *, "XN", XN
!   print *, "XNclose" , XNclose
!print *, "R0close ", R0close

	do i=1, XN+1

!if( i <XNclose +2) then
!starti = floor ( sqrt(sx(i) *XNclose**2/ R0close ) +1 );
!else
!starti = floor((sx(i) -R0close)*XNfar/ (R0 -R0close) +XNclose+1)
!endif


!	starti = floor(sx(i)/dx);

!for non-uniform interface
if (sx(i)<=R0cl) then
starti = floor((sx(i)*(XNclose**2.0)/R0cl)**(1.0/2.0)+1.0)
else
starti = floor((sx(i)-R0cl)*XNfar/(R0-R0cl)+XNclose +1.0)
endif

		starti = max(starti, 1);
		starti = min(starti, sizecyl-3);
	 
		!horizontal velocity component
		call interpbridge(upperRangeR, zcoordinateR, cylindervelR(1:upperRangeR, starti), abs(sy(i)), myinterp1)
		call interpbridge(upperRangeR, zcoordinateR, cylindervelR(1:upperRangeR, starti+1), abs(sy(i)), myinterp2)
		call interpbridge(upperRangeR, zcoordinateR, cylindervelR(1:upperRangeR, starti+2), abs(sy(i)), myinterp3)
		call interpbridge(upperRangeR, zcoordinateR, cylindervelR(1:upperRangeR, starti+3), abs(sy(i)), myinterp4)
		
!	 	tempx = (/ ((starti-1+j)*dx, j=0,3) /)
        tempx = (/(startx(starti+j), j=0,3)/) !try for non-uniform interface
		tempy = (/ myinterp1, myinterp2, myinterp3, myinterp4 /)
	 	tempy=tempy*sign(sy(i))
	 	
		call interpbridge(4, tempx, tempy, max(min(sx(i), R0), real(0.0, kind=8)), k1(i))
		
		!vertical velocity component
		 
		call interpbridge(upperRangeZ, zcoordinateZ, cylindervelZ(1:upperRangeZ, starti), abs(sy(i)), myinterp1)
		call interpbridge(upperRangeZ, zcoordinateZ, cylindervelZ(1:upperRangeZ, starti+1), abs(sy(i)), myinterp2)
		call interpbridge(upperRangeZ, zcoordinateZ, cylindervelZ(1:upperRangeZ, starti+2), abs(sy(i)), myinterp3)
		call interpbridge(upperRangeZ, zcoordinateZ, cylindervelZ(1:upperRangeZ, starti+3), abs(sy(i)), myinterp4)
		 
		tempy = (/ myinterp1, myinterp2, myinterp3, myinterp4 /)
		
		call interpbridge(4, tempx, tempy, max(min(sx(i), R0), real(0.0, kind=8)), k2(i))
		
	end do
	
	!print *, sx(XN+1), sy(XN+1), k1(XN+1), k2(XN+1)
	!print *, "stokes 2"


!third reflection minus stokes part



u3z =  -(R**2*(12.62665286929866*R0**2*(sx**2 + sy**2)**3*(sx**2 + 2*sy**2) +&
   1.276699789481672*R**6*(3*sx**4 - 24*sx**2*sy**2 + 8*sy**4) +&
   2*R**2*(sx**2 + sy**2)**2*(2.10444214488311*R0**2*(sx**2 - 2*sy**2) -&
 6.260028903991107*(sx**4 + 3*sx**2*sy**2 + 2*sy**4)) -&
   R**4*(sx**2 + sy**2)*(-0.14001148948167197*(sx**4 + 20*sx**2*sy**2 - 16*sy**4) -&
      2.18001729431815*(-2*sx**4 + 2*sx**2*sy**2 + 4*sy**4) +&
      1.1366883*(3*sx**4 - 24*sx**2*sy**2 + 8*sy**4)))*U)/&
(8.*R0**3*(sx**2 + sy**2)**4.5)



u3r =    -(R**2*sx*sy*(-6.38349894740836*R**6*&
      (3*sx**2 - 4*sy**2) +&
     12.62665286929866*R0**2*(sx**2 + sy**2)**3 +&
     2*R**2*(sx**2 + sy**2)**2*&
      (-6.31332643464933*R0**2 -&
        6.260028903991107*(sx**2 + sy**2)) +&
     R**4*(sx**2 + sy**2)*&
      (1.5401263842983917*sx**2 -&
        3.3602757475601273*sy**2 +&
        1.1366883*(23*sx**2 - 12*sy**2) +&
        13.080103765908902*(sx**2 + sy**2)))*U)/&
(8.*R0**3*(sx**2 + sy**2)**4.5)

!Stokes Flow with reflections

su = U*((-0.75*R*sx*sy/myr**3+0.75*R**3*sx*sy/myr**5)-k1) + u3r
sv = U*(1+(-0.75*R/myr-0.75*R*sy**2/myr**3-0.25*R**3/myr**3+0.75*sy**2*R**3/myr**5)-k2)+u3z;

!Stokes Free Space 

!su = U*(1)*((-0.75*R*sx*sy/myr**3+0.75*R**3*sx*sy/myr**5))
!sv = U*(1+(1)*(-0.75*R/myr-0.75*R*sy**2/myr**3-0.25*R**3/myr**3+0.75*sy**2*R**3/myr**5));


end subroutine stokes


!subroutine stokes(su, sv)
!	use globalinfo
!	implicit none
!	
!	real (kind=8) :: su(XN+1), sv(XN+1), FZ(NZ), FR(NR)
!	real (kind=8), external :: sign
!	integer myi
!  	real    ( kind = 4 ) wsavez(4*NZ+15), wsaver(4*NR+15)
!  	complex ( kind = 4 ) tempFR(NR), tempFZ(NZ)
!
!	
!	do myi = 1, XN+1
!		
!		!r-component of velocity
!		call cffti ( NR, wsaver )
!		call uhatR ( sx(myi), WR, FR )
!		FR(1) = 0.0		
!		tempFR = real(FR, kind=4)
!  		call cfftb ( NR, tempFR, wsaver )
!  		FR = real(tempFR/NR*LR/pi, kind=8)
!  		call interpbridge(NR, myzrangeNR, FR, abs(sy(myi)), su(myi))
!  		su(myi) = -1.0*sign(su(myi))*su(myi)
!  		
!		!z-component of velocity
!		call cffti ( NZ, wsavez )
!		call uhatZ ( sx(myi), WZ+epsilon, FZ )
!		tempFZ = real(FZ, kind=4)
!  		call cfftb ( NZ, tempFZ, wsavez )
!  		FZ = real(exp(imagi*epsilon*abs(myk))*real(tempFZ/NZ, kind=8)*(LZ-epsilon)/pi+3.0*R/pi*epsilon*((-2.0*R**2/(3.0*R0**2)+1.0)*sx(myi)**2/R0**2+log(epsilon*0.5*R0)-1.0), kind=8)
!  		call interpbridge(NZ, myzrangeNZ, FZ, abs(sy(myi)), sv(myi))
!  		
!    end do
!
!end subroutine stokes

subroutine Hfunc(Num, lambda, ReturnH)
	use globalinfo
	implicit none
	
	integer index, Num
	real (kind=8) :: lambda(Num), besselk0(Num), besselk1(Num),&
	  	besseli1(Num), besseli2(Num), besseli0(Num), ReturnH(Num),&
		BESSK, BESSI
	
	do index = 1, Num
		besselk0(index) = BESSK(0,R0*lambda(index))
		besselk1(index) = BESSK(1,R0*lambda(index))
		besseli0(index) = BESSI(0,R0*lambda(index))
		besseli1(index) = BESSI(1,R0*lambda(index))
		besseli2(index) = BESSI(2,R0*lambda(index))
	end do
	 
	ReturnH = R*(3.0-(6.0+R**2*lambda**2)*(real(besselk0)*besseli2+besseli1*real(besselk1)))/(besseli0*besseli2-besseli1**2);

!	print *, "besselk0", besselk0(1)
!	print *, "besselk1", besselk1(1)
!	print *, "besseli0", besseli0(1)
!	print *, "besseli1", besseli1(1)
!	print *, "besseli2", besseli2(1)
!	print *, "ReturnH", ReturnH(1)
end subroutine

subroutine Gfunc(Num, lambda, ReturnG)
	use globalinfo
	implicit none
	
	integer index, Num
	real (kind=8) :: lambda(Num), besselk1(Num), besselk2(Num), besseli1(Num), besseli2(Num), besseli0(Num), ReturnG(Num),&
		BESSK, BESSI

	do index = 1, Num
		besselk1(index) = BESSK(1,R0*lambda(index))
		besselk2(index) = BESSK(2,R0*lambda(index))
		besseli0(index) = BESSI(0,R0*lambda(index))
		besseli1(index) = BESSI(1,R0*lambda(index))
		besseli2(index) = BESSI(2,R0*lambda(index))
	end do
	
	ReturnG = R*(-3.0+R**2*lambda**2*(real(besselk1)*besseli1+besseli0*real(besselk2)))/(besseli0*besseli2-besseli1**2);
end

!subroutine uhatR(myR, lambda, F)
!	use globalinfo
!	implicit none
!	
!	integer index
!	real (kind=8) :: myR, lambda(NR), F(NR), besseli0(NR),&
!		besseli1(NR), myG(NR), myH(NR), BESSI
!	
!	do index = 1, NR
!		besseli0(index) = BESSI(0,myR*lambda(index))
!		besseli1(index) = BESSI(1,myR*lambda(index))
!	end do
!		
!	call Gfunc(NR, lambda, myG)
!	call Hfunc(NR, lambda, myH)
!	
!	F = 0.5*(0.5*myR*lambda*(myH+myG)*besseli0-myG*besseli1)
!end subroutine
!
!subroutine uhatZ(myR, lambda, F)
!	use globalinfo
!	implicit none
!	
!	integer index
!	real (kind=8) :: myR, lambda(NZ), F(NZ), besseli0(NZ),&
!		besseli1(NZ), myG(NZ), myH(NZ), BESSI
!	
!	do index = 1, NZ
!		besseli0(index) = BESSI(0,myR*lambda(index))
!		besseli1(index) = BESSI(1,myR*lambda(index))
!	end do
!	
!	call Gfunc(NZ, lambda, myG)
!	call Hfunc(NZ, lambda, myH)
!	
!	F = 0.5*(0.5*myR*lambda*(myH+myG)*besseli1+myH*besseli0)
!end subroutine
!
function sign(val)
	implicit none
	
	real(kind=8) val, sign
	
	if (val< 0) then
		sign = -1.0
	else
		sign = 1.0
	endif
end function
