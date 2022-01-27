 ! ifort fulltime_2009_05_15_Fortran.f90 stress.f90 pertvel.f90 cylindervel.f90 intlib3.f90 interp.f90  ellipke2.f90 rk4.f90 Besselk.f90 Besseli.f90 fft.f90 -o output
! bsub -q week -o zresult -R blade output
! bjobs  [Job  ID]


!OUTPUTS INTERPOLATED INTERFACE IN SPHEREVEL and W

!Full solution 2d, 3dlog, matched with far field
!Mapped Interface: Starts Non-Uniform
!Xnclose Xnfar 
!initialize x interface different
!check ustokes cylindervel 
!to not addpoints - take away fillgaps()


!*******Includes third reflection =check cylinder vel *******!



MODULE globalinfo

!experimental parameters
real (kind=8), parameter :: rhot                =1.3612!1.33971!1.37741;
!density of top fluid
real (kind=8), parameter :: rhob                = 1.36583!1.36239!1.37891;
!density of bottom fluid
real (kind=8), parameter :: rhos                = 1.36712!1.46755;
!density of the sphere
real (kind=8), parameter :: mu                  = 10.09!17.59;
!dynamic viscosity of fluid
real (kind=8) :: U                              = 0.0;
!initial velocity of fluid
real (kind=8), parameter :: y0                  = -16 !14.61;
!initial position of interface
real (kind=8), parameter :: R                   = 0.641
!radius of sphere
real (kind=8), parameter :: R0                  = 5.4;
!lateral boundary of the cylinder
real (kind=8), parameter :: maxTime             = 6000;
!time to run to
real (kind=8), parameter :: g                   = 981.0;
!gravity
real (kind=8), parameter :: pi 				    = 3.14159265;

!numerical parameters
real (kind=8), parameter :: dt                  = 4;	!time step
real (kind =8), parameter :: numtrapz         = 0.005;     !h for trapezoidal trapz1trapz1
real (kind=8), parameter :: integthres          = 0.1E-5; !simpson integration
real (kind=8), parameter :: singthres          = 0.1E-6; !threshold singularity x=y
real (kind=8), parameter :: logsing          = 0.1; !threshold log "sing" x=y
real (kind=8), parameter :: logtrapzbig          = 0.5;!h trapz for log around "sing"
real (kind=8), parameter :: logtrapz          = 0.1; ! htrapz for log
real (kind = 8 ) :: FFR =2*R ! radius of full solution, make FFR = 0 for using far field all the time
real (kind =8) ,parameter::  dx = 0.05;


!=======================NON-UNIFORM INTERFACE========================================!

integer (kind=4) ::  RorZ !, XN=ceiling(R0/dx)
real (kind = 8 ) :: R0cl = 2*R
integer (kind=4),parameter  ::       XNfar=20,      XNclose=60
!ceiling(R0/dx)
integer (kind=4) :: XN = XNclose+XNfar
!real (kind=8)::		 dx= R0/(XNclose+XNfar); !initial grid size interface
!integer (kind=4) ::  RorZ, XN=ceiling(R0/dx)




!integer (kind =4), parameter :: numtrapz         = 50;!slow! N pts 4 grid points for trapezoidal (CHECK that trapzN used)
	!dependent parameters
	real (kind=8) :: ms                  		   = 4.0/3.0*pi*R**3*rhos; 				!mass of the sphere
	real (kind=8) :: oneoversixpiamuK     		   = (1-2.10444*(R/R0) +2.08877*((R/R0)**3))* 1.0/(6.0*pi*R*mu)
	real (kind=8) :: buoyancytop         		   = -4.0/3.0*pi*R**3*g*rhot;           !buoyant force when sphere is above the interface
	real (kind=8) :: buoyancybottom         	   = -4.0/3.0*pi*R**3*g*rhob;           !buoyant force when sphere is below the interface
	real (kind=8) :: buoyancyCoeff1      		   = -pi*g/3.0*(rhob-rhot);           	!first coefficient for the buoyant force when sphere is in interface
	real (kind=8) :: buoyancyCoeff2       		   = -2.0*pi*g/3.0*R**3*(rhob+rhot);    !second coefficient for the buoyant force when sphere is in interface
	real (kind=8) :: myt        		   		   = 0
	
	!interface
real (kind=8), dimension(:), allocatable :: x, y, sx, sy, su ,sv ,wu ,wv, cinterpx, cinterpy,cwinterpx, cwinterpy
real (kind=8) yend, xflagb, xflagl, px, py, myrho, myzeta
	integer (kind=4) :: flagb, flagu, flagl,flagt, cinternalcount, cwinternalcount
	
	!fourier components
	real (kind=8), parameter :: epsilon		 		= 0.001
	real (kind=8), parameter ::	LZ                  = 10.0;          !domain window for z-component of cylinder vel
	real (kind=8), parameter ::	LR                  = 10.0;    !domain window for R-component of cylinder vel
	integer (kind=4), parameter :: NZ               = 2**16;                        	!number of points of discretization
	integer (kind=4), parameter :: NR               = 2**14;
	real (kind=8), parameter :: u3coeff             = -2.1044428*R/R0+2.1800173*R**3/R0**3;
	real (kind=8), parameter :: upperZ 				= 40.0;
	integer (kind=4), parameter :: upperRangeZ 		= ceiling((LZ-epsilon)*upperZ/(2.0*pi)+1.0),& 
					upperRangeR 					= ceiling(LR*upperZ/(2.0*pi)+1.0);
!real (kind=8) :: cylindervelR(upperRangeR, ceiling(R0/dx)+1), cylindervelZ(upperRangeZ, ceiling(R0/dx) +1)
	real (kind=8) :: cylindervelR(upperRangeR, XNclose+XNfar+1), cylindervelZ(upperRangeZ, XNclose+XNfar+1)
	real (kind=8) :: zcoordinateZ(upperRangeZ), zcoordinateR(upperRangeR)

	complex, parameter :: MINUS_ONE 				= -1.0
	!real (kind=8) :: myzrangeNR(NR), myzrangeNZ(NZ)
	complex :: imagi								= SQRT(MINUS_ONE)

	real (kind=8) WZ(NZ), WR(NR), myHZ(NZ), myGZ(NZ),&
	myHR(NR), myGR(NR), firstpartZ(NZ), firstpartR(NR), myk(NZ)

real (kind=8) AreaReflux,AreaSpherePortion,AreaEntrain, startx(XNclose+XNfar+1), starty(XNclose+XNfar+1),wforceE, wforceR, wforce, ArchBouyancy,ArchBE, stresspert, stresspertA,ArchBR, stresspertE,stresspertReflux,stresspertcoeff 
END MODULE globalinfo

program fulltime_2009_05_02_Fortran
	use globalinfo
	implicit none
	
	integer (kind=4) i, ierr, flag, ix, iy, iv
	real (kind=8) :: velocity(ceiling(maxTime/dt)+1), stresspertvect(ceiling(maxTime/dt)+1), index, abserr=0.001,& 
    relerr=0.001
	real (kind=8), dimension(:), allocatable :: V, VP
Character(len=65) :: filename
	external rhoode
	
!initialize interface
ALLOCATE(x(XNfar+XNclose+1), STAT=ierr)
!ALLOCATE(x(XN+1), STAT=ierr)
	IF (ierr /= 0) PRINT*, "x : Allocation failed"
	
ALLOCATE(y(XNfar +XNclose+1), STAT=ierr)
!ALLOCATE(y(XN+1), STAT=ierr)
	IF (ierr /= 0) PRINT*, "y : Allocation failed"

DO i=1,XNclose + XNfar+1

if (i<XNclose+2) then
x (i) = (i-1)**(2) *(R0cl/XNclose**(2))
else 
x( i ) = (i - XNclose-1)* (R0 -R0cl ) / XNfar + R0cl
endif 

end do 


!       x = (/(i*R0/(XN), i=0,XN)/)
       y = y0

startx= x;
starty = y;


!IF (ALLOCATED(cinterpx)) DEALLOCATE(cinterpx, STAT=ierr)
!IF (ALLOCATED(cinterpy)) DEALLOCATE(cinterpy, STAT=ierr)


!IF (ALLOCATED(cwinterpx)) DEALLOCATE(cwinterpx, STAT=ierr)
!IF (ALLOCATED(cwinterpy)) DEALLOCATE(cwinterpy, STAT=ierr)


!initialize interface interpolated spherevel
ALLOCATE(cinterpx(1000000), STAT=ierr)
!ALLOCATE(x(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "cintertpx : Allocation failed"

ALLOCATE(cinterpy(1000000), STAT=ierr)
!ALLOCATE(y(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "cintertpy : Allocation failed"

!initialize interface interpolated from w
ALLOCATE(cwinterpx(1000000), STAT=ierr)
!ALLOCATE(x(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "cintertpx : Allocation failed"

ALLOCATE(cwinterpy(1000000), STAT=ierr)
!ALLOCATE(y(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "cintertpy : Allocation failed"


cinternalcount=XN+1;
cwinternalcount=XN+1;




do ix=1, max(cinternalcount,cwinternalcount)
cinterpx(ix)=x(ix);
cinterpy(ix)= y(ix);
cwinterpx(ix)=x(ix);
cwinterpy(ix)= y(ix);
end do


open (unit =9,file = 'VolumeTrack.dat')
write(9,*) "Volume of Entrainment, Volume of Reflux, Volume of Portion of Sphere"


open (unit =8,file = 'WForce.dat')
write(8,*) " ArchBER, Wforce = 6pimuAstresscoeff(wFE - wFR), SphARchBoyancy, ArchBR,wforceR,wforceE"




!	x = (/(i**2*R0/(XN**2), i=0,XN)/)
!	y = y0

!	print *, x(1:100)	
	
	!initialize cylinder velocity 
	call cylindervelinit()

	
	!reset the interface
	!XN = 199.0;
!	IF (ALLOCATED(x)) DEALLOCATE(x,STAT=ierr)
!	IF (ALLOCATED(y)) DEALLOCATE(y,STAT=ierr)
!	ALLOCATE(x(XN+1), STAT=ierr)
!	IF (ierr /= 0) PRINT*, "x : Allocation failed"
!	ALLOCATE(y(XN+1), STAT=ierr)
!	IF (ierr /= 0) PRINT*, "y : Allocation failed"
!	x = (/(i**2*R0/XN**2, i=0,XN)/)
!	y = y0



!**********************  TIME LOOP *******************!
	
	do index = 0,ceiling(maxTime/dt)


		ALLOCATE(V(2*(XN+1)), STAT=ierr)
		IF (ierr /= 0) PRINT*, "V : Allocation failed"

		ALLOCATE(VP(2*(XN+1)), STAT=ierr)
		IF (ierr /= 0) PRINT*, "VP : Allocation failed"


!initialize stokes flow
ALLOCATE(su(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "x : Allocation failed"

ALLOCATE(sv(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "y : Allocation failed"

su = 0
sv = 0



!initialize w flow
ALLOCATE(wu(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "x : Allocation failed"

ALLOCATE(wv(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "y : Allocation failed"

wu = 0
wv = 0
		
		V(1:XN+1) = x
		V(XN+2:2*(XN+1)) = y


!===================== Write  interface ====================== !

WRITE (filename, fmt='(a,f10.2,a)') 'interface',index+1,'.dat'
!print *, filename


open (unit =2,file = filename,form='formatted')
write(2,*) "x,y at time=", myt, "rhos", rhos

do ix=1, XN+1
write(2,*), x(ix), ",", y(ix),","
end do




!===================== Write interpolated interface ================== !


WRITE (filename, fmt='(a,f10.2,a)') 'interpolation',index +1,'.dat'
!print *, "after writing filename"



open (unit =6,file = filename,form='formatted')
!open (unit =6,file = 'interpolation.dat',form='formatted')

!print *, 'before loop ', max(cinternalcount,cwinternalcount)
write(6,*) "interpolated interface x,y at time=", myt

do ix=1, max(cinternalcount,cwinternalcount)
write(6,*), cinterpx(ix), ",", cinterpy(ix),",", cwinterpx(ix), ",", cwinterpy(ix),","
end do



!====================== ODE SOLVER ====================== !



		flag = 1
		call r8_rkf45 (rhoode, 2*(XN+1), V, VP, index*dt, (index+1.0)*dt, relerr, abserr, flag )
		
		x = V(1:XN+1)
		y = V(XN+2:2*(XN+1))



! Write data files



!goto 87
WRITE (filename, fmt='(a,f10.2,a)') 'stokes',index+1,'.dat'

!print *, filename

open (unit =4,file = filename,form='formatted')
write(4,*) "us,sv at time=", myt

do ix=1, XN+1
write(4,*), su(ix), ",", sv(ix), ","
end do


WRITE (filename, fmt='(a,f10.2,a)') 'wpert',index+1,'.dat'

!print *, filename

open (unit =5,file = filename,form='formatted')
write(5,*) "wu,wv at time=", myt

do ix=1, XN+1
write(5,*), wu(ix), ",", wv(ix), ","
end do




!print *, filename

open (unit =9,file = 'VolumeTrack.dat')
write(9,*), AreaEntrain, ",", AreaReflux, ",", AreaSpherePortion


open (unit =8,file = 'WForce.dat')
write(8,*), ArchBE, ",", wforce, ",", ArchBouyancy, ",", ArchBR,",",wforceR,",", wforceE


open (unit =11,file = 'sphereVel.dat')
write(11,*), U, ",", myt, ",", yend





!87 continue


IF (ALLOCATED(V)) DEALLOCATE(V,STAT=ierr)
IF (ALLOCATED(VP)) DEALLOCATE(VP,STAT=ierr)


IF (ALLOCATED(wu)) DEALLOCATE(wu,STAT=ierr)
IF (ALLOCATED(wv)) DEALLOCATE(wv,STAT=ierr)

IF (ALLOCATED(su)) DEALLOCATE(su,STAT=ierr)
IF (ALLOCATED(sv)) DEALLOCATE(sv,STAT=ierr)




		call fillgaps()


		velocity(index+1)=U
		stresspertvect(index+1)=stresspert

	print *, U, ",", myt, ",", yend



	end do
	print *, "velocity"
	do iv = 1, ceiling(maxTime/dt)+1
		print *, velocity(iv), ","
	end do
	
	print *, "stress"
	do iv = 1, ceiling(maxTime/dt)+1
		print *, stresspertvect(iv), ","
	end do
	
	print *, " "
	print *, "*******************************************************************"
	print *, "x"
	do ix=1, XN+1
		print *, x(ix), ","
	end do
	print *, "y"
	do iy=1, XN+1
		print *, y(iy), ","
	end do
	print *, "*******************************************************************"
	print *, " "
end program fulltime_2009_05_02_Fortran

subroutine rhoode(T, V, VP)
	use globalinfo
	implicit none
	
	real (kind=8) :: T, sr(XN+1), V(2*(XN+1)), VP(2*(XN+1)) !, wu(XN+1), wv(XN+1), su(XN+1), sv(XN+1)
	integer (kind=4) ierr, i, ix, iy
	
	myt = T
	
	ALLOCATE(sx(XN+1), STAT=ierr)
	IF (ierr /= 0) PRINT*, "sx : Allocation failed"
	
	ALLOCATE(sy(XN+1), STAT=ierr)
	IF (ierr /= 0) PRINT*, "sy : Allocation failed"

	wu = 0.0
	wv = 0.0
	
	sx = V(1:XN+1)
	sy = V(XN+2:2*(XN+1))

!	print *, sx, ",", sy, ","
!	print *, " "
!	print *, sy
!	print *, " "
	
!	if(sy(XN+1) >= 2.0) then
!		print *, " "
!		print *, "*******************************************************************"
!		print *, "x"
!		do ix=1, XN+1
!			print *, sx(ix), ","
!		end do
!		print *, "y"
!		do iy=1, XN+1
!			print *, sy(iy), ","
!		end do
!		print *, "*******************************************************************"
!		print *, " "
!		stop
!	endif
	!print *, 'sph vel'
	call sphvel()
!	print *, U
	
!	print *, 'sph done'
	call specialpositions(sx, sy)
!print *, flagb
!	print *, 'spec pos done'
	call wTN() !wu, wv)
!	print *, 'w done'
	call stokes() !su, sv)
!	print *, 'stokes done'
!	print *, myt
!	print *, "w", wu, ",", wv
	
	
	sr = sqrt(sx**2+sy**2)
	do i = 1, XN+1
		if (sr(i) <= R) then
			wu(i) = 0.0
			wv(i) = 0.0
			su(i) = 0.0
			sv(i) = 0.0
		endif
	end do
	
	VP(1:XN+1) = su+wu
	VP(XN+2:2*(XN+1)) = sv+wv
	
	IF (ALLOCATED(sx)) DEALLOCATE(sx, STAT=ierr)
	IF (ALLOCATED(sy)) DEALLOCATE(sy, STAT=ierr)
	
end subroutine rhoode

subroutine sphvel()
	use globalinfo
	implicit none
	
	real (kind=8) :: buoyancy, stressbelowsphere, stresssidesphere, stressabovesphere, stressbackflow, &
stressbelowsphereA, stresssidesphereA, stressabovesphereA, stressbackflowA, &
stressbelowsphereE, stresssidesphereE, stressabovesphereE, stressbackflowE,cindex


	real (kind=8), external :: stresstail1D, stressIntegrandFlat1D, stressIntegrandsphere1D, stressIntegrand1D,&
stresstail1DA, stressIntegrandFlat1DA, stressIntegrandsphere1DA, stressIntegrand1DA, &
stresstail1DE, stressIntegrandFlat1DE, stressIntegrandsphere1DE, stressIntegrand1DE


Character(len=45):: filename
	integer i,ix, ierr
	
	!for a two layer fluid only
	if (sy(XN+1)>=R) then
		buoyancy = buoyancybottom;
	elseif (abs(sy(XN+1))<R) then
		buoyancy = buoyancyCoeff1*(3.0*R**2*sy(XN+1)-sy(XN+1)**3)+buoyancyCoeff2;
	else
		buoyancy = buoyancytop;
	endif

!========= Initialize Interpolated  ===== !



IF (ALLOCATED(cinterpx)) DEALLOCATE(cinterpx, STAT=ierr)
IF (ALLOCATED(cinterpy)) DEALLOCATE(cinterpy, STAT=ierr)

!initialize interface
ALLOCATE(cinterpx(1000000), STAT=ierr)
!ALLOCATE(x(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "cintertpx : Allocation failed"

ALLOCATE(cinterpy(1000000), STAT=ierr)
!ALLOCATE(y(XN+1), STAT=ierr)
IF (ierr /= 0) PRINT*, "cintertpy : Allocation failed"

cinternalcount=0.0;
cinterpx =0;
cinterpy = 0;


	!calculate stress force
	stresspert = 0.0;
	if (maxval(sy) > minval(sy)) then
		!print *, "sv 1"
		call specialpositions(sx, sy)
		!print *, "sv 2"
		
		stressbelowsphere = 0.0;
		stresssidesphere = 0.0;
		stressabovesphere = 0.0;
		stressbackflow = 0.0;
		if  (flagu /= 0) then

		call simp(stressIntegrand1D, sy(1), max(-R, sy(1)), integthres, stressbelowsphere)
			call simp(stressIntegrandsphere1D, max(-R, sy(1)), R, integthres, stresssidesphere)
			call simp(stressIntegrand1D, R, yend, integthres, stressabovesphere)
!
!call trapz1(stressIntegrand1D, sy(1), max(-R, sy(1)), numtrapz, stressbelowsphere)
!call trapz1(stressIntegrandsphere1D, max(-R, sy(1)), R, numtrapz, stresssidesphere)
!call trapz1(stressIntegrand1D, R, yend, numtrapz, stressabovesphere)


		elseif  (flagl /= 0) then

call simp(stressIntegrand1D, sy(1), max(-R, sy(1)), integthres, stressbelowsphere)
call simp(stressIntegrandsphere1D, max(-R, sy(1)), yend, integthres, stresssidesphere)

!call trapz1(stressIntegrand1D, sy(1), max(-R, sy(1)), numtrapz, stressbelowsphere)
!call trapz1(stressIntegrandsphere1D, max(-R, sy(1)), yend, numtrapz, stresssidesphere)
!

else
call simp(stressIntegrandFlat1D, sx(1), xflagb, integthres, stressbelowsphere)
!call trapz1(stressIntegrandFlat1D, sx(1), xflagb, numtrapz, stressbelowsphere)

endif

if (flagb /= XN+1) then
!call simp(stresstail1D, xflagb, sx(XN+1), integthres, stressbackflow)
call trapz1(stresstail1D, xflagb, sx(XN+1),real(0.01,kind=8), stressbackflow)

endif
stresspertcoeff=-0.25*g*(rhot-rhob)*R*2.0*pi
stresspert  = stresspertcoeff*(stressbelowsphere + stresssidesphere + stressabovesphere -stressbackflow)
endif






!calculate archimedean force of fluid
stresspertA = 0.0;
if (maxval(sy) > minval(sy)) then
!print *, "sv 1"
call specialpositions(sx, sy)
!print *, "sv 2"

stressbelowsphereA = 0.0;
stresssidesphereA = 0.0;
stressabovesphereA = 0.0;
stressbackflowA = 0.0;
if  (flagu /= 0) then
!		print *, "sy(1)", sy(1)
call simp(stressIntegrand1DA, sy(1), max(-R, sy(1)), integthres, stressbelowsphereA)
call simp(stressIntegrandsphere1DA, max(-R, sy(1)), R, integthres, stresssidesphereA)
call simp(stressIntegrand1DA, R, yend, integthres, stressabovesphereA)

!			call trapz1(stressIntegrand1D, sy(1), max(-R, sy(1)), 0.01, stressbelowsphere)
!                        call trapz1(stressIntegrandsphere1D, max(-R, sy(1)), R, 0.01, stresssidesphere)
!                        call trapz1(stressIntegrand1D, R, yend, 0.01, stressabovesphere)
!print *, "Case 1"

elseif  (flagl /= 0) then
!print *, "sv 3"
call simp(stressIntegrand1DA, sy(1), max(-R, sy(1)), integthres, stressbelowsphereA)
!                call trapz1(stressIntegrand1D, sy(1), max(-R, sy(1)), 0.01, stressbelowsphere)

!print *, "sv 4"
call simp(stressIntegrandsphere1DA, max(-R, sy(1)), yend, integthres, stresssidesphereA)
!        call trapz1(stressIntegrandsphere1D, max(-R, sy(1)), yend, 0.01, stresssidesphere)


!print *, "Case 2"
else
call simp(stressIntegrandFlat1DA, sx(1), xflagb, integthres, stressbelowsphereA)
!                        call trapz1(stressIntegrandFlat1D, sx(1), xflagb, 0.01, stressbelowsphere)

!print *, "Case 3"
endif

if (flagb /= XN+1) then
!print *, "sv 6"
call simp(stresstail1DA, xflagb, sx(XN+1), integthres, stressbackflowA)
!	call trapz1(stresstail1D, xflagb, sx(XN+1), 0.01, stressbackflow)

!print *, "sv 7"
endif
!	print *, stressbelowsphere , stresssidesphere , stressabovesphere, stressbackflow
stresspertA  = -g*(rhob-rhot)*(stressbelowsphereA + stresssidesphereA + stressabovesphereA-stressbackflowA)
endif






!calculate density anomaly force as it scales with shell size epsilon = -sy(1) - A


!calculate archimedean force of fluid
stresspertE = 0.0;
if (maxval(sy) > minval(sy)) then


stressbelowsphereE = 0.0;
stresssidesphereE = 0.0;
stressabovesphereE = 0.0;
stressbackflowE = 0.0;


if  (flagu /= 0) then
!		print *, "sy(1)", sy(1)
call simp(stressIntegrand1DE, sy(1), max(-R, sy(1)), integthres, stressbelowsphereE)
call simp(stressIntegrandsphere1DE, max(-R, sy(1)), R, integthres, stresssidesphereE)
call simp(stressIntegrand1DE, R, yend, integthres, stressabovesphereE)



elseif  (flagl /= 0) then
!print *, "sv 3"
call simp(stressIntegrand1DE, sy(1), max(-R, sy(1)), integthres, stressbelowsphereE)

call simp(stressIntegrandsphere1DE, max(-R, sy(1)), yend, integthres, stresssidesphereE)


!print *, "Case 2"
else
call simp(stressIntegrandFlat1DE, sx(1), xflagb, integthres, stressbelowsphereE)
!                        call trapz1(stressIntegrandFlat1D, sx(1), xflagb, 0.01, stressbelowsphere)

!print *, "Case 3"
endif

if (flagb /= XN+1) then
!print *, "sv 6"
call simp(stresstail1DE, xflagb, sx(XN+1), integthres, stressbackflowE)

endif

stresspertE  = g*(rhob-rhot)*(stressbelowsphereE + stresssidesphereE + stressabovesphereE-stressbackflowE)
endif





wforceE= stressbelowsphere + stresssidesphere + stressabovesphere
wforceR=stressbackflow
wforce = oneoversixpiamuK*stresspert

ArchBouyancy =oneoversixpiamuK*(g*ms + buoyancy)
ArchBE=-stresspertA *oneoversixpiamuK
ArchBR =-g*(rhob-rhot)*(-stressbackflowA)*oneoversixpiamuK
stresspertReflux= -g* ( rhob-rhot)*stressbackflowE


	U = oneoversixpiamuK*(g*ms + buoyancy + stresspert)
	!stresspert = oneoversixpiamuK*(stresspertcoeff*stresssidesphere)/4.102740976860280e-01

end subroutine sphvel

subroutine specialpositions(myx, myy)
	!determine special positions on the interface
	use globalinfo
	implicit none
	 
	real (kind=8) buoyancy
	real (kind=8) :: myx(XN+1), myy(XN+1)
	integer (kind=4) temp(1), tempmaxi
	
	flagl = 0
	flagu = 0
	 
	!find position of backflow
	yend = myy(XN+1)
	
	if(maxval(myy) > yend) then
		flagb = XN+1
		do tempmaxi = 1, XN
			if (myy(tempmaxi+1) > yend .and. myy(tempmaxi) <= yend) then
				flagb = tempmaxi
				exit
			endif
		end do

		tempmaxi = min(flagb-2, XN-3)
		tempmaxi = max(tempmaxi, 1)
		!print *, 'sp 1'
		call interpbridge(5, myy(tempmaxi:tempmaxi+4), myx(tempmaxi:tempmaxi+4),  yend, xflagb)
		!print *, 'sp 2'
	else
	   flagb = XN+1;
	   xflagb = myx(XN+1);
	endif 
	 
	!find x position of bottom of sphere
	
	if (yend >= -R) then	
		temp = minloc(abs(myy+R))
		flagl = temp(1)
		if (myy(1) >= -R) then
			xflagl = 0
		else
			!print *, 'sp 3'
			call interpbridge(min(flagb+2, XN+1), myy(1:min(flagb+2, XN+1)), myx(1:min(flagb+2, XN+1)), -R, xflagl)
			!print *, 'sp 4'
		endif
	endif
	 
     !find x position of top of sphere
	
	if (yend >= R) then	
    
    		!flagu=is 1 if interface is past sphere top
		flagu = 1;
        
		temp = minloc(abs(myy-R))
		flagt = temp(1)
	endif
end subroutine specialpositions






subroutine fillgaps()
	use globalinfo
	implicit none
	
	real (kind=8), dimension(:), allocatable :: newx, newy
	real (kind=8) dist, newpt
	integer (kind=4) internalcount, ierr, xi, posi
		
	!initialize new interface
	ALLOCATE(newx(2*(XN+1)), STAT=ierr)
	IF (ierr /= 0) PRINT*, "newx : Allocation failed"
	
	ALLOCATE(newy(2*(XN+1)), STAT=ierr)
	IF (ierr /= 0) PRINT*, "newy : Allocation failed"
	
	call specialpositions(x, y)
    
    internalcount   = 0.0
    do xi=1,XN
        internalcount       = internalcount+1.0
        newx(internalcount) = x(xi)
        newy(internalcount) = y(xi)
        
        dist = sqrt((x(xi+1)-x(xi))**2+(y(xi+1)-y(xi))**2)

        if  ((sqrt(x(xi)**2+y(xi)**2) < (2*R) .and. dist > dx ) .or. dist >R/2) then
       
            !cubic interpolation to find point to fill gap.
            internalcount           = internalcount+1
            if(x(xi) > 2.0*R) then
                newx(internalcount) = 0.5*(x(xi+1)+x(xi))
                posi = min(XN+1.0, xi+3.0)
                !print *, "fillgap1"
                call interpbridge( 7, x(posi-6.0:posi), y(posi-6.0:posi), 0.5*(x(xi+1)+x(xi)), newpt)
                !print *, "fillgap2"
                newy(internalcount) = newpt
            else
                posi = max(xi-3.0, 1.0)
                if(flagl /= 0.0 .and. xi > flagl) then
                    newy(internalcount) = 0.5*(y(xi+1)+y(xi))

                    !print *, "fillgap3"
                    if (y(xi+1) >= y(xi)) then 
	                    call interpbridge( 7, y(posi:posi+6.0), x(posi:posi+6.0), 0.5*(y(xi+1)+y(xi)), newpt)
	                else
	                	call interpbridge( 7, y(posi+6.0:posi:-1.0), x(posi+6.0:posi:-1.0), 0.5*(y(xi+1)+y(xi)), newpt)
	                endif
	                !print *, "fillgap4"
					if (newpt <= max(x(xi+1), x(xi)) .and. newpt >= min(x(xi+1), x(xi))) then
	                    newx(internalcount) = newpt
	                elseif (((0.5*(y(xi+1)+y(xi)))**2 + (0.5*(x(xi+1)+x(xi)))**2) > R**2) then
	                	newx(internalcount) = 0.5*(x(xi+1)+x(xi))
	                else
	                	newx(internalcount) = sqrt(R**2 - (0.5*(y(xi+1)+y(xi)))**2)
	                endif
                else
                	newx(internalcount) = 0.5*(x(xi+1)+x(xi))
                	call interpbridge( 7, x(posi:posi+6.0), y(posi:posi+6.0), 0.5*(x(xi+1)+x(xi)), newpt)
!                    newx(internalcount) = 0.5*(x(xi+1)+x(xi))
                    !print *, "fillgap5"
!                    call interpbridge( 7, x(posi:posi+6.0), y(posi:posi+6.0), 0.5*(x(xi+1)+x(xi)), newpt)
                    !print *, "fillgap6"

					if (newpt <= max(y(xi+1), y(xi)) .and. newpt >= min(y(xi+1), y(xi))) then
	                    newy(internalcount) = newpt
	                elseif (((0.5*(y(xi+1)+y(xi)))**2 + (0.5*(x(xi+1)+x(xi)))**2) > R**2) then
	                	newy(internalcount) = 0.5*(y(xi+1)+y(xi))
	                else
	                	newy(internalcount) = sqrt(R**2 - (0.5*(x(xi+1)+x(xi)))**2)
	                endif
                endif
            endif
        endif
    end do
    newx(internalcount+1) = x(XN+1)
    newy(internalcount+1) = y(XN+1)

    IF (ALLOCATED(x)) DEALLOCATE(x,STAT=ierr)
	IF (ALLOCATED(y)) DEALLOCATE(y,STAT=ierr)
	
	XN = internalcount
	ALLOCATE(x(XN+1), STAT=ierr)
	IF (ierr /= 0) PRINT*, "fillgap - x : Allocation failed"
	
	ALLOCATE(y(XN+1), STAT=ierr)
	IF (ierr /= 0) PRINT*, "fillgap - y : Allocation failed"
	
	x   = newx(1:internalcount+1)
    y   = newy(1:internalcount+1)
    
    IF (ALLOCATED(newx)) DEALLOCATE(newx,STAT=ierr)
	IF (ALLOCATED(newy)) DEALLOCATE(newy,STAT=ierr)
end subroutine fillgaps

subroutine interpbridge(N, interpx, interpy, xval, yval)
	use globalinfo
	implicit none
	
	integer (kind=4) :: N, setmin(1), mini, maxi, tempi, tempj, interpchecki=1
	real (kind=8) :: interpx(N), interpy(N), d(N), checkorder(N-1)
	real (kind=8) xval, yval, checkmin, checkmax

	if (N == 1) then
		yval = interpy(1)
	else
		checkorder = interpx(2:N) - interpx(1:N-1)
		checkmin = minval(checkorder)
		interpchecki = 1
		if (checkmin <= 0) then
			mini = 1
			maxi = 1
			
			do interpchecki=1, N
				do tempi = maxi, N-1
					if (checkorder(tempi) > 0) then
						mini = tempi
						maxi = N
						do tempj = tempi, N-1
							if (checkorder(tempj) < 0) then
								maxi = tempj
								exit
							endif
						end do
						exit
					endif
				end do
				if (xval<= interpx(maxi) .and. xval >= interpx(mini)) then
					exit
				endif
			end do
		else
			mini = 1
			maxi = N
		endif
		
		if (xval > interpx(maxi) .or. xval < interpx(mini) .or. interpchecki == N) then
			if (xval>interpx(maxi)) then
				print *, "too large"
			elseif (xval < interpx(mini)) then
				print *, "too small"
			else
				print *, "interpchecki", interpchecki, N
			endif
			print *, "out of domain error"
			print *, "time", myt, "xval", xval, "flagl", flagl, "flagb", flagb, "flagu", flagu
			print *, "interpx"
			do interpchecki = 1, N
				print *, interpx(interpchecki)
			end do
			print *, "interpy"
			do interpchecki = 1, N
				print *, interpy(interpchecki)
			end do
			print *, "x"
			do interpchecki = 1, XN+1
				print *, sx(interpchecki)
			end do
			print *, "y"
			do interpchecki = 1, XN+1
				print *, sy(interpchecki)
			end do

			print *, interpx(maxi), xval, interpx(mini)
			stop
		endif		
	
		call spline_pchip_set (maxi-mini+1, interpx(mini:maxi), interpy(mini:maxi), d)
		call spline_pchip_val (maxi-mini+1, interpx(mini:maxi), interpy(mini:maxi), d, 1, xval, yval)
	endif
end subroutine
