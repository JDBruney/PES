
MODULE globalinfotest
	
!experimental parameters
	real (kind=8), parameter :: rhot                =1.33971!1.37741; 							
!density of top fluid
	real (kind=8), parameter :: rhob                = 1.36239!1.37891;								
!density of bottom fluid
	real (kind=8), parameter :: rhos                = 1.3651!1.46755; 							
!density of the sphere
	real (kind=8), parameter :: mu                  =6.23!17.59;							
!dynamic viscosity of fluid
	real (kind=8) :: U                   		 = 0.0;								
!initial velocity of fluid
	real (kind=8), parameter :: y0                  = -14.61;  							
!initial position of interface
	real (kind=8), parameter :: R                   = 0.635;								
!radius of sphere
	real (kind=8), parameter :: R0                  = 5.4;								
!lateral boundary of the cylinder
	real (kind=8), parameter :: maxTime             = 1800.0;								
!time to run to
	real (kind=8), parameter :: g                   = 981.0;							
!gravity
	real (kind=8), parameter :: pi 				    = 3.14159265;
	
	!numerical parameters
	real (kind=8), parameter :: dt                  = 4.0;	!time step

integer (kind=4) ::  RorZ !, XN=ceiling(R0/dx)   
real (kind = 8 ) :: R0close = 1.5*R/2
integer (kind=4),parameter  ::       XNfar=20,      XNclose=20
!ceiling(R0/dx)
integer (kind=4) :: XN = XNclose+XNfar
!real (kind=8)::		 dx= R0/(XNclose+XNfar); !initial grid size interface
real (kind =8) ,parameter::  dx = -13452; 
!integer (kind=4) ::  RorZ, XN=ceiling(R0/dx)
	real (kind=8), parameter :: integthres          = 0.1E-6; 				 		
!integration accuracy
  real (kind =8), parameter :: numtrapz         = 0.1;     !h for  grid points for trapezoidal
	
	!dependent parameters
	real (kind=8) :: ms                  		   = 4.0/3.0*pi*R**3*rhos; 				!mass of the sphere
	real (kind=8) :: oneoversixpiamuK     		   = 1.0/(6.0*pi*R*mu*(1.0+2.10444*(R/R0)+4.4286677*(R&
    /R0)**2+7.2309626238083045*(R/R0)**3))
	real (kind=8) :: stresspertcoeff     		   = -0.25*g*(rhot-rhob)*R*2.0*pi;    	!coefficient of the perturbation stress
	real (kind=8) :: stresspert          		   = 0.0;								!perturbation stress
	real (kind=8) :: buoyancytop         		   = -4.0/3.0*pi*R**3*g*rhot;           !buoyant force when sphere is above the interface
	real (kind=8) :: buoyancybottom         	   = -4.0/3.0*pi*R**3*g*rhob;           !buoyant force when sphere is below the interface
	real (kind=8) :: buoyancyCoeff1      		   = -pi*g/3.0*(rhob-rhot);           	!first coefficient for the buoyant force when sphere is in interface
	real (kind=8) :: buoyancyCoeff2       		   = -2.0*pi*g/3.0*R**3*(rhob+rhot);    !second coefficient for the buoyant force when sphere is in interface
	real (kind=8) :: drhogover8mu        		   = (rhob-rhot)*g/(8.0*mu);          	!coefficient for the perturbation flow
	real (kind=8) :: myt        		   		   = 0
	
	!interface
	real (kind=8), dimension(:), allocatable :: x, y, sx, sy
	real (kind=8) yend, xflagb, xflagl, px, py, myrho, myzeta
	integer (kind=4) :: flagb, flagu, flagl,flagt
	
	!fourier components
	real (kind=8), parameter :: epsilon		 		= 0.001		
	real (kind=8), parameter ::	LZ                  = 1.0;                            	!domain window for z-component of cylinder vel
	real (kind=8), parameter ::	LR                  = 1.0;                          	!domain window for R-component of cylinder vel
	integer (kind=4), parameter :: NZ               = 2**16;                        	!number of points of discretization
	integer (kind=4), parameter :: NR               = 2**14;
	real (kind=8), parameter :: u3coeff             = -2.10444*R/R0!+1.04335*R^3/R0^3;
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
    
END MODULE globalinfotest

Program main
!=====================================================
! Evaluate each function, test accuracy and time
!=====================================================
use globalinfotest
implicit none

!double precision f,u 

!Real*8 f, u
Real*8   w2IntegrandR,  w2IntegrandZ
external  w2IntegrandR,  w2IntegrandZ
Real*8 start, finish,theta
Real*8 ellipticE, ellipticK,ellipticE1, ellipticK1,wu(XNfar +XNclose +1),wv(XNfar +XNclose +1)
DOUBLE PRECISION k, k1, kbar, tempK, tempE, DRF, DRD, ex, ey, ez
integer ier
integer (kind=4) i, ierr, flag, ix, iy, iv
	real (kind=8) :: velocity(ceiling(maxTime/dt)+1), stresspertvect(ceiling(maxTime/dt)+1), index, abserr=0.001,& 
    relerr=0.001
	real (kind=8), dimension(:), allocatable :: V, VP
 
    
       !! Initialize interface
       
          	
	!initialize interface
ALLOCATE(x(XNfar+XNclose+1), STAT=ierr)
!ALLOCATE(x(XN+1), STAT=ierr)
	IF (ierr /= 0) PRINT*, "x : Allocation failed"
	
ALLOCATE(y(XNfar +XNclose+1), STAT=ierr)
!ALLOCATE(y(XN+1), STAT=ierr)
	IF (ierr /= 0) PRINT*, "y : Allocation failed"

DO i=1,XNclose + XNfar+1

if (i<XNclose+2) then
x (i) = (i-1)**(2) *(R0close/XNclose**(2)) 
else 
x( i ) = (i - XNclose-1)* (R0 -R0close ) / XNfar + R0close 
endif 

end do 


!       x = (/(i*R0/(XN), i=0,XN)/)
!        y = y0
 
       
             
                   
                   
y=(/-0.700822540750432,&
-0.700494057840482,&
-0.699566221138486,&
-0.697297046203183,&
-0.6916031241635811,&
-0.679672565084966,&
-0.658188529566053,&
-0.623449366459655,&
-0.571435823447248,&
-0.498653975459863,&
-0.401843336484698,&
-0.27979663305419,&
-0.13197490236099,&
0.0416124700333768,&
0.204844518325956,&
0.283208511276595,&
0.302975848716101,&
0.305095834373256,&
0.302574970103542,&
0.299210242567293,&
0.296375141523533,&
0.288753308808814,&
0.287521662083106,&
0.287573101904101,&
0.287616780063165,&
0.287631780188805,&
0.287615779650907,&
0.287596998209435,&
0.287575995820051,&
0.28756395482624,&
0.287555903110654,&
0.28753982635487,&
0.287535873403897,&
0.287520338239892,&
0.287501583681356,&
0.28749722605293,&
0.287457833257147,&
0.287436504562885,&
0.287428905170233,&
0.28717466506654,&
0.287227816336151/)



 sy=y;
 sx=x;   

   
call cpu_time(start)

call specialpositions(x, y)

call wTN(wu,wv)

call cpu_time (finish)


write (*,*) wu(1)
write (*,*) wv(1)
write(*,*)  finish -start

end

subroutine wTN(wu, wv)
	use globalinfotest
	implicit none


	real (kind = 8) wu(XN+1), wv(XN+1), wbackflow(XN+1), wsidesphere(XN+1), wbelowsphere(XN+1), wabovesphere(XN+1),&
		tempbackflow, tempsidesphere, tempbelowsphere, tempabovesphere
	real (kind=8), external :: w2IntegrandBackflow, w2IntegrandBelowSphere, w2IntegrandPartialSphere,&
		w2IntegrandZetaSphere, w2IntegrandZetaVert, w2IntegrandR, w2IntegrandZ
	integer (kind=4) wi
	 
	wbackflow		= real(0.0, kind=8)
	wsidesphere 	= real(0.0, kind=8)
	wbelowsphere 	= real(0.0, kind=8)
	wabovesphere 	= real(0.0, kind=8)
	
	do RorZ = 0, 1

		if (minval(sy) < maxval(sy)) then
      
			do wi = 1, XN+1
            
				px = sx(wi)
				py = sy(wi)
         

				tempbackflow	= real(0.0, kind=8)
				tempsidesphere	= real(0.0, kind=8)
				tempbelowsphere	= real(0.0, kind=8)
				tempabovesphere	= real(0.0, kind=8)
                
             
                
				if(flagb <XN+1.0 ) then
					
					!call trapz1(w2IntegrandBackflow, xflagb, sx(XN+1), numtrapz, tempbackflow)
				    call trapznu(w2IntegrandBackflow, flagb, XN+1, sx, tempbackflow)

                endif
				if (flagu /= 0.0) then
				
					!call trapz1(w2IntegrandPartialSphere, real(0.0, kind=8), xflagl, numtrapz, tempbelowsphere)
                    call trapznu(w2IntegrandPartialSphere, 1, flagl, sx,tempbelowsphere)

                
					!call trapz1(w2IntegrandZetaSphere, max(-R, sy(1)), R, numtrapz, tempsidesphere)
					call trapznu(w2IntegrandZetaSphere, flagl, flagt, sy, tempsidesphere)

                   
					!call trapz1(w2IntegrandZetaVert, R, yend, numtrapz, tempabovesphere)
					call trapznu(w2IntegrandZetaVert, flagt, flagb, sy, tempabovesphere)
               

				
                    elseif (flagl /= 0.0) then
                    
                   
					
					!call trapz1(w2IntegrandPartialSphere, real(0.0, kind=8), xflagl,numtrapz,tempbelowsphere)
                    call trapznu(w2IntegrandPartialSphere, 1, flagl, sx, tempbelowsphere)

                  
                    
					!call trapz1(w2IntegrandZetaSphere, max(-R, sy(1)), yend, numtrapz, tempsidesphere)
					call trapznu(w2IntegrandZetaSphere, flagl, flagb, sy, tempsidesphere)
				else
                    
            !call trapz1(w2IntegrandBelowSphere, max(sx(1), real(0.0, kind=8)), xflagb, numtrapz, tempbelowsphere)	
            call trapznu(w2IntegrandBelowSphere, 1, flagb, sx, tempbelowsphere)	
       
	
				endif
				
				wbackflow(wi) 		= tempbackflow
				wbelowsphere(wi) 	= tempbelowsphere
				wsidesphere(wi) 	= tempsidesphere
				wabovesphere(wi) 	= tempabovesphere
			end do
		endif
		if (RorZ > 0) then
			wu = (-wbackflow+wbelowsphere+wsidesphere+wabovesphere)*drhogover8mu
		else
			wv = (-wbackflow+wbelowsphere+wsidesphere+wabovesphere)*drhogover8mu
			!print *, sx(2), sy(2)
		
!			print *, "wbackflow", wbackflow(1)
!			print *, "wbelowsphere",wbelowsphere(1)
!			print *, "wsidesphere",wsidesphere(1)
!			print *, "wabovesphere",wabovesphere(1)
		endif

	end do

!print *, "wu" , wu
!print *, "wv", wv                       
	
	!print *, "flow components", wbelowsphere(3), wsidesphere(3)
	if(abs(wbelowsphere(1)) > 1000.0) then
print *, "stopped at pertvelTN", wbelowsphere(1)	
	stop
	endif
	!print *, "w7"
end subroutine wTN

function w2IntegrandBackflow(rhoi)
	use globalinfotest
	implicit none
	
    integer(kind=4) rhoi
	real (kind=8)  rho,w2IntegrandBackflow
	real (kind=8), external :: w2IntegrandZeta
	myrho=sx(rhoi);
	
	!call trapz1(w2IntegrandZeta, yend, zcoord, numtrapz, w2IntegrandBackflow)
    call trapznu(w2IntegrandZeta, flagb, rhoi, sy, w2IntegrandBackflow)

end

function w2IntegrandZetaSphere(zetai)
	use globalinfotest
	implicit none
    
	integer(kind=4) zetai,xloweri,temp(1)
	real (kind=8)   zeta,xlower, w2IntegrandZetaSphere
	real (kind=8), external :: w2IntegrandRho
	
    

    xlower  = sqrt(R**2 - sy(zetai)**2)
    temp = minloc(abs(sx-xlower))
    xloweri = temp(1)

myzeta=sy(zetai)
    !call trapz1(w2IntegrandRho, xlower, xupper, numtrapz, w2IntegrandZetaSphere)
    call trapznu(w2IntegrandRho, xloweri, zetai, sx, w2IntegrandZetaSphere)

end

function w2IntegrandPartialSphere(rhoi)
	use globalinfotest
	implicit none
	
       integer(kind=4) rhoi
	real (kind=8) rho, zcoord, w2IntegrandPartialSphere
	real (kind=8), external :: w2IntegrandZeta
myrho=sx(rhoi);

!call trapz1(w2IntegrandZeta, zcoord, -R, numtrapz, w2IntegrandPartialSphere)
call trapznu(w2IntegrandZeta, rhoi, flagl, sy, w2IntegrandPartialSphere)

end

function w2IntegrandZetaVert(zetai)
	use globalinfotest
	implicit none
	
    integer(kind=4) zetai
	real (kind=8)  zeta,w2IntegrandZetaVert
	real (kind=8), external :: w2IntegrandRho

myzeta=sy(zetai);
	!call trapz1(w2IntegrandRho, real(0.0, kind=8), xupper, numtrapz,w2IntegrandZetaVert)
     call trapznu(w2IntegrandRho, 1, zetai, sx,w2IntegrandZetaVert)

end


function w2IntegrandBelowSphere(rhoi)
	use globalinfotest
	implicit none
	
	integer(kind=4) rhoi
    real (kind=8) rho, zcoord, w2IntegrandBelowSphere
	real (kind=8), external :: w2IntegrandZeta
myrho=sx(rhoi);


    !call trapz1(w2IntegrandZeta, zcoord, yend, numtrapz, w2IntegrandBelowSphere)
    call trapznu(w2IntegrandZeta, rhoi, flagb, sy, w2IntegrandBelowSphere)


end




!=====================================================
!Integrands
!=====================================================

function w2IntegrandZeta(zetai)
	use globalinfotest
	implicit none
	
    integer(kind=4) rhoi, zetai
	real (kind=8) zeta, w2IntegrandZeta
	real (kind=8), external :: w2IntegrandR, w2IntegrandZ
    

    
    Real*8 ellipticE, ellipticK,ellipticE1, ellipticK1
	DOUBLE PRECISION k,k1, kbar, tempK, tempE, DRF, DRD, ex, ey, ez
	integer ier
    


    myzeta = sy(zetai);

    !myrho= sx(rhoi)
    

if ( ((py-myzeta)**2+(px-myrho)**2)>integthres*1.0E-01) then

!==========================================================
!inputs for ellipticK, ellipticE,ellipticK1, and ellipticE1
!==========================================================
    	k       = 4.0*px*myrho/((py-myzeta)**2+(px-myrho)**2)

		kbar    = k/(k+1.0)
	
		ex = 0.0
		ey = 1.0-kbar
		ez = 1.0
					
		tempK 	= DRF(ex, ey, ez, ier)
		tempE 	= tempK-1.0/3.0*kbar*DRD(ex, ey, ez, ier)

			 
		ellipticE	= sqrt(1.0+k)*tempE
		ellipticK	= (1.0/sqrt(1.0+k))*tempK
        
endif

        
   k1 = -((-4.0)*R**2.0*px*myrho*(R**4.0+(-2.0)* &
  R**2.0*(px*myrho+py*myzeta)+(px**2.0+py**2.0)*(myrho**2.0+myzeta**2.0))**(-1.0))
        
        		kbar    = k1/(k1+1.0)
	
		ex = 0.0
		ey = 1.0-kbar
		ez = 1.0
					
		tempK 	= DRF(ex, ey, ez, ier)
		tempE 	= tempK-1.0/3.0*kbar*DRD(ex, ey, ez, ier)

			 
		ellipticE1	= sqrt(1.0+k1)*tempE
		ellipticK1	= (1.0/sqrt(1.0+k1))*tempK
        
        	if (RorZ == 1.0) then
		w2IntegrandZeta =w2IntegrandR(ellipticK,ellipticE,ellipticK1,ellipticE1)	
        else
		w2IntegrandZeta =w2IntegrandZ(ellipticK,ellipticE,ellipticK1,ellipticE1)	
        endif
        
!endif
 !else   

		!w2IntegrandZeta = 0
	
    
	!print *, "integrand value", w2IntegrandZeta
	!print *, "eval pts", zeta, px, py
	!print *, "eval pt and result", myzeta,  w2IntegrandZeta, myt
!endif    
end 



function w2IntegrandRho(rhoi)
	use globalinfotest
	implicit none
	
	integer( kind=4 ) rhoi,zetai
	real (kind=8) rho, w2IntegrandRho
	real (kind=8), external :: w2IntegrandR, w2IntegrandZ
    
    
    Real*8 ellipticE, ellipticK,ellipticE1, ellipticK1
	DOUBLE PRECISION k, k1, kbar, tempK, tempE, DRF, DRD, ex, ey, ez
	integer ier

	

	
myrho = sy(rhoi)
!myzeta= sy (zetai)

 if ( ((py-myzeta)**2+(px-myrho)**2)>integthres*1.0E-01) then    
    
!=====================================================
! inputs for ellipticK, ellipticE,ellipticK1, and ellipticE1
!=====================================================
    	k       = 4.0*px*myrho/((py-myzeta)**2+(px-myrho)**2)
        
		kbar    = k/(k+1.0)
	
		ex = 0.0
		ey = 1.0-kbar
		ez = 1.0
					
		tempK 	= DRF(ex, ey, ez, ier)
		tempE 	= tempK-1.0/3.0*kbar*DRD(ex, ey, ez, ier)

			 
		ellipticE	= sqrt(1.0+k)*tempE
		ellipticK	= (1.0/sqrt(1.0+k))*tempK
  endif      
        
  k1 = -((-4.0)*R**2*px*myrho*(R**4.0+(-2.0)* &
  R**2.0*(px*myrho+py*myzeta)+(px**2.0+py**2.0)*(myrho**2.0+myzeta**2.0))**(-1.0));
        
        		kbar    = k1/(k1+1.0)
	
		ex = 0.0
		ey = 1.0-kbar
		ez = 1.0
					
		tempK 	= DRF(ex, ey, ez, ier)
		tempE 	= tempK-1.0/3.0*kbar*DRD(ex, ey, ez, ier)

			 
		ellipticE1	= sqrt(1.0+k1)*tempE
		ellipticK1	= (1.0/sqrt(1.0+k1))*tempK
!print *, "k1", k1        

        if (RorZ == 1.0) then
		w2IntegrandRho = w2IntegrandR(ellipticK,ellipticE,ellipticK1,ellipticE1)	
        else
		w2IntegrandRho = w2IntegrandZ(ellipticK,ellipticE,ellipticK1,ellipticE1)
        	endif
        
!   endif     

! else  
 
!if(w2IntegrandRho /= w2IntegrandRho) then
!print*, "NaN occurs when", "x=", px, "py", py, "rho", myrho, "zeta" , 
!myzeta
!endif
		!w2IntegrandRho = 0
        	
!endif     
end 




function w2IntegrandR(ellipticK,ellipticE,ellipticK1, ellipticE1)
! this is the horizontal component of the velocity
	use globalinfotest
	implicit none
 
 	real (kind=8) w2IntegrandR,ellipticE, ellipticK,ellipticE1, ellipticK1,ILogR
    real (kind=8), external :: I1R,I2R,I3R,I4R,I5R,I6R,I7R,I8R,LogTermThetaR

w2IntegrandR=0

if (px > 0.0) then



if (abs (myzeta*px + py * myrho ) > 0.1) then
!call simp2(LogTermThetaR, real(0.0,kind=8), real(2.0*pi,kind=8),integthres, ILogR)
call trapz1(LogTermThetaR,  real(0.0,kind=8),real(2.0*pi,kind=8), numtrapz, ILogR)
else
call trapz1(LogTermThetaR,  real(0.0,kind=8),real(2.0*pi,kind=8), real(0.5,kind=8), ILogR)
endif





			w2IntegrandR = I1R(ellipticK,ellipticE)+I2R(ellipticK1,ellipticE1)&
            +I3R(ellipticK1,ellipticE1)+I4R(ellipticK1,ellipticE1)&
            +I5R(ellipticK1,ellipticE1)+I6R(ellipticK1,ellipticE1)&
            +I7R(ellipticK1,ellipticE1)+I8R(ellipticK1,ellipticE1)+ILogR
	
		w2IntegrandR = w2IntegrandR*myrho/pi
endif
!print *, "w2IntegrandR", w2IntegrandR

end




function w2IntegrandZ(ellipticK, ellipticE,ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
 	real (kind=8) w2IntegrandZ,ellipticE, ellipticK,ellipticE1, ellipticK1,ILogZ
    real (kind=8), external :: I1Z,I2Z,I3Z,I4Z,I5Z,I6Z,I7Z,I8Z,LogTermThetaZ



if (px ==0.0 .and. myrho ==0.0 ) then




call trapz1(LogTermThetaZ,real(0.0,kind=8),real(2.0*pi,kind=8),numtrapz,ILogZ)

if (ILogZ /= ILogZ) then
ILogZ = 0
!print *, "aparent singularity"

endif

else


if (abs(myzeta*px +py * myrho) > 0.1 )  then  

call trapz1(LogTermThetaZ,real(0.0,kind=8),real(2.0*pi,kind=8),numtrapz,ILogZ)

!print*, "no singularity"
else
call trapz1(LogTermThetaZ,real(0.0,kind=8),real(2.0*pi,kind=8),real(0.5,kind=8), ILogZ)
!print *, "aparent singularity"

endif
endif



		w2IntegrandZ = I1Z(ellipticK,ellipticE)+I2Z(ellipticK1,ellipticE1)&
        +I3Z(ellipticK1,ellipticE1)+I4Z(ellipticK1,ellipticE1)&
        +I5Z(ellipticK1,ellipticE1)+I6Z(ellipticK1,ellipticE1) &
        +I7Z(ellipticK1,ellipticE1)+I8Z(ellipticK1,ellipticE1)+ILogZ
	
		w2IntegrandZ = w2IntegrandZ*myrho/pi
      
    !  print *, "w2IntegrandZ", w2IntegrandZ

	
end



function I1R(ellipticK,ellipticE)
use globalinfotest
implicit none

Real*8 I1R, ellipticE, ellipticK

I1R=0
if ( ((py-myzeta)**2+(px-myrho)**2)>integthres*1.0E-01) then


I1R= 2.0*px**(-1.0)*(px**2.0+(-2.0)*px*myrho+myrho**2.0+(py+(-1.0)*myzeta)**2.0)**(-1.0/2.0)* &
(px**2.0+2.0*px*myrho+myrho**2.0+(py+(-1.0)*myzeta)**2.0)**(-1.0)*(py+(-1.0)*myzeta)*( &
(px**2.0+(-1.0)*myrho**2.0+(-1.0)*(py+(-1.0)*myzeta)**2.0)*ellipticE+(px**2.0+2.0* &
px*myrho+myrho**2.0+(py+(-1.0)*myzeta)**2.0)*ellipticK)

endif

end




function I1Z(ellipticK,ellipticE)
use globalinfotest
implicit none

Real*8 I1Z, ellipticE, ellipticK

I1Z=0

if ( ((py-myzeta)**2+(px-myrho)**2)>integthres*1.0E-01) then


I1Z=4.0*(py+(-1.0)*myzeta)**2.0*(px**2.0+(-2.0)*px*myrho+myrho**2.0+((-1.0)*py+myzeta)**2.0) &
**(-1.0/2.0)*(px**2.0+2.0*px*myrho+myrho**2.0+((-1.0)*py+myzeta)**2.0)**(-1.0)* &
ellipticE+4.0*(px**2.0+(-2.0)*px*myrho+myrho**2.0+((-1.0)*py+myzeta)**2.0)**( &
-1.0/2.0)*ellipticK;
endif

end




function I2R (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 

 Real*8 I2R, ellipticE1, ellipticK1
 Real*8 ::  a
 
 
 a=R;
	        
I2R=     (-2*(a**2/(myrho**2 + myzeta**2))**2.5*&
       (a**2*myzeta - (myrho**2 + myzeta**2)*py)*&
       ((a**4 - 2*a**2*myzeta*py - (myrho**2 + myzeta**2)*(px**2 - py**2))*&
          ellipticE1 - &
         (a**4 + 2*a**2*(myrho*px - myzeta*py) + &
            (myrho**2 + myzeta**2)*(px**2 + py**2))*&
        ellipticK1))/&
     (a**2*px*(a**4 + 2*a**2*(myrho*px - myzeta*py) + &
         (myrho**2 + myzeta**2)*(px**2 + py**2))*&
       Sqrt((a**4 - 2*a**2*(myrho*px + myzeta*py) + &
           (myrho**2 + myzeta**2)*(px**2 + py**2))/(myrho**2 + myzeta**2)))
end





function I2Z(ellipticK1,ellipticE1)
	use globalinfotest
	implicit none
 

 	Real*8 I2Z, ellipticE1, ellipticK1
    Real*8 ::  a
    a=R;
	     
I2Z= 4.0*(a**2.0*(myrho**2.0+myzeta**2.0)**(-1.0))**(1.0/2.0)*((myrho**2.0+myzeta**2.0)**(-1.0) &
*(a**4.0+(-2.0)*a**2.0*(px*myrho+py*myzeta)+(px**2.0+py**2.0)*(myrho**2.0+myzeta**2.0) &
))**(-1.0/2.0)*((-1.0)*a**2.0*(myrho**2.0+myzeta**2.0)**(-2.0)*(myrho**2.0*py+myzeta* &
((-1.0)*a**2.0+py*myzeta))**2.0*(a**4.0+2.0*a**2.0*(px*myrho+(-1.0)*py*myzeta)+( &
px**2.0+py**2.0)*(myrho**2.0+myzeta**2.0))**(-1.0)*ellipticE1+(-1.0)*ellipticK1)

end





function I3R (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
 	
  	Real*8 I3R, ellipticE1, ellipticK1,coeff
      Real*8 ::  a
      
       a=R;
    
coeff=(a**2+(-1.0)*myrho**2+(-1.0)*myzeta**2)*(myrho**2+myzeta**2)**(-1.0/2.0);
  



I3R=      coeff*(  (-2*myzeta*((a**4 - 2*a**2*(myrho*px + myzeta*py) + &
           (myrho**2 + myzeta**2)*(px**2 + py**2))*   &
             ellipticE1 - &
   (a**4 - 2*a**2*myzeta*py + &
             (myrho**2 + myzeta**2)*(px**2 + py**2))*&
           ellipticK1   ))/&
      (a*(myrho**2 + myzeta**2)**2*px*&
        Sqrt((a**4 - 2*a**2*(myrho*px + myzeta*py) + &
            (myrho**2 + myzeta**2)*(px**2 + py**2))/&
          (myrho**2 + myzeta**2))))
return
end



function I3Z (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
   	Real*8 I3Z, ellipticE1, ellipticK1,coeff
       real*8 ::  a
       a=R;

coeff=(a**2+(-1.0)*myrho**2+(-1.0)*myzeta**2)*(myrho**2+myzeta**2)**(-1.0/2.0);

I3Z=coeff*(4.0*a*myzeta**2.0*(myrho**2.0+myzeta**2.0)**(-2.0)*((myrho**2.0+myzeta**2.0)**(-1.0)*( &
  a**4.0+(-2.0)*a**2.0*(px*myrho+py*myzeta)+(px**2.0+py**2.0)*(myrho**2.0+myzeta**2.0))) &
  **(-1.0/2.0)*ellipticK1)
  return

end



function I4R(ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
   	Real*8 I4R, ellipticE1, ellipticK1,coeff,a
    a=R;
    

coeff=(a**2+(-1.0)*myrho**2+(-1.0)*myzeta**2)*(myrho**2+myzeta**2)**(-1.0/2.0);

I4R=-coeff* ((2*a*(-((a**2*myzeta)/(myrho**2 + myzeta**2)) + py)*&
        ((a**4 - 2*a**2*myzeta*py + &
             (myrho**2 + myzeta**2)*(px**2 + py**2))*&
           (a**4 - 2*a**2*(myrho*px + myzeta*py) + &
             (myrho**2 + myzeta**2)*(px**2 + py**2))*&
           ellipticE1 - &
          (a**8 - 4*a**6*myzeta*py - &
             4*a**2*myzeta*(myrho**2 + myzeta**2)*py*&
              (px**2 + py**2) + &
             (myrho**2 + myzeta**2)**2*(px**2 + py**2)**2 + &
             a**4*(-2*myrho**2*(px**2 - py**2) + &
                2*myzeta**2*(px**2 + 3*py**2)))*&
           ellipticK1))/&
      ((myrho**2 + myzeta**2)**2*px*&
        (a**4 + 2*a**2*(myrho*px - myzeta*py) + &
          (myrho**2 + myzeta**2)*(px**2 + py**2))*&
        ((a**4 - 2*a**2*(myrho*px + myzeta*py) + &
             (myrho**2 + myzeta**2)*(px**2 + py**2))/&
           (myrho**2 + myzeta**2))**1.5))
           return

end






function I4Z(ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
    	Real*8 I4Z, ellipticE1, ellipticK1,coeff
       real*8 ::  a
       a=R;

	
coeff=(a**2+(-1.0)*myrho**2+(-1.0)*myzeta**2)*(myrho**2+myzeta**2)**(-1.0/2.0);

I4Z=-coeff* ((4*a**3*myzeta*(-(a**2*myzeta) + (myrho**2 + myzeta**2)*py)*&
        ellipticE1)/((myrho**2 + myzeta**2)**2*(a**4 + 2*a**2*(myrho*px &
        - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
        Sqrt((a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 +&
         myzeta**2)*(px**2 + py**2))/(myrho**2 + myzeta**2))))
return
end



function I5R(ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
   
 
   	Real*8  ellipticE1, ellipticK1,coeff,I5R
    
       real*8 ::  a
       a=R;


coeff=(a**2+(-1.0)*myrho**2+(-1.0)*myzeta**2)*(myrho**2+myzeta**2)**(-1.0/2.0);

I5R= -coeff*((-2*a**3*myzeta*((a**4 - 2*a**2*myzeta*py - (myrho**2 + myzeta**2)*(px**2 &
- py**2))*ellipticE1 - &
   (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticK1))/&
     ((myrho**2 + myzeta**2)**2*px*(a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + &
     myzeta**2)*(px**2 + py**2))*&
       sqrt((a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 +&
        py**2))/(myrho**2 + myzeta**2))))
return
end



function I5Z (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
   	Real*8 I5Z, ellipticE1, ellipticK1,coeff
       real*8 ::  a
       a=R;

	coeff=(a**2+(-1.0)*myrho**2+(-1.0)*myzeta**2)*(myrho**2+myzeta**2)**(-1.0/2.0);

I5Z= -coeff*( (4*a**3*myzeta*(-(a**2*myzeta) + (myrho**2 + myzeta**2)*py)*ellipticE1 )/&
     ((myrho**2 + myzeta**2)**2*(a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
       Sqrt((a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))/(myrho**2 + myzeta**2))) )
end






function I6R (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
 
   	Real*8 I6R, ellipticE1, ellipticK1,coeff
  real*8 ::  a
       a=R;

	
	    coeff=(a**2+(-1.0)*myrho**2+(-1.0)*myzeta**2)*(myrho**2+myzeta**2)**(-1.0/2.0);

I6R= coeff* ((4*myzeta*((a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))* &
          ((a**8 - 4*a**6*myzeta*py - 4*a**2*myzeta*(myrho**2 + myzeta**2)*py*(px**2 + py**2) + &
               (myrho**2 + myzeta**2)**2*(px**2 + py**2)**2 + 2*a**4*(myrho**2*py**2 + myzeta**2*(px**2 + 3*py**2)))*ellipticE1- &
            (a**4 - 2*a**2*myzeta*py + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
             (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticK1 ) + &
         a**2*(-a**2 + myzeta*py)*((a**4 - 2*a**2*myzeta*py + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
             (a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticE1- &
            (a**8 - 4*a**6*myzeta*py - 4*a**2*myzeta*(myrho**2 + myzeta**2)*py*(px**2 + py**2) + &
               (myrho**2 + myzeta**2)**2*(px**2 + py**2)**2 + a**4*(-2*myrho**2*(px**2 - py**2) + 2*myzeta**2*(px**2 + 3*py**2)))*&
             ellipticK1 )))/(a*(myrho**2 + myzeta**2)**3*px*(a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2&
              + myzeta**2)*(px**2 + py**2))*((a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*&
              (px**2 +py**2))/(myrho**2 + myzeta**2))**1.5))


end



function I6Z (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
 
  	Real*8 I6Z, ellipticE1, ellipticK1,coeff
  real*8 ::  a
       a=R;

	
	    coeff=(a**2+(-1.0)*myrho**2+(-1.0)*myzeta**2)*(myrho**2+myzeta**2)**(-1.0/2.0);

I6Z= coeff* (   (-4*a*myzeta**2*((a**4 - (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticE1 &
           + (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticK1&
         ))/((myrho**2 + myzeta**2)**2*(a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
       Sqrt((a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))/(myrho**2 + myzeta**2))))

end






function I7R (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
   	Real*8 I7R, ellipticE1, ellipticK1,coeff
      real*8 ::  a
       a=R;

	coeff=(-1.0/2.0)*(a**2.0+(-1.0)*px**2.0+(-1.0)*py**2.0)*(myrho**2.0+myzeta**2.0)**(-3.0/2.0)*(( &
  -1.0)*a**2.0+myrho**2.0+myzeta**2.0)

I7R=- coeff* (  (2*((-3*myzeta*(a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
            (a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
            (-((a**4 - 2*a**2*myzeta*py - (myrho**2 + myzeta**2)*(px**2 - py**2))*ellipticE1 ) + &
              (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticK1 ))/&
          (myrho**2 + myzeta**2) - a**2*(-((a**2*myzeta)/(myrho**2 + myzeta**2)) + py)*&
          (-((a**8 - 4*a**6*myzeta*py - 4*a**2*myzeta*(myrho**2 + myzeta**2)*py*(px**2 + py**2) + &
                 (myrho**2 + myzeta**2)**2*(px**2 + py**2)**2 + 2*a**4*(myrho**2*(7*px**2 + py**2) +& 
                 myzeta**2*(px**2 + 3*py**2)))*ellipticE1 ) + (a**4 - 2*a**2*myzeta*py + (myrho**2 + &
                 myzeta**2)*(px**2 + py**2))*(a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 +&
                  myzeta**2)*(px**2 + py**2))*ellipticK1 + 2*(myrho**2 + myzeta**2)*px**2*(4*(a**4 -&
                   2*a**2*myzeta*py + (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticE1 - &
               (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
                ellipticK1 ))))/(a*px*(a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)&
                *(px**2 + py**2))**2*((a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)&
                *(px**2 + py**2))/(myrho**2 + myzeta**2))**1.5))

end



function I7Z (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
   	Real*8 I7Z, ellipticE1, ellipticK1,coeff
       real*8 ::  a
       a=R;

		coeff=(-1.0/2.0)*(a**2.0+(-1.0)*px**2.0+(-1.0)*py**2.0)*(myrho**2.0+myzeta**2.0)**(-3.0/2.0)*(( &
  -1.0)*a**2.0+myrho**2.0+myzeta**2.0)

I7Z=- coeff* ((4*(myrho**2 + myzeta**2)**2*(((a**2*(myrho**2 + 4*myzeta**2) - 3*myzeta*(myrho**2 + myzeta**2)*py)*&
           (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
           (a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticE1 )/&
         (a*(myrho**2 + myzeta**2)**3) - (a*(-((a**2*myzeta)/(myrho**2 + myzeta**2)) + py)**2*&
           (4*(a**4 - 2*a**2*myzeta*py + (myrho**2 + myzeta**2)*(px**2 + py**2))*ellipticE1- &
             (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
              ellipticK1 ))/(myrho**2 + myzeta**2)))/((a**4 + 2*a**2*(myrho*px-myzeta*py) + &
              (myrho**2 + myzeta**2)*(px**2 + py**2))**2*((a**4 - 2*a**2*(myrho*px + myzeta*py)&
               + (myrho**2 + myzeta**2)*(px**2 + py**2))/(myrho**2 + myzeta**2))**1.5))

end






function I8R (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none
 
    	Real*8 I8R, ellipticE1, ellipticK1,coeff
 real*8 ::  a
       a=R;
    
	coeff=(-1.0/2.0)*(a**2.0+(-1.0)*px**2.0+(-1.0)*py**2.0)*(myrho**2.0+myzeta**2.0)**(-3.0/2.0)*(( &
  -1.0)*a**2.0+myrho**2.0+myzeta**2.0)

I8R=- coeff* ( (4*myzeta*Sqrt(myrho**2 + myzeta**2)*((a**12 - 5*a**10*myzeta*py - a**8*(myrho**2 + &
        5*myzeta**2)*(px**2 - 2*py**2) + (myrho**2 + myzeta**2)**3*px**2*(px**2 + py**2)**2 + & 
            a**2*myzeta*(myrho**2 + myzeta**2)**2*py*(3*px**4 + 2*px**2*py**2 - py**4) + &
            2*a**6*myzeta*py*(myzeta**2*(7*px**2 - 5*py**2) + myrho**2*(9*px**2 - 3*py**2)) - &
            a**4*(myzeta**4*(5*px**4 + 12*px**2*py**2 - 5*py**4) + 6*myrho**2*myzeta**2*(px**4&
             + 4*px**2*py**2 - py**4) + myrho**4*(px**4 + 12*px**2*py**2 - py**4)))*ellipticE1- &
         (a**12 + a**10*(2*myrho*px - 5*myzeta*py) + (myrho**2 + myzeta**2)**3*px**2*(px**2 + py**2)**2 + &
            a**2*(myrho**2 + myzeta**2)**2*(px**2 + py**2)*(2*myrho*px**3 - myzeta*py*(3*px**2 + py**2)) - &
            a**4*(myrho**2 + myzeta**2)*(px**2 + py**2)*(2*myrho*myzeta*px*py + myrho**2*(px**2 - py**2)&
             - myzeta**2*(px**2 + 5*py**2)) + a**8*(-6*myrho*myzeta*px*py - myrho**2*(px**2 - 2*py**2) +&
              myzeta**2*(px**2 + 10*py**2)) - 2*a**6*(-3*myrho*myzeta**2*px*py**2 + 3*myrho**2*myzeta*py**3 &
              + myzeta**3*py*(2*px**2 + 5*py**2) + myrho**3*(2*px**3 - px*py**2)))*ellipticK1 ))/&
     (a*px*(a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))**2*&
       (a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))**1.5))
       
      
end




function I8Z (ellipticK1, ellipticE1)
	use globalinfotest
	implicit none

    	Real*8 I8Z, ellipticE1, ellipticK1,coeff
       real*8 ::  a
       a=R;

    
	coeff=(-1.0/2.0)*(a**2.0+(-1.0)*px**2.0+(-1.0)*py**2.0)*(myrho**2.0+myzeta**2.0)**(-3.0/2.0)*(( &
  -1.0)*a**2.0+myrho**2.0+myzeta**2.0)

I8Z=- coeff* ( (4*a*myzeta*((-2*myzeta*(a**4 + 2*a**2*(myrho*px - myzeta*py) + &
              (myrho**2 + myzeta**2)*(px**2 + py**2))*&
            (a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
            ellipticE1 )/&
          (myrho**2 + myzeta**2) + (-((a**2*myzeta)/(myrho**2 + myzeta**2)) + py)*&
          (2*(-a**2 + myzeta*py)*(4*(a**4 - 2*a**2*myzeta*py + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
                ellipticE1 - &
               (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
                ellipticK1 ) + &
            ((a**8 - 4*a**6*myzeta*py - 4*a**2*myzeta*(myrho**2 + myzeta**2)*py*(px**2 + py**2) + &
                  (myrho**2 + myzeta**2)**2*(px**2 + py**2)**2 + &
                  2*a**4*(myrho**2*(7*px**2 + py**2) + myzeta**2*(px**2 + 3*py**2)))*&
                ellipticE1 - &
               (a**4 - 2*a**2*myzeta*py + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
                (a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))*&
                ellipticK1 )/a**2&
            )))/((a**4 + 2*a**2*(myrho*px - myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))**2*&
       ((a**4 - 2*a**2*(myrho*px + myzeta*py) + (myrho**2 + myzeta**2)*(px**2 + py**2))/&
          (myrho**2 + myzeta**2))**1.5))
    
end






function LogTermThetaR(theta)
!----------------------------------------
!  Integrand coming from log term to integrate 3D
!----------------------------------------
	use globalinfotest
	implicit none
    
     
real*8   phi,x2,coeffn,LogTermThetaR
real*8  theta,pxx,y1,y2,y3,x1,x3,pxy,ys1,ys2,ys3,pxs,pxxys
 real*8 ::  a
       a=R;




y1 = myrho *cos(theta)
y2= myrho* sin(theta) 
y3 = myzeta
x1 = px 
x2= 0 
x3 = py

 
pxy=(myrho**2.d0+myzeta**2.d0)**(1.d0/2.d0)
pxx=(px**2.d0+py**2.d0)**(1.d0/2.d0)


ys1 = a**2/pxy**2.d0*y1
ys2= a**2/pxy**2.d0*y2
ys3 = a**2/pxy**2.d0*y3 




coeffn=    (-3.0*(a**2 - x1**2 - x2**2 - x3**2)*(a**2 - y1**2 - y2**2 - y3**2))/&
       (2.0*a*(y1**2 + y2**2 + y3**2))



pxs=(ys1**2.d0+ys2**2.d0+ys3**2.d0)**(1.d0/2.d0)

pxxys=((x1-ys1)**2 + (x2-ys2)**2 +(x3-ys3)**2)**(1.d0/2.d0)



LogTermThetaR= coeffn *((((pxs*x1 + pxx*ys1)*(pxs*x3 + pxx*ys3))/ &
        (pxx*(pxs*pxx + x1*ys1 + x2*ys2 + x3*ys3)**2) - &
       (x1*ys3)/(pxx*(pxs*pxx + x1*ys1 + x2*ys2 + x3*ys3)) + &
       ((pxs - pxxys)*(pxs*(x1 - ys1) + pxxys*ys1)*(pxs*(x3 - ys3) + pxxys*ys3))/&
        (pxxys**2*(-pxs**2 + pxs*pxxys + x1*ys1 + x2*ys2 + x3*ys3)**2) + &
       ((x1 - ys1)*(pxs**2*(x3 - ys3) + pxxys**2*ys3))/&
        (pxxys**3*(-pxs**2 + pxs*pxxys + x1*ys1 + x2*ys2 + x3*ys3)))/pxs)
return
end


function LogTermThetaZ(theta)
!----------------------------------------
!  Integrand coming from log term to integrate 3D
!----------------------------------------
	use globalinfotest
	implicit none
    
real*8   phi,x2,coeffn,LogTermThetaZ
real*8  theta,pxx,y1,y2,y3,x1,x3,pxy,ys1,ys2,ys3,pxs,pxxys
 real*8 ::  a
       a=R;


phi=0
y1 = myrho *cos(theta)
y2= myrho* sin(theta) 
y3 = myzeta
x1 = px
x2= 0 
x3 = py

 
pxy=(myrho**2.d0+myzeta**2.d0)**(1.d0/2.d0)
pxx=(px**2.d0+py**2.d0)**(1.d0/2.d0)

ys1 = a**2/pxy**2.d0*y1
ys2= a**2/pxy**2.d0*y2
ys3 = a**2/pxy**2.d0*y3 




coeffn=    (-3.0*(a**2 - x1**2 - x2**2 - x3**2)*(a**2 - y1**2 - y2**2 - y3**2))/&
       (2.0*a*(y1**2 + y2**2 + y3**2))


pxs=(ys1**2+ys2**2+ys3**2)**(1.d0/2.d0)

pxxys=((x1-ys1)**2 + (x2-ys2)**2 +(x3-ys3)**2)**(1.d0/2.d0)


LogTermThetaZ=coeffn*(((pxs*x3 + pxx*ys3)**2/(pxx*(pxs*pxx + x1*ys1 + x2*ys2 + x3*ys3)**2) - &
       (pxs*pxx + x3*ys3)/(pxx*(pxs*pxx + x1*ys1 + x2*ys2 + x3*ys3)) + &
       ((pxs - pxxys)*(pxs*(x3 - ys3) + pxxys*ys3)**2)/ &
        (pxxys**2*(-pxs**2 + pxs*pxxys + x1*ys1 + x2*ys2 + x3*ys3)**2) + &
       (pxs*pxxys**2*(-pxs + pxxys) + pxs**2*(x3 - ys3)**2 + pxxys**2*(x3 - ys3)*ys3)/&
        (pxxys**3*(-pxs**2 + pxs*pxxys + x1*ys1 + x2*ys2 + x3*ys3)))/pxs )



return
end

























  DOUBLE PRECISION FUNCTION DRF (X, Y, Z, IER)
!
!! DRF computes the incomplete or complete elliptic integral of 1st kind.
!
!***PURPOSE  Compute the incomplete or complete elliptic integral of the
!            1st kind.  For X, Y, and Z non-negative and at most one of
!            them zero, RF(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -1/2
!                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X, Y or Z is zero, the integral is complete.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      DOUBLE PRECISION (RF-S, DRF-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     DRF
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the first kind
!          Standard FORTRAN function routine
!          Double precision version
!          The routine calculates an approximation result to
!          DRF(X,Y,Z) = Integral from zero to infinity of
!
!                               -1/2     -1/2     -1/2
!                     (1/2)(t+X)    (t+Y)    (t+Z)    dt,
!
!          where X, Y, and Z are nonnegative and at most one of them
!          is zero.  If one of them  is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling sequence
!          DRF( X, Y, Z, IER )
!
!          Parameters On entry
!          Values assigned by the calling routine
!
!          X      - Double precision, nonnegative variable
!
!          Y      - Double precision, nonnegative variable
!
!          Z      - Double precision, nonnegative variable
!
!
!
!          On Return    (values assigned by the DRF routine)
!
!          DRF     - Double precision approximation to the integral
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine. It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!          X, Y, Z are unaltered.
!
!
!   3.    Error Messages
!
!
!         Value of IER assigned by the DRF routine
!
!                  Value assigned         Error Message Printed
!                  IER = 1                MIN(X,Y,Z)  <  0.0D0
!                      = 2                MIN(X+Y,X+Z,Y+Z)  <  LOLIM
!                      = 3                MAX(X,Y,Z)  >  UPLIM
!
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                   Not less than 5 * (machine minimum).
!
!          UPLIM  - Upper limit of valid arguments
!
!                   Not greater than (machine maximum) / 5.
!
!
!                     Acceptable values for:   LOLIM      UPLIM
!                     IBM 360/370 SERIES   :   3.0D-78     1.0D+75
!                     CDC 6000/7000 SERIES :   1.0D-292    1.0D+321
!                     UNIVAC 1100 SERIES   :   1.0D-307    1.0D+307
!                     CRAY                 :   2.3D-2466   1.09D+2465
!                     VAX 11 SERIES        :   1.5D-38     3.0D+37
!
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!
!
!          ERRTOL - Relative error due to truncation is less than
!                   ERRTOL ** 6 / (4 * (1-ERRTOL)  .
!
!
!
!        The accuracy of the computed approximation to the integral
!        can be controlled by choosing the value of ERRTOL.
!        Truncation of a Taylor series after terms of fifth order
!        introduces an error less than the amount shown in the
!        second column of the following table for each value of
!        ERRTOL in the first column.  In addition to the truncation
!        error there will be round-off error, but in practice the
!        total error from both sources is usually less than the
!        amount given in the table.
!
!
!
!
!
!          Sample choices:  ERRTOL   Relative Truncation
!                                    error less than
!                           1.0D-3    3.0D-19
!                           3.0D-3    2.0D-16
!                           1.0D-2    3.0D-13
!                           3.0D-2    2.0D-10
!                           1.0D-1    3.0D-7
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   DRF Special Comments
!
!
!
!          Check by addition theorem: DRF(X,X+Z,X+W) + DRF(Y,Y+Z,Y+W)
!          = DRF(0,Z,W), where X,Y,Z,W are positive and X * Y = Z * W.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral DRF(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, Z are unaltered.
!
!
!
!          ********************************************************
!
!          WARNING: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!   Special double precision functions via DRF
!
!
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 1st kind
!
!                  -----------------------------------------
!
!
!
!                                             2         2   2
!                  F(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1)
!
!
!                                  2
!                  K(K) = DRF(0,1-K ,1)
!
!
!                         PI/2     2   2      -1/2
!                       = INT  (1-K SIN (PHI) )   D PHI
!                          0
!
!
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 1st kind
!
!                  -----------------------------------------
!
!
!                                          22    2
!                  EL1(X,KC) = X DRF(1,1+KC X ,1+X )
!
!
!                  Lemniscate constant A
!
!                  -----------------------------------------
!
!
!                       1      4 -1/2
!                  A = INT (1-S )    DS = DRF(0,1,2) = DRF(0,2,1)
!                       0
!
!
!
!    -------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Changed calls to XERMSG to standard form, and some
!           editorial changes.  (RWC))
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DRF
  CHARACTER*16 XERN3, XERN4, XERN5, XERN6
  INTEGER IER
  DOUBLE PRECISION LOLIM, UPLIM, EPSLON, ERRTOL, D1MACH
  DOUBLE PRECISION C1, C2, C3, E2, E3, LAMDA
  DOUBLE PRECISION MU, S, X, XN, XNDEV
  DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
   ZNROOT
  LOGICAL FIRST
  SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DRF
!
  if (FIRST) THEN
     ERRTOL = (4.0D0*D1MACH(3))**(1.0D0/6.0D0)
     LOLIM  = 5.0D0 * D1MACH(1)
     UPLIM  = D1MACH(2)/5.0D0
!
     C1 = 1.0D0/24.0D0
     C2 = 3.0D0/44.0D0
     C3 = 1.0D0/14.0D0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  DRF = 0.0D0
  if (MIN(X,Y,Z) < 0.0D0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     call XERMSG ('SLATEC', 'DRF', &
        'MIN(X,Y,Z) < 0 WHERE X = ' // XERN3 // ' Y = ' // XERN4 // &
        ' AND Z = ' // XERN5, 1, 1)
     return
  end if
!
  if (MAX(X,Y,Z) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'DRF', &
        'MAX(X,Y,Z) > UPLIM WHERE X = '  // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6, 3, 1)
     return
  end if
!
  if (MIN(X+Y,X+Z,Y+Z) < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'DRF', &
        'MIN(X+Y,X+Z,Y+Z) < LOLIM WHERE X = ' // XERN3 // &
        ' Y = ' // XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' // &
        XERN6, 2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
  ZN = Z
!
   30 MU = (XN+YN+ZN)/3.0D0
  XNDEV = 2.0D0 - (MU+XN)/MU
  YNDEV = 2.0D0 - (MU+YN)/MU
  ZNDEV = 2.0D0 - (MU+ZN)/MU
  EPSLON = MAX(ABS(XNDEV),ABS(YNDEV),ABS(ZNDEV))
  if (EPSLON < ERRTOL) go to 40
  XNROOT = SQRT(XN)
  YNROOT = SQRT(YN)
  ZNROOT = SQRT(ZN)
  LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
  XN = (XN+LAMDA)*0.250D0
  YN = (YN+LAMDA)*0.250D0
  ZN = (ZN+LAMDA)*0.250D0
  go to 30
!
   40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
  E3 = XNDEV*YNDEV*ZNDEV
  S  = 1.0D0 + (C1*E2-0.10D0-C2*E3)*E2 + C3*E3
  DRF = S/SQRT(MU)
!
  return
end

FUNCTION RF (X, Y, Z, IER)
!
!! RF computes the incomplete or complete elliptic integral of the 1st kind.  
!
!  For X, Y, and Z non-negative and at most one of
!            them zero, RF(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -1/2
!                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X, Y or Z is zero, the integral is complete.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      SINGLE PRECISION (RF-S, DRF-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     RF
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the first kind
!          Standard FORTRAN function routine
!          Single precision version
!          The routine calculates an approximation result to
!          RF(X,Y,Z) = Integral from zero to infinity of
!
!                               -1/2     -1/2     -1/2
!                     (1/2)(t+X)    (t+Y)    (t+Z)    dt,
!
!          where X, Y, and Z are nonnegative and at most one of them
!          is zero.  If one of them is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling Sequence
!          RF( X, Y, Z, IER )
!
!          Parameters on Entry
!          Values assigned by the calling routine
!
!          X      - Single precision, nonnegative variable
!
!          Y      - Single precision, nonnegative variable
!
!          Z      - Single precision, nonnegative variable
!
!
!
!          On Return     (values assigned by the RF routine)
!
!          RF     - Single precision approximation to the integral
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine.  It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!          X, Y, Z are unaltered.
!
!
!   3.    Error Messages
!
!         Value of IER assigned by the RF routine
!
!                  Value assigned         Error Message Printed
!                  IER = 1                MIN(X,Y,Z)  <  0.0E0
!                      = 2                MIN(X+Y,X+Z,Y+Z)  <  LOLIM
!                      = 3                MAX(X,Y,Z)  >  UPLIM
!
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                   Not less than 5 * (machine minimum).
!
!          UPLIM  - Upper limit of valid arguments
!
!                   Not greater than (machine maximum) / 5.
!
!
!                     Acceptable Values For:   LOLIM      UPLIM
!                     IBM 360/370 SERIES   :   3.0E-78     1.0E+75
!                     CDC 6000/7000 SERIES :   1.0E-292    1.0E+321
!                     UNIVAC 1100 SERIES   :   1.0E-37     1.0E+37
!                     CRAY                 :   2.3E-2466   1.09E+2465
!                     VAX 11 SERIES        :   1.5E-38     3.0E+37
!
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!
!
!          ERRTOL - Relative error due to truncation is less than
!                   ERRTOL ** 6 / (4 * (1-ERRTOL)  .
!
!
!
!              The accuracy of the computed approximation to the inte-
!              gral can be controlled by choosing the value of ERRTOL.
!              Truncation of a Taylor series after terms of fifth order
!              introduces an error less than the amount shown in the
!              second column of the following table for each value of
!              ERRTOL in the first column.  In addition to the trunca-
!              tion error there will be round-off error, but in prac-
!              tice the total error from both sources is usually less
!              than the amount given in the table.
!
!
!
!
!
!          Sample Choices:  ERRTOL   Relative Truncation
!                                    error less than
!                           1.0E-3    3.0E-19
!                           3.0E-3    2.0E-16
!                           1.0E-2    3.0E-13
!                           3.0E-2    2.0E-10
!                           1.0E-1    3.0E-7
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   RF Special Comments
!
!
!
!          Check by addition theorem: RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W)
!          = RF(0,Z,W), where X,Y,Z,W are positive and X * Y = Z * W.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral RF(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, and Z are unaltered.
!
!
!
!          ********************************************************
!
!          Warning: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!   Special Functions via RF
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 1st kind
!                  ----------------------------------------------
!
!
!                                            2         2   2
!                  F(PHI,K) = SIN(PHI) RF(COS (PHI),1-K SIN (PHI),1)
!
!
!                                 2
!                  K(K) = RF(0,1-K ,1)
!
!                         PI/2     2   2      -1/2
!                       = INT  (1-K SIN (PHI) )   D PHI
!                          0
!
!
!
!
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 1st kind
!                  ----------------------------------------------
!
!
!                                         22    2
!                  EL1(X,KC) = X RF(1,1+KC X ,1+X )
!
!
!
!
!                  Lemniscate constant A
!                  ---------------------
!
!
!                       1      4 -1/2
!                  A = INT (1-S )    DS = RF(0,1,2) = RF(0,2,1)
!                       0
!
!
!    -------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Changed calls to XERMSG to standard form, and some
!           editorial changes.  (RWC))
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RF
  real RF
  CHARACTER*16 XERN3, XERN4, XERN5, XERN6
  INTEGER IER
  REAL LOLIM, UPLIM, EPSLON, ERRTOL
  REAL C1, C2, C3, E2, E3, LAMDA
  REAL MU, S, X, XN, XNDEV
  REAL XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
   ZNROOT
  LOGICAL FIRST
  SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  RF
!
  if (FIRST) THEN
     ERRTOL = (4.0E0*R1MACH(3))**(1.0E0/6.0E0)
     LOLIM  = 5.0E0 * R1MACH(1)
     UPLIM  = R1MACH(2)/5.0E0
!
     C1 = 1.0E0/24.0E0
     C2 = 3.0E0/44.0E0
     C3 = 1.0E0/14.0E0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  RF = 0.0E0
  if (MIN(X,Y,Z) < 0.0E0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     call XERMSG ('SLATEC', 'RF', &
        'MIN(X,Y,Z) < 0 WHERE X = ' // XERN3 // ' Y = ' // XERN4 // &
        ' AND Z = ' // XERN5, 1, 1)
     return
  end if
!
  if (MAX(X,Y,Z) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'RF', &
        'MAX(X,Y,Z) > UPLIM WHERE X = '  // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6, 3, 1)
     return
  end if
!
  if (MIN(X+Y,X+Z,Y+Z) < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'RF', &
        'MIN(X+Y,X+Z,Y+Z) < LOLIM WHERE X = ' // XERN3 // &
        ' Y = ' // XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' // &
        XERN6, 2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
  ZN = Z
!
   30 MU = (XN+YN+ZN)/3.0E0
  XNDEV = 2.0E0 - (MU+XN)/MU
  YNDEV = 2.0E0 - (MU+YN)/MU
  ZNDEV = 2.0E0 - (MU+ZN)/MU
  EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
  if (EPSLON < ERRTOL) go to 40
  XNROOT =  SQRT(XN)
  YNROOT =  SQRT(YN)
  ZNROOT =  SQRT(ZN)
  LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
  XN = (XN+LAMDA)*0.250E0
  YN = (YN+LAMDA)*0.250E0
  ZN = (ZN+LAMDA)*0.250E0
  go to 30
!
   40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
  E3 = XNDEV*YNDEV*ZNDEV
  S  = 1.0E0 + (C1*E2-0.10E0-C2*E3)*E2 + C3*E3
  RF = S/SQRT(MU)
!
  return
end

DOUBLE PRECISION FUNCTION DRD (X, Y, Z, IER)
!
!! DRD computes the incomplete or complete elliptic integral of 2nd kind.
!
!***PURPOSE  Compute the incomplete or complete elliptic integral of
!            the 2nd kind. For X and Y nonnegative, X+Y and Z positive,
!            DRD(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -3/2
!                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X or Y is zero, the integral is complete.
!
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      DOUBLE PRECISION (RD-S, DRD-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     DRD
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the second kind
!          Standard FORTRAN function routine
!          Double precision version
!          The routine calculates an approximation result to
!          DRD(X,Y,Z) = Integral from zero to infinity of
!                              -1/2     -1/2     -3/2
!                    (3/2)(t+X)    (t+Y)    (t+Z)    dt,
!          where X and Y are nonnegative, X + Y is positive, and Z is
!          positive.  If X or Y is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling Sequence
!
!          DRD( X, Y, Z, IER )
!
!          Parameters On Entry
!          Values assigned by the calling routine
!
!          X      - Double precision, nonnegative variable
!
!          Y      - Double precision, nonnegative variable
!
!                   X + Y is positive
!
!          Z      - Double precision, positive variable
!
!
!
!          On Return    (values assigned by the DRD routine)
!
!          DRD     - Double precision approximation to the integral
!
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine. It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!
!          X, Y, Z are unaltered.
!
!   3.    Error Messages
!
!         Value of IER assigned by the DRD routine
!
!                  Value assigned         Error message printed
!                  IER = 1                MIN(X,Y)  <  0.0D0
!                      = 2                MIN(X + Y, Z )  <  LOLIM
!                      = 3                MAX(X,Y,Z)  >  UPLIM
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y, and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                    Not less  than 2 / (machine maximum) ** (2/3).
!
!          UPLIM  - Upper limit of valid arguments
!
!                 Not greater than (0.1D0 * ERRTOL / machine
!                 minimum) ** (2/3), where ERRTOL is described below.
!                 In the following table it is assumed that ERRTOL will
!                 never be chosen smaller than 1.0D-5.
!
!
!                    Acceptable values for:   LOLIM      UPLIM
!                    IBM 360/370 SERIES   :   6.0D-51     1.0D+48
!                    CDC 6000/7000 SERIES :   5.0D-215    2.0D+191
!                    UNIVAC 1100 SERIES   :   1.0D-205    2.0D+201
!                    CRAY                 :   3.0D-1644   1.69D+1640
!                    VAX 11 SERIES        :   1.0D-25     4.5D+21
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!          ERRTOL    Relative error due to truncation is less than
!                    3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
!
!
!
!        The accuracy of the computed approximation to the integral
!        can be controlled by choosing the value of ERRTOL.
!        Truncation of a Taylor series after terms of fifth order
!        introduces an error less than the amount shown in the
!        second column of the following table for each value of
!        ERRTOL in the first column.  In addition to the truncation
!        error there will be round-off error, but in practice the
!        total error from both sources is usually less than the
!        amount given in the table.
!
!
!
!
!          Sample choices:  ERRTOL   Relative truncation
!                                    error less than
!                           1.0D-3    4.0D-18
!                           3.0D-3    3.0D-15
!                           1.0D-2    4.0D-12
!                           3.0D-2    3.0D-9
!                           1.0D-1    4.0D-6
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   DRD Special Comments
!
!
!
!          Check: DRD(X,Y,Z) + DRD(Y,Z,X) + DRD(Z,X,Y)
!          = 3 / SQRT(X * Y * Z), where X, Y, and Z are positive.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral DRD(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, Z are unaltered.
!
!
!
!          ********************************************************
!
!          WARNING: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!    -------------------------------------------------------------------
!
!
!   Special double precision functions via DRD and DRF
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind
!
!                  -----------------------------------------
!
!
!                                             2         2   2
!                  E(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1) -
!
!                     2      3             2         2   2
!                  -(K/3) SIN (PHI) DRD(COS (PHI),1-K SIN (PHI),1)
!
!
!                                  2        2            2
!                  E(K) = DRF(0,1-K ,1) - (K/3) DRD(0,1-K ,1)
!
!                         PI/2     2   2      1/2
!                       = INT  (1-K SIN (PHI) )  D PHI
!                          0
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind
!
!                  -----------------------------------------
!
!                                               22    2
!                  EL2(X,KC,A,B) = AX DRF(1,1+KC X ,1+X ) +
!
!                                              3          22    2
!                                 +(1/3)(B-A) X DRD(1,1+KC X ,1+X )
!
!
!
!
!                  Legendre form of alternative ELLIPTIC INTEGRAL
!                  of 2nd kind
!
!                  -----------------------------------------
!
!
!
!                            Q     2       2   2  -1/2
!                  D(Q,K) = INT SIN P  (1-K SIN P)     DP
!                            0
!
!
!
!                                     3          2     2   2
!                  D(Q,K) = (1/3) (SIN Q) DRD(COS Q,1-K SIN Q,1)
!
!
!
!
!                  Lemniscate constant  B
!
!                  -----------------------------------------
!
!
!
!
!                       1    2    4 -1/2
!                  B = INT  S (1-S )    DS
!                       0
!
!
!                  B = (1/3) DRD (0,2,1)
!
!
!                  Heuman's LAMBDA function
!
!                  -----------------------------------------
!
!
!
!                  (PI/2) LAMBDA0(A,B) =
!
!                                    2                2
!                 = SIN(B) (DRF(0,COS (A),1)-(1/3) SIN (A) *
!
!                            2               2         2       2
!                  *DRD(0,COS (A),1)) DRF(COS (B),1-COS (A) SIN (B),1)
!
!                            2       3             2
!                  -(1/3) COS (A) SIN (B) DRF(0,COS (A),1) *
!
!                           2         2       2
!                   *DRD(COS (B),1-COS (A) SIN (B),1)
!
!
!
!                  Jacobi ZETA function
!
!                  -----------------------------------------
!
!                             2                 2       2   2
!                  Z(B,K) = (K/3) SIN(B) DRF(COS (B),1-K SIN (B),1)
!
!
!                                       2             2
!                             *DRD(0,1-K ,1)/DRF(0,1-K ,1)
!
!                               2       3           2       2   2
!                            -(K /3) SIN (B) DRD(COS (B),1-K SIN (B),1)
!
!
! ---------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Modify calls to XERMSG to put in standard form.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DRD
  CHARACTER*16 XERN3, XERN4, XERN5, XERN6
  INTEGER IER
  DOUBLE PRECISION LOLIM, TUPLIM, UPLIM, EPSLON, ERRTOL, D1MACH
  DOUBLE PRECISION C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA
  DOUBLE PRECISION MU, POWER4, SIGMA, S1, S2, X, XN, XNDEV
  DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
   ZNROOT
  LOGICAL FIRST
  SAVE ERRTOL, LOLIM, UPLIM, C1, C2, C3, C4, FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DRD
  if (FIRST) THEN
     ERRTOL = (D1MACH(3)/3.0D0)**(1.0D0/6.0D0)
     LOLIM  = 2.0D0/(D1MACH(2))**(2.0D0/3.0D0)
     TUPLIM = D1MACH(1)**(1.0E0/3.0E0)
     TUPLIM = (0.10D0*ERRTOL)**(1.0E0/3.0E0)/TUPLIM
     UPLIM  = TUPLIM**2.0D0
!
     C1 = 3.0D0/14.0D0
     C2 = 1.0D0/6.0D0
     C3 = 9.0D0/22.0D0
     C4 = 3.0D0/26.0D0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  DRD = 0.0D0
  if (  MIN(X,Y) < 0.0D0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     call XERMSG ('SLATEC', 'DRD', &
        'MIN(X,Y) < 0 WHERE X = ' // XERN3 // ' AND Y = ' // &
        XERN4, 1, 1)
     return
  end if
!
  if (MAX(X,Y,Z) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'DRD', &
        'MAX(X,Y,Z) > UPLIM WHERE X = ' // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6, &
        3, 1)
     return
  end if
!
  if (MIN(X+Y,Z) < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'DRD', &
        'MIN(X+Y,Z) < LOLIM WHERE X = ' // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' // XERN6, &
        2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
  ZN = Z
  SIGMA = 0.0D0
  POWER4 = 1.0D0
!
   30 MU = (XN+YN+3.0D0*ZN)*0.20D0
  XNDEV = (MU-XN)/MU
  YNDEV = (MU-YN)/MU
  ZNDEV = (MU-ZN)/MU
  EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
  if (EPSLON < ERRTOL) go to 40
  XNROOT = SQRT(XN)
  YNROOT = SQRT(YN)
  ZNROOT = SQRT(ZN)
  LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
  SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
  POWER4 = POWER4*0.250D0
  XN = (XN+LAMDA)*0.250D0
  YN = (YN+LAMDA)*0.250D0
  ZN = (ZN+LAMDA)*0.250D0
  go to 30
!
   40 EA = XNDEV*YNDEV
  EB = ZNDEV*ZNDEV
  EC = EA - EB
  ED = EA - 6.0D0*EB
  EF = ED + EC + EC
  S1 = ED*(-C1+0.250D0*C3*ED-1.50D0*C4*ZNDEV*EF)
  S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
  DRD = 3.0D0*SIGMA + POWER4*(1.0D0+S1+S2)/(MU*SQRT(MU))
!
  return
end

subroutine XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)

!! XERMSG processes XERROR messages.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERMSG-A)
!***KEYWORDS  ERROR MESSAGE, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
!   XERMSG processes a diagnostic message in a manner determined by the
!   value of LEVEL and the current value of the library error control
!   flag, KONTRL.  See subroutine XSETF for details.
!
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
!
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
!
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
!
!                   call XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
!
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
!
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   call XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
!
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to count out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
!
!    NERR     An integer value that is chosen by the library routine's
!             author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer.
!
!    LEVEL    An integer value in the range 0 to 2 that indicates the
!             level (severity) of the error.  Their meanings are
!
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
!
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
!
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
!
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
!
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!***REVISION HISTORY  (YYMMDD)
!   880101  DATE WRITTEN
!   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!           THERE ARE TWO BASIC CHANGES.
!           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!               OF LOWER CASE.
!   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!           THE PRINCIPAL CHANGES ARE
!           1.  CLARIFY COMMENTS IN THE PROLOGUES
!           2.  RENAME XRPRNT TO XERPRN
!           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!               CHARACTER FOR NEW RECORDS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           CLEAN UP THE CODING.
!   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!           PREFIX.
!   891013  REVISED TO CORRECT COMMENTS.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!           XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERMSG
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
  CHARACTER*8 XLIBR, XSUBR
  CHARACTER*72  TEMP
  CHARACTER*20  LFIRST
!***FIRST EXECUTABLE STATEMENT  XERMSG
  LKNTRL = J4SAVE (2, 0, .FALSE.)
  MAXMES = J4SAVE (4, 0, .FALSE.)
!
!       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!          SHOULD BE PRINTED.
!
!       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!
  if (NERR < -9999999 .OR. NERR > 99999999 .OR. NERR == 0 .OR. &
     LEVEL < -1 .OR. LEVEL > 2) THEN
     call XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // &
        'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// &
        'JOB ABORT DUE TO FATAL ERROR.', 72)
     call XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
     call XERHLT (' ***XERMSG -- INVALID INPUT')
     return
  end if
!
!       RECORD THE MESSAGE.
!
  I = J4SAVE (1, NERR, .TRUE.)
  call XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!       HANDLE PRINT-ONCE WARNING MESSAGES.
!
  if (LEVEL == -1 .AND. KOUNT > 1) RETURN
!
!       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!
  XLIBR  = LIBRAR
  XSUBR  = SUBROU
  LFIRST = MESSG
  LERR   = NERR
  LLEVEL = LEVEL
  call XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
!
  LKNTRL = MAX(-2, MIN(2,LKNTRL))
  MKNTRL = ABS(LKNTRL)
!
!       SKIP PRINTING if THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!       ZERO AND THE ERROR IS NOT FATAL.
!
  if (LEVEL < 2 .AND. LKNTRL == 0) go to 30
  if (LEVEL == 0 .AND. KOUNT > MAXMES) go to 30
  if (LEVEL == 1 .AND. KOUNT > MAXMES .AND. MKNTRL == 1) go to 30
  if (LEVEL == 2 .AND. KOUNT > MAX(1,MAXMES)) go to 30
!
!       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY if CONTROL FLAG
!       IS NOT ZERO.
!
  if (LKNTRL  /=  0) THEN
     TEMP(1:21) = 'MESSAGE FROM ROUTINE '
     I = MIN(LEN(SUBROU), 16)
     TEMP(22:21+I) = SUBROU(1:I)
     TEMP(22+I:33+I) = ' IN LIBRARY '
     LTEMP = 33 + I
     I = MIN(LEN(LIBRAR), 16)
     TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
     TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
     LTEMP = LTEMP + I + 1
     call XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  end if
!
!       if LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!       1.  LEVEL OF THE MESSAGE
!              'INFORMATIVE MESSAGE'
!              'POTENTIALLY RECOVERABLE ERROR'
!              'FATAL ERROR'
!       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!              'PROG CONTINUES'
!              'PROG ABORTED'
!       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!              'TRACEBACK REQUESTED'
!              'TRACEBACK NOT REQUESTED'
!       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!       EXCEED 74 CHARACTERS.
!       WE SKIP THE NEXT BLOCK if THE INTRODUCTORY LINE IS NOT NEEDED.
!
  if (LKNTRL  >  0) THEN
!
!       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
     if (LEVEL  <=  0) THEN
        TEMP(1:20) = 'INFORMATIVE MESSAGE,'
        LTEMP = 20
     ELSEIF (LEVEL  ==  1) THEN
        TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
        LTEMP = 30
     ELSE
        TEMP(1:12) = 'FATAL ERROR,'
        LTEMP = 12
     ENDIF
!
!       THEN WHETHER THE PROGRAM WILL CONTINUE.
!
     if ((MKNTRL == 2 .AND. LEVEL >= 1) .OR. &
         (MKNTRL == 1 .AND. LEVEL == 2)) THEN
        TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
        LTEMP = LTEMP + 14
     ELSE
        TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
        LTEMP = LTEMP + 16
     ENDIF
!
!       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
     if (LKNTRL  >  0) THEN
        TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
        LTEMP = LTEMP + 20
     ELSE
        TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
        LTEMP = LTEMP + 24
     ENDIF
     call XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  end if
!
!       NOW SEND OUT THE MESSAGE.
!
  call XERPRN (' *  ', -1, MESSG, 72)
!
!       if LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!          TRACEBACK.
!
  if (LKNTRL  >  0) THEN
     WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
     DO 10 I=16,22
        if (TEMP(I:I)  /=  ' ') go to 20
   10    CONTINUE
!
   20    call XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
     call FDUMP
  end if
!
!       if LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
  if (LKNTRL  /=  0) THEN
     call XERPRN (' *  ', -1, ' ', 72)
     call XERPRN (' ***', -1, 'END OF MESSAGE', 72)
     call XERPRN ('    ',  0, ' ', 72)
  end if
!
!       if THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 if (LEVEL <= 0 .OR. (LEVEL == 1 .AND. MKNTRL <= 1)) RETURN
!
!       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!       SUMMARY if THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
  if (LKNTRL > 0 .AND. KOUNT < MAX(1,MAXMES)) THEN
     if (LEVEL  ==  1) THEN
        call XERPRN &
           (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
     ELSE
        call XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
     ENDIF
     call XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
     call XERHLT (' ')
  ELSE
     call XERHLT (MESSG)
  end if
  return
end

subroutine XERPRN (PREFIX, NPREF, MESSG, NWRAP)
!
!! XERPRN prints XERROR error messages processed by XERMSG.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERPRN-A)
!***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
! This routine sends one or more lines to each of the (up to five)
! logical units to which error messages are to be sent.  This routine
! is called several times by XERMSG, sometimes with a single line to
! print and sometimes with a (potentially very long) message that may
! wrap around into multiple lines.
!
! PREFIX  Input argument of type CHARACTER.  This argument contains
!         characters to be put at the beginning of each line before
!         the body of the message.  No more than 16 characters of
!         PREFIX will be used.
!
! NPREF   Input argument of type INTEGER.  This argument is the number
!         of characters to use from PREFIX.  If it is negative, the
!         intrinsic function LEN is used to determine its length.  If
!         it is zero, PREFIX is not used.  If it exceeds 16 or if
!         LEN(PREFIX) exceeds 16, only the first 16 characters will be
!         used.  If NPREF is positive and the length of PREFIX is less
!         than NPREF, a copy of PREFIX extended with blanks to length
!         NPREF will be used.
!
! MESSG   Input argument of type CHARACTER.  This is the text of a
!         message to be printed.  If it is a long message, it will be
!         broken into pieces for printing on multiple lines.  Each line
!         will start with the appropriate prefix and be followed by a
!         piece of the message.  NWRAP is the number of characters per
!         piece; that is, after each NWRAP characters, we break and
!         start a new line.  In addition the characters '$$' embedded
!         in MESSG are a sentinel for a new line.  The counting of
!         characters up to NWRAP starts over for each new line.  The
!         value of NWRAP typically used by XERMSG is 72 since many
!         older error messages in the SLATEC Library are laid out to
!         rely on wrap-around every 72 characters.
!
! NWRAP   Input argument of type INTEGER.  This gives the maximum size
!         piece into which to break MESSG for printing on multiple
!         lines.  An embedded '$$' ends a line, and the count restarts
!         at the following character.  If a line break does not occur
!         on a blank (it would split a word) that word is moved to the
!         next line.  Values of NWRAP less than 16 will be treated as
!         16.  Values of NWRAP greater than 132 will be treated as 132.
!         The actual line length will be NPREF + NWRAP after NPREF has
!         been adjusted to fall between 0 and 16 and NWRAP has been
!         adjusted to fall between 16 and 132.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   880621  DATE WRITTEN
!   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
!           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
!           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
!           SLASH CHARACTER IN FORMAT STATEMENTS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
!           LINES TO BE PRINTED.
!   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
!           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
!   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Added code to break messages between words.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERPRN
  CHARACTER*(*) PREFIX, MESSG
  INTEGER NPREF, NWRAP
  CHARACTER*148 CBUFF
  INTEGER IU(5), NUNIT
  CHARACTER*2 NEWLIN
  PARAMETER (NEWLIN = '$$')
!***FIRST EXECUTABLE STATEMENT  XERPRN
  call XGETUA(IU,NUNIT)
!
!       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
!       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
!       ERROR MESSAGE UNIT.
!
  N = I1MACH(4)
  DO 10 I=1,NUNIT
     if (IU(I)  ==  0) IU(I) = N
   10 CONTINUE
!
!       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
!       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
!       THE REST OF THIS ROUTINE.
!
  if ( NPREF  <  0 ) THEN
     LPREF = LEN(PREFIX)
  ELSE
     LPREF = NPREF
  end if
  LPREF = MIN(16, LPREF)
  if (LPREF  /=  0) CBUFF(1:LPREF) = PREFIX
!
!       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
!       TIME FROM MESSG TO PRINT ON ONE LINE.
!
  LWRAP = MAX(16, MIN(132, NWRAP))
!
!       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
!
  LENMSG = LEN(MESSG)
  N = LENMSG
  DO 20 I=1,N
     if (MESSG(LENMSG:LENMSG)  /=  ' ') go to 30
     LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
!
!       if THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
!
  if (LENMSG  ==  0) THEN
     CBUFF(LPREF+1:LPREF+1) = ' '
     DO 40 I=1,NUNIT
        WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
     return
  end if
!
!       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
!       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
!       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
!       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
!
!       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
!       INDEX INTRINSIC FUNCTION RETURNS ZERO if THERE IS NO OCCURRENCE
!       OR if THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
!       OF THE SECOND ARGUMENT.
!
!       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
!       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
!       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
!       POSITION NEXTC.
!
!       LPIECE  ==  0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
!                       REMAINDER OF THE CHARACTER STRING.  LPIECE
!                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
!                       WHICHEVER IS LESS.
!
!       LPIECE  ==  1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
!                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
!                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
!                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
!                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
!                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
!                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
!                       SHOULD BE INCREMENTED BY 2.
!
!       LPIECE  >  LWRAP+1  REDUCE LPIECE TO LWRAP.
!
!       ELSE            THIS LAST CASE MEANS 2  <=  LPIECE  <=  LWRAP+1
!                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
!                       PROPERLY HANDLES THE END CASE WHERE LPIECE  ==
!                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
!                       AT THE END OF A LINE.
!
  NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
  if (LPIECE  ==  0) THEN
!
!       THERE WAS NO NEW LINE SENTINEL FOUND.
!
     IDELTA = 0
     LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
     if (LPIECE  <  LENMSG+1-NEXTC) THEN
        DO 52 I=LPIECE+1,2,-1
           if (MESSG(NEXTC+I-1:NEXTC+I-1)  ==  ' ') THEN
              LPIECE = I-1
              IDELTA = 1
              GOTO 54
           ENDIF
   52       CONTINUE
     ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC = NEXTC + LPIECE + IDELTA
  ELSEIF (LPIECE  ==  1) THEN
!
!       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
!       DON'T PRINT A BLANK LINE.
!
     NEXTC = NEXTC + 2
     go to 50
  ELSEIF (LPIECE  >  LWRAP+1) THEN
!
!       LPIECE SHOULD BE SET DOWN TO LWRAP.
!
     IDELTA = 0
     LPIECE = LWRAP
     DO 56 I=LPIECE+1,2,-1
        if (MESSG(NEXTC+I-1:NEXTC+I-1)  ==  ' ') THEN
           LPIECE = I-1
           IDELTA = 1
           GOTO 58
        ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC = NEXTC + LPIECE + IDELTA
  ELSE
!
!       if WE ARRIVE HERE, IT MEANS 2  <=  LPIECE  <=  LWRAP+1.
!       WE SHOULD DECREMENT LPIECE BY ONE.
!
     LPIECE = LPIECE - 1
     CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC  = NEXTC + LPIECE + 2
  end if
!
!       PRINT
!
  DO 60 I=1,NUNIT
     WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
!
  if (NEXTC  <=  LENMSG) go to 50
  return
end
subroutine XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, &
     ICOUNT)
!
!! XERSVE records that an XERROR error has occurred.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (XERSVE-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!
!        call XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!
! *Arguments:
!
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or zero if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
!
! *Description:
!
!   Record that this error occurred and possibly dump and clear the
!   tables.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   800319  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900413  Routine modified to remove reference to KFLAG.  (WRB)
!   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!           sequence, use IF-THEN-ELSE, make number of saved entries
!           easily changeable, changed routine name from XERSAV to
!           XERSVE.  (RWC)
!   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERSVE
  PARAMETER (LENTAB=10)
  INTEGER LUN(5)
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
  CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
  CHARACTER*20 MESTAB(LENTAB), MES
  DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
  SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
  DATA KOUNTX/0/, NMSG/0/
!***FIRST EXECUTABLE STATEMENT  XERSVE
!
  if (KFLAG <= 0) THEN
!
!        Dump the table.
!
     if (NMSG == 0) RETURN
!
!        Print to each unit.
!
     call XGETUA (LUN, NUNIT)
     DO 20 KUNIT = 1,NUNIT
        IUNIT = LUN(KUNIT)
        if (IUNIT == 0) IUNIT = I1MACH(4)
!
!           Print the table header.
!
        WRITE (IUNIT,9000)
!
!           Print body of table.
!
        DO 10 I = 1,NMSG
           WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), &
              NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
!
!           Print number of other errors.
!
        if (KOUNTX /= 0) WRITE (IUNIT,9020) KOUNTX
        WRITE (IUNIT,9030)
   20    CONTINUE
!
!        Clear the error tables.
!
     if (KFLAG == 0) THEN
        NMSG = 0
        KOUNTX = 0
     ENDIF
  ELSE
!
!        PROCESS A MESSAGE...
!        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!
     LIB = LIBRAR
     SUB = SUBROU
     MES = MESSG
     DO 30 I = 1,NMSG
        if (LIB == LIBTAB(I) .AND. SUB == SUBTAB(I) .AND. &
           MES == MESTAB(I) .AND. NERR == NERTAB(I) .AND. &
           LEVEL == LEVTAB(I)) THEN
              KOUNT(I) = KOUNT(I) + 1
              ICOUNT = KOUNT(I)
              return
        ENDIF
   30    CONTINUE
!
     if (NMSG < LENTAB) THEN
!
!           Empty slot found for new message.
!
        NMSG = NMSG + 1
        LIBTAB(I) = LIB
        SUBTAB(I) = SUB
        MESTAB(I) = MES
        NERTAB(I) = NERR
        LEVTAB(I) = LEVEL
        KOUNT (I) = 1
        ICOUNT    = 1
     ELSE
!
!           Table is full.
!
        KOUNTX = KOUNTX+1
        ICOUNT = 0
     ENDIF
  end if
  return
!
!     Formats.
!
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / &
     ' LIBRARY    SUBROUTINE MESSAGE START             NERR', &
     '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
end

subroutine XERHLT (MESSG)
!
!! XERHLT aborts program execution after printing XERROR error message.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERHLT-A)
!***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        ***Note*** machine dependent routine
!        XERHLT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     Description of Parameters
!        MESSG is as in XERMSG.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to delete length of character
!           and changed routine name from XERABT to XERHLT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERHLT
  CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERHLT
  STOP
end

subroutine XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
!
!! XERCNT allows user control over handling of XERROR errors.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERCNT-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCNT.
!        If the user has provided his own version of XERCNT, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     Description of Parameters
!
!      --Input--
!        LIBRAR - the library that the routine is in.
!        SUBROU - the subroutine that XERMSG is being called from
!        MESSG  - the first 20 characters of the error message.
!        NERR   - same as in the call to XERMSG.
!        LEVEL  - same as in the call to XERMSG.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
!           names, changed routine name from XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERCNT
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
!***FIRST EXECUTABLE STATEMENT  XERCNT
  return
end

subroutine FDUMP
!
!! FDUMP makes a symbolic dump (should be locally written).
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (FDUMP-A)
!***KEYWORDS  ERROR, XERMSG
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
  return
end

subroutine XGETUA (IUNITA, N)
!
!! XGETUA returns unit numbers to which XERROR messages are sent.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XGETUA-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.
!
!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N  <  5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XGETUA
  DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XGETUA
  N = J4SAVE(5,0,.FALSE.)
  DO 30 I=1,N
     INDEX = I+4
     if (I == 1) INDEX = 3
     IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
  return
end

function J4SAVE (IWHICH, IVALUE, ISET)
!
!! J4SAVE saves or recalls global variables needed by error handling routines.
!
!***LIBRARY   SLATEC (XERROR)
!***TYPE      INTEGER (J4SAVE-I)
!***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
!***SEE ALSO  XERMSG
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900205  Minor modifications to prologue.  (WRB)
!   900402  Added TYPE section.  (WRB)
!   910411  Added KEYWORDS section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  J4SAVE
  LOGICAL ISET
  INTEGER IPARAM(9)
  SAVE IPARAM
  DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
  DATA IPARAM(5)/1/
  DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
  J4SAVE = IPARAM(IWHICH)
  if (ISET) IPARAM(IWHICH) = IVALUE
  return
end

FUNCTION D1MACH (I)
!
!! D1MACH returns floating point machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        D = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of D above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH( 3) = B**(-T), the smallest relative spacing.
!   D1MACH( 4) = B**(1-T), the largest relative spacing.
!   D1MACH( 5) = LOG10(B)
!
!   Assume double precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0  <=  X(I)  <  B for I=1,...,T, 0  <  X(1), and
!   EMIN  <=  E  <=  EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(14) = T, the number of base-B digits.
!   I1MACH(15) = EMIN, the smallest exponent E.
!   I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   900911  Added SUN 386i constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!***END PROLOGUE  D1MACH
!
  double precision d1mach
  INTEGER SMALL(4)
  INTEGER LARGE(4)
  INTEGER RIGHT(4)
  INTEGER DIVER(4)
  INTEGER LOG10(4)
!
  DOUBLE PRECISION DMACH(5)
  SAVE DMACH
!
  EQUIVALENCE (DMACH(1),SMALL(1))
  EQUIVALENCE (DMACH(2),LARGE(1))
  EQUIVALENCE (DMACH(3),RIGHT(1))
  EQUIVALENCE (DMACH(4),DIVER(1))
  EQUIVALENCE (DMACH(5),LOG10(1))
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
!     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
!     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
!     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA SMALL(1) / ZC00800000 /
!     DATA SMALL(2) / Z000000000 /
!     DATA LARGE(1) / ZDFFFFFFFF /
!     DATA LARGE(2) / ZFFFFFFFFF /
!     DATA RIGHT(1) / ZCC5800000 /
!     DATA RIGHT(2) / Z000000000 /
!     DATA DIVER(1) / ZCC6800000 /
!     DATA DIVER(2) / Z000000000 /
!     DATA LOG10(1) / ZD00E730E7 /
!     DATA LOG10(2) / ZC77800DC0 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O0000000000000000 /
!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O0007777777777777 /
!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /
!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /
!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O7770000000000000 /
!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O7777777777777777 /
!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /
!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /
!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA SMALL(1) / Z"3001800000000000" /
!     DATA SMALL(2) / Z"3001000000000000" /
!     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
!     DATA LARGE(2) / Z"4FFE000000000000" /
!     DATA RIGHT(1) / Z"3FD2800000000000" /
!     DATA RIGHT(2) / Z"3FD2000000000000" /
!     DATA DIVER(1) / Z"3FD3800000000000" /
!     DATA DIVER(2) / Z"3FD3000000000000" /
!     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
!     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA SMALL(1) / 00564000000000000000B /
!     DATA SMALL(2) / 00000000000000000000B /
!     DATA LARGE(1) / 37757777777777777777B /
!     DATA LARGE(2) / 37157777777777777777B /
!     DATA RIGHT(1) / 15624000000000000000B /
!     DATA RIGHT(2) / 00000000000000000000B /
!     DATA DIVER(1) / 15634000000000000000B /
!     DATA DIVER(2) / 00000000000000000000B /
!     DATA LOG10(1) / 17164642023241175717B /
!     DATA LOG10(2) / 16367571421742254654B /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn OR -pd8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CC0000000000000' /
!     DATA DMACH(4) / Z'3CD0000000000000' /
!     DATA DMACH(5) / Z'3FF34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'00010000000000000000000000000000' /
!     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3F900000000000000000000000000000' /
!     DATA DMACH(4) / Z'3F910000000000000000000000000000' /
!     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
!
!     MACHINE CONSTANTS FOR THE CRAY
!
!     DATA SMALL(1) / 201354000000000000000B /
!     DATA SMALL(2) / 000000000000000000000B /
!     DATA LARGE(1) / 577767777777777777777B /
!     DATA LARGE(2) / 000007777777777777774B /
!     DATA RIGHT(1) / 376434000000000000000B /
!     DATA RIGHT(2) / 000000000000000000000B /
!     DATA DIVER(1) / 376444000000000000000B /
!     DATA DIVER(2) / 000000000000000000000B /
!     DATA LOG10(1) / 377774642023241175717B /
!     DATA LOG10(2) / 000007571421742254654B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC DMACH(5)
!
!     DATA SMALL /    20K, 3*0 /
!     DATA LARGE / 77777K, 3*177777K /
!     DATA RIGHT / 31420K, 3*0 /
!     DATA DIVER / 32020K, 3*0 /
!     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA DMACH(1) / '0000000000000010'X /
!     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
!     DATA DMACH(3) / '0000000000003CC0'X /
!     DATA DMACH(4) / '0000000000003CD0'X /
!     DATA DMACH(5) / '79FF509F44133FF3'X /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FORMAT
!
      DATA DMACH(1) / Z'0010000000000000' /
      DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
      DATA DMACH(3) / Z'3CA0000000000000' /
      DATA DMACH(4) / Z'3CB0000000000000' /
      DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
!     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
!     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
!     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
!     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING D_FLOATING
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1), SMALL(2) /        128,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
!     DATA DIVER(1), DIVER(2) /       9472,           0 /
!     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
!
!     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING G_FLOATING
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1), SMALL(2) /         16,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
!     DATA DIVER(1), DIVER(2) /      15568,           0 /
!     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
!
!     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!     (ASSUMING Real*16 IS THE DEFAULT DOUBLE PRECISION)
!
!     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
!     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
!     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
!     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
!     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
!     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
!     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
!     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
!     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
!     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
!     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
!     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
!     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
!     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
!     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) /  40000B,       0 /
!     DATA SMALL(3), SMALL(4) /       0,       1 /
!     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
!     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
!     DATA RIGHT(3), RIGHT(4) /       0,    225B /
!     DATA DIVER(1), DIVER(2) /  40000B,       0 /
!     DATA DIVER(3), DIVER(4) /       0,    227B /
!     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
!     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
!     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
!     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
!     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
!     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
!     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
!     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
!
!     DATA SMALL(1) / 2.23D-308  /
!     DATA LARGE(1) / 1.79D+308  /
!     DATA RIGHT(1) / 1.11D-16   /
!     DATA DIVER(1) / 2.22D-16   /
!     DATA LOG10(1) / 0.301029995663981195D0 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
!     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
!     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
!     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
!     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
!     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
!     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
!     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
!     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /    8388608,           0 /
!     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
!     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
!     DATA DIVER(1), DIVER(2) /  620756992,           0 /
!     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
!
!     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
!     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
!     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
!     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
!     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /    128,      0 /
!     DATA SMALL(3), SMALL(4) /      0,      0 /
!     DATA LARGE(1), LARGE(2) /  32767,     -1 /
!     DATA LARGE(3), LARGE(4) /     -1,     -1 /
!     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
!     DATA RIGHT(3), RIGHT(4) /      0,      0 /
!     DATA DIVER(1), DIVER(2) /   9472,      0 /
!     DATA DIVER(3), DIVER(4) /      0,      0 /
!     DATA LOG10(1), LOG10(2) /  16282,   8346 /
!     DATA LOG10(3), LOG10(4) / -31493, -12296 /
!
!     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!     DATA SMALL(3), SMALL(4) / O000000, O000000 /
!     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!     DATA LARGE(3), LARGE(4) / O177777, O177777 /
!     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
!     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
!     DATA DIVER(1), DIVER(2) / O022400, O000000 /
!     DATA DIVER(3), DIVER(4) / O000000, O000000 /
!     DATA LOG10(1), LOG10(2) / O037632, O020232 /
!     DATA LOG10(3), LOG10(4) / O102373, O147770 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'00010000000000000000000000000000' /
!     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
!     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
!     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
!
!     MACHINE CONSTANTS FOR THE SUN 386i
!
!     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
!     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
!     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
!     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!
!     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
!     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
!     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
!     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
!     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
!
!***FIRST EXECUTABLE STATEMENT  D1MACH
!
  if ( I < 1 .OR. I > 5 ) then
    call XERMSG ('SLATEC', 'D1MACH', 'I OUT OF BOUNDS', 1, 2)
  end if

  D1MACH = DMACH(I)

  return
end


FUNCTION I1MACH (I)
!
!! I1MACH returns integer machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      INTEGER (I1MACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!        K = I1MACH(I)
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit.
!     I1MACH( 4) = the standard error message unit.
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0  <=  X(I)  <  A for I=0,...,S-1.
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0  <=  X(I)  <  B for I=1,...,T,
!                0  <  X(1), and EMIN  <=  E  <=  EMAX.
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   891012  Added VAX G-floating constants.  (WRB)
!   891012  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
!           (RWC)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added Convex -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
!           options.  (DWL, RWC and WRB).
!***END PROLOGUE  I1MACH
!
  integer i1mach
  INTEGER IMACH(16),OUTPUT
  SAVE IMACH
  EQUIVALENCE (IMACH(4),OUTPUT)
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT COMPILER
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        129 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1025 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA IMACH( 1) /          7 /
!     DATA IMACH( 2) /          2 /
!     DATA IMACH( 3) /          2 /
!     DATA IMACH( 4) /          2 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         33 /
!     DATA IMACH( 9) / Z1FFFFFFFF /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -256 /
!     DATA IMACH(13) /        255 /
!     DATA IMACH(14) /         60 /
!     DATA IMACH(15) /       -256 /
!     DATA IMACH(16) /        255 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         48 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /          8 /
!     DATA IMACH(11) /         13 /
!     DATA IMACH(12) /        -50 /
!     DATA IMACH(13) /         76 /
!     DATA IMACH(14) /         26 /
!     DATA IMACH(15) /        -50 /
!     DATA IMACH(16) /         76 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         48 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /          8 /
!     DATA IMACH(11) /         13 /
!     DATA IMACH(12) /        -50 /
!     DATA IMACH(13) /         76 /
!     DATA IMACH(14) /         26 /
!     DATA IMACH(15) /     -32754 /
!     DATA IMACH(16) /      32780 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -4095 /
!     DATA IMACH(13) /       4094 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -4095 /
!     DATA IMACH(16) /       4094 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /    6LOUTPUT/
!     DATA IMACH( 5) /         60 /
!     DATA IMACH( 6) /         10 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         48 /
!     DATA IMACH( 9) / 00007777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /       -929 /
!     DATA IMACH(13) /       1070 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /       -929 /
!     DATA IMACH(16) /       1069 /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / Z'7FFFFFFF' /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1023 /
!     DATA IMACH(13) /       1023 /
!     DATA IMACH(14) /        113 /
!     DATA IMACH(15) /     -16383 /
!     DATA IMACH(16) /      16383 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -pd8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1023 /
!     DATA IMACH(13) /       1023 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CRAY
!     USING THE 46 BIT INTEGER COMPILER OPTION
!
!     DATA IMACH( 1) /        100 /
!     DATA IMACH( 2) /        101 /
!     DATA IMACH( 3) /        102 /
!     DATA IMACH( 4) /        101 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         46 /
!     DATA IMACH( 9) / 1777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -8189 /
!     DATA IMACH(13) /       8190 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -8099 /
!     DATA IMACH(16) /       8190 /
!
!     MACHINE CONSTANTS FOR THE CRAY
!     USING THE 64 BIT INTEGER COMPILER OPTION
!
!     DATA IMACH( 1) /        100 /
!     DATA IMACH( 2) /        101 /
!     DATA IMACH( 3) /        102 /
!     DATA IMACH( 4) /        101 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 777777777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -8189 /
!     DATA IMACH(13) /       8190 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -8099 /
!     DATA IMACH(16) /       8190 /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!
!     DATA IMACH( 1) /         11 /
!     DATA IMACH( 2) /         12 /
!     DATA IMACH( 3) /          8 /
!     DATA IMACH( 4) /         10 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /         16 /
!     DATA IMACH(11) /          6 /
!     DATA IMACH(12) /        -64 /
!     DATA IMACH(13) /         63 /
!     DATA IMACH(14) /         14 /
!     DATA IMACH(15) /        -64 /
!     DATA IMACH(16) /         63 /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FLOAT
!
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          6 /
      DATA IMACH( 4) /          6 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -125 /
      DATA IMACH(13) /        128 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1021 /
      DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING D_FLOATING
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING G_FLOATING
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         32 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         24 /
!     DATA IMACH( 6) /          3 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         23 /
!     DATA IMACH( 9) /    8388607 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         38 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /         43 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         63 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          4 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         39 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          4 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         55 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          7 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         32 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1015 /
!     DATA IMACH(16) /       1017 /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) /  Z7FFFFFFF /
!     DATA IMACH(10) /         16 /
!     DATA IMACH(11) /          6 /
!     DATA IMACH(12) /        -64 /
!     DATA IMACH(13) /         63 /
!     DATA IMACH(14) /         14 /
!     DATA IMACH(15) /        -64 /
!     DATA IMACH(16) /         63 /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          5 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         54 /
!     DATA IMACH(15) /       -101 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          5 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         62 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGER ARITHMETIC.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGER ARITHMETIC.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1021 /
!     DATA IMACH(13) /       1024 /
!     DATA IMACH(14) /        113 /
!     DATA IMACH(15) /     -16381 /
!     DATA IMACH(16) /      16384 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          1 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         60 /
!     DATA IMACH(15) /      -1024 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA IMACH( 1) /          1 /
!     DATA IMACH( 2) /          1 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!***FIRST EXECUTABLE STATEMENT  I1MACH
!
  if ( I < 1 .OR. I > 16 ) then
    WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
    STOP
  end if

  I1MACH = IMACH(I)

  return
end

FUNCTION R1MACH (I)
!
!! R1MACH returns floating point machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   R1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = R1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   R1MACH(3) = B**(-T), the smallest relative spacing.
!   R1MACH(4) = B**(1-T), the largest relative spacing.
!   R1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0  <=  X(I)  <  B for I=1,...,T, 0  <  X(1), and
!   EMIN  <=  E  <=  EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!***END PROLOGUE  R1MACH
!
  real r1mach
  INTEGER SMALL(2)
  INTEGER LARGE(2)
  INTEGER RIGHT(2)
  INTEGER DIVER(2)
  INTEGER LOG10(2)
!
  REAL RMACH(5)
  SAVE RMACH
!
  EQUIVALENCE (RMACH(1),SMALL(1))
  EQUIVALENCE (RMACH(2),LARGE(1))
  EQUIVALENCE (RMACH(3),RIGHT(1))
  EQUIVALENCE (RMACH(4),DIVER(1))
  EQUIVALENCE (RMACH(5),LOG10(1))
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7F7FFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7EFFFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA SMALL(1) / 16#00800000 /
!     DATA LARGE(1) / 16#7FFFFFFF /
!     DATA RIGHT(1) / 16#33800000 /
!     DATA DIVER(1) / 16#34000000 /
!     DATA LOG10(1) / 16#3E9A209B /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA RMACH(1) / Z400800000 /
!     DATA RMACH(2) / Z5FFFFFFFF /
!     DATA RMACH(3) / Z4E9800000 /
!     DATA RMACH(4) / Z4EA800000 /
!     DATA RMACH(5) / Z500E730E8 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
!
!     DATA RMACH(1) / O1771000000000000 /
!     DATA RMACH(2) / O0777777777777777 /
!     DATA RMACH(3) / O1311000000000000 /
!     DATA RMACH(4) / O1301000000000000 /
!     DATA RMACH(5) / O1157163034761675 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA RMACH(1) / Z"3001800000000000" /
!     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
!     DATA RMACH(3) / Z"3FD2800000000000" /
!     DATA RMACH(4) / Z"3FD3800000000000" /
!     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA RMACH(1) / 00564000000000000000B /
!     DATA RMACH(2) / 37767777777777777776B /
!     DATA RMACH(3) / 16414000000000000000B /
!     DATA RMACH(4) / 16424000000000000000B /
!     DATA RMACH(5) / 17164642023241175720B /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7F7FFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn COMPILER OPTION
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7FFFFFFF' /
!     DATA RMACH(3) / Z'34800000' /
!     DATA RMACH(4) / Z'35000000' /
!     DATA RMACH(5) / Z'3F9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 OR -pd8 COMPILER OPTION
!
!     DATA RMACH(1) / Z'0010000000000000' /
!     DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' /
!     DATA RMACH(3) / Z'3CC0000000000000' /
!     DATA RMACH(4) / Z'3CD0000000000000' /
!     DATA RMACH(5) / Z'3FF34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CRAY
!
!     DATA RMACH(1) / 200034000000000000000B /
!     DATA RMACH(2) / 577767777777777777776B /
!     DATA RMACH(3) / 377224000000000000000B /
!     DATA RMACH(4) / 377234000000000000000B /
!     DATA RMACH(5) / 377774642023241175720B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC RMACH(5)
!
!     DATA SMALL /    20K,       0 /
!     DATA LARGE / 77777K, 177777K /
!     DATA RIGHT / 35420K,       0 /
!     DATA DIVER / 36020K,       0 /
!     DATA LOG10 / 40423K,  42023K /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA RMACH(1) / '00000080'X /
!     DATA RMACH(2) / 'FFFF7FFF'X /
!     DATA RMACH(3) / '00003480'X /
!     DATA RMACH(4) / '00003500'X /
!     DATA RMACH(5) / '209B3F9A'X /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FLOAT
!
      DATA RMACH(1) / Z'00800000' /
      DATA RMACH(2) / Z'7F7FFFFF' /
      DATA RMACH(3) / Z'33800000' /
      DATA RMACH(4) / Z'34000000' /
      DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1) /       128 /
!     DATA LARGE(1) /    -32769 /
!     DATA RIGHT(1) /     13440 /
!     DATA DIVER(1) /     13568 /
!     DATA LOG10(1) / 547045274 /
!
!     DATA SMALL(1) / Z00000080 /
!     DATA LARGE(1) / ZFFFF7FFF /
!     DATA RIGHT(1) / Z00003480 /
!     DATA DIVER(1) / Z00003500 /
!     DATA LOG10(1) / Z209B3F9A /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!     (ASSUMING REAL*4 IS THE DEFAULT REAL)
!
!     DATA SMALL(1) / '00800000'X /
!     DATA LARGE(1) / '7F7FFFFF'X /
!     DATA RIGHT(1) / '33800000'X /
!     DATA DIVER(1) / '34000000'X /
!     DATA LOG10(1) / '3E9A209B'X /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
!     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
!     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
!     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA RMACH(1) / O402400000000 /
!     DATA RMACH(2) / O376777777777 /
!     DATA RMACH(3) / O714400000000 /
!     DATA RMACH(4) / O716400000000 /
!     DATA RMACH(5) / O776464202324 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA SMALL(1) / 00004000000B /
!     DATA LARGE(1) / 17677777777B /
!     DATA RIGHT(1) / 06340000000B /
!     DATA DIVER(1) / 06400000000B /
!     DATA LOG10(1) / 07646420233B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA RMACH(1) / Z00100000 /
!     DATA RMACH(2) / Z7FFFFFFF /
!     DATA RMACH(3) / Z3B100000 /
!     DATA RMACH(4) / Z3C100000 /
!     DATA RMACH(5) / Z41134413 /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!
!     DATA SMALL(1) / 1.18E-38      /
!     DATA LARGE(1) / 3.40E+38      /
!     DATA RIGHT(1) / 0.595E-07     /
!     DATA DIVER(1) / 1.19E-07      /
!     DATA LOG10(1) / 0.30102999566 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
!
!     DATA RMACH(1) / "000400000000 /
!     DATA RMACH(2) / "377777777777 /
!     DATA RMACH(3) / "146400000000 /
!     DATA RMACH(4) / "147400000000 /
!     DATA RMACH(5) / "177464202324 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1) /    8388608 /
!     DATA LARGE(1) / 2147483647 /
!     DATA RIGHT(1) /  880803840 /
!     DATA DIVER(1) /  889192448 /
!     DATA LOG10(1) / 1067065499 /
!
!     DATA RMACH(1) / O00040000000 /
!     DATA RMACH(2) / O17777777777 /
!     DATA RMACH(3) / O06440000000 /
!     DATA RMACH(4) / O06500000000 /
!     DATA RMACH(5) / O07746420233 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /   128,     0 /
!     DATA LARGE(1), LARGE(2) / 32767,    -1 /
!     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
!     DATA DIVER(1), DIVER(2) / 13568,     0 /
!     DATA LOG10(1), LOG10(2) / 16282,  8347 /
!
!     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
!     DATA DIVER(1), DIVER(2) / O032400, O000000 /
!     DATA LOG10(1), LOG10(2) / O037632, O020233 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA RMACH(1) / Z'0010000000000000' /
!     DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA RMACH(3) / Z'3CA0000000000000' /
!     DATA RMACH(4) / Z'3CB0000000000000' /
!     DATA RMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
!
!     DATA RMACH(1) / O000400000000 /
!     DATA RMACH(2) / O377777777777 /
!     DATA RMACH(3) / O146400000000 /
!     DATA RMACH(4) / O147400000000 /
!     DATA RMACH(5) / O177464202324 /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA SMALL(1), SMALL(2) /     0,    256/
!     DATA LARGE(1), LARGE(2) /    -1,   -129/
!     DATA RIGHT(1), RIGHT(2) /     0,  26880/
!     DATA DIVER(1), DIVER(2) /     0,  27136/
!     DATA LOG10(1), LOG10(2) /  8347,  32538/
!
!***FIRST EXECUTABLE STATEMENT  R1MACH
!
  if ( I < 1 .OR. I > 5 ) then
    call XERMSG ('SLATEC', 'R1MACH', 'I OUT OF BOUNDS', 1, 2)
  end if

  R1MACH = RMACH(I)

  return
end

!subroutine avint ( ftab, xtab, ntab, a, b, result )
!!
!!***********************************************************************
!!
!!! AVINT estimates the integral of unevenly spaced data.
!!
!!
!!  Discussion:
!!
!!    The method uses overlapping parabolas and smoothing.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!    P E Hennion,
!!    Algorithm 77,
!!    Interpolation, Differentiation and Integration,
!!    Communications of the Association for Computing Machinery,
!!    Volume 5, page 96, 1962.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real FTAB(NTAB), the function values,
!!    FTAB(I) = F(XTAB(I)).
!!
!!    Input, real XTAB(NTAB), the abscissas at which the
!!    function values are given.  The XTAB's must be distinct
!!    and in ascending order.
!!
!!    Input, integer NTAB, the number of entries in FTAB and
!!    XTAB.  NTAB must be at least 3.
!!
!!    Input, real A, the lower limit of integration.  A should
!!    be, but need not be, near one endpoint of the interval
!!    (X(1), X(NTAB)).
!!
!!    Input, real B, the upper limit of integration.  B should
!!    be, but need not be, near one endpoint of the interval
!!    (X(1), X(NTAB)).
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!  implicit none
!!
!  integer ntab
!!
!  real a
!  real atemp
!  real b
!  real btemp
!  real ca
!  real cb
!  real cc
!  real ctemp
!  real ftab(ntab)
!  integer i
!  integer ihi
!  integer ilo
!  integer ind
!  real result
!  real sum1
!  real syl
!  real term1
!  real term2
!  real term3
!  real x1
!  real x2
!  real x3
!  real xtab(ntab)
!!
!  if ( ntab < 3 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'AVINT - Fatal error!'
!    write ( *, '(a,i6)' ) '  NTAB is less than 3.  NTAB = ', ntab
!    stop
!  end if
! 
!  do i = 2, ntab
! 
!    if ( xtab(i) <= xtab(i-1) ) then
!      write ( *, '(a)' ) ' '
!      write ( *, '(a)' ) 'AVINT - Fatal error!'
!      write ( *, '(a)' ) '  XTAB(I) is not greater than XTAB(I-1).'
!      write ( *, '(a,i6)' ) '  Here, I = ', I
!      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
!      write ( *, '(a,g14.6)' ) '  XTAB(I) =   ', xtab(i)
!      stop
!    end if
! 
!  end do
! 
!  result = 0.0E+00
! 
!  if ( a == b ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'AVINT - Warning!'
!    write ( *, '(a)' ) '  A = B, integral=0.'
!    return
!  end if
!!
!!  If A > B, temporarily switch A and B, and store sign.
!!
!  if ( a > b ) then
!    syl = b
!    b = a
!    a = syl
!    ind = -1
!  else
!    syl = a
!    ind = 1
!  end if
!!
!!  Bracket A and B between XTAB(ILO) and XTAB(IHI).
!!
!  ilo = 1
!  ihi = ntab
!
!  do i = 1, ntab
!    if ( xtab(i) >= a ) then
!      exit
!    end if
!    ilo = ilo + 1
!  end do
!
!  ilo = max ( 2, ilo )
!  ilo = min ( ilo, ntab-1 )
!
!  do i = 1, ntab
!    if ( b >= xtab(i) ) then
!      exit
!    end if
!    ihi = ihi - 1
!  end do
!  
!  ihi = min ( ihi, ntab-1 )
!  ihi = max ( ilo, ihi-1 )
!!
!!  Carry out approximate integration from XTAB(ILO) to XTAB(IHI).
!!
!  sum1 = 0.0E+00
! 
!  do i = ilo, ihi
! 
!    x1 = xtab(i-1)
!    x2 = xtab(i)
!    x3 = xtab(i+1)
! 
!    term1 = ftab(i-1) / ((x1-x2)*(x1-x3))
!    term2 = ftab(i) / ((x2-x1)*(x2-x3))
!    term3 = ftab(i+1) / ((x3-x1)*(x3-x2))
! 
!    atemp = term1 + term2 + term3
!    btemp = -(x2+x3)*term1-(x1+x3)*term2-(x1+x2)*term3
!    ctemp = x2*x3*term1+x1*x3*term2+x1*x2*term3
! 
!    if ( i <= ilo ) then
!      ca = atemp
!      cb = btemp
!      cc = ctemp
!    else
!      ca = 0.5E+00 * ( atemp + ca )
!      cb = 0.5E+00 * ( btemp + cb )
!      cc = 0.5E+00 * ( ctemp + cc )
!    end if
! 
!    sum1 = sum1 &
!          + ca * ( x2**3 - syl**3 ) / 3.0E+00 &
!          + cb * 0.5E+00 * ( x2**2 - syl**2 ) &
!          + cc * ( x2 - syl )
! 
!    ca = atemp
!    cb = btemp
!    cc = ctemp
! 
!    syl = x2
! 
!  end do
! 
!  result = sum1 &
!        + ca * ( b**3 - syl**3 ) / 3.0E+00 &
!        + cb * 0.5E+00 * ( b**2 - syl**2 ) &
!        + cc * ( b - syl )
!!
!!  Restore original values of A and B, reverse sign of integral
!!  because of earlier switch.
!!
!  if ( ind /= 1 ) then
!    ind = 1
!    syl = b
!    b = a
!    a = syl
!    result = -result
!  end if
! 
!  return
!end
!subroutine cadre ( func, a, b, abserr, relerr, error, result, ind )
!!
!!***********************************************************************
!!
!!! CADRE estimates the integral of F(X) from A to B.
!!
!!
!!  Discussion:
!!
!!    CADRE is the Cautious Adaptive Romberg Extrapolator.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!    Carl DeBoor and J R Rice,
!!    CADRE: An algorithm for numerical quadrature,
!!    Mathematic Software, pages 417-449,
!!    Academic Press, New York, 1971.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real, external FUNC, the name of the function to be integrated.
!!    The user must declare the name an external parameter in the calling
!!    program, write a function routine of the form FUNCTION FUNC(X) which
!!    evaluates the function at X, and pass the name of the function
!!    in FUNC.
!!
!!    Input, real A, the lower limit of integration.
!!
!!    Input, real B, the upper limit of integration.
!!
!!    Input, real ABSERR, the absolute error tolerance.
!!
!!    Input, real RELERR, the relative error tolerance.
!!
!!    Output, real ERROR, an estimate of the absolute error.
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!!    Output, integer IND, reliability indicator.
!!    If IND <= 2, RESULT is very reliable.  Higher values of
!!    IND indicate less reliable values of RESULT.
!!
!  implicit none
!!
!  integer, parameter :: mxstge = 30
!  integer, parameter :: maxtbl = 10
!  integer, parameter :: maxts = 2049
!!
!  real a
!  real abserr
!  real ait(maxtbl)
!  logical aitken
!  real, parameter :: aitlow = 1.1E+00
!  real, parameter :: aittol = 0.1E+00
!  real astep
!  real b
!  real beg
!  real begin(mxstge)
!  real bma
!  real curest
!  real dif(maxtbl)
!  real diff
!  real end
!  real ergoal
!  real erra
!  real errer
!  real error
!  real errr
!  real est(mxstge)
!  real fbeg
!  real fbeg2
!  real fend
!  real fextm1
!  real fextrp
!  real finis(mxstge)
!  real fn
!  real fnsize
!  real, external :: func
!  logical h2conv
!  real h2next
!  real h2tfex
!  real, parameter :: h2tol = 0.15E+00
!  real hovn
!  integer i
!  integer ibeg
!  integer ibegs(mxstge)
!  integer iend
!  integer ii
!  integer iii
!  integer ind
!  integer istage
!  integer istep
!  integer istep2
!  integer it
!  integer l
!  integer lm1
!  integer n
!  integer n2
!  integer nnleft
!  real prever
!  real r(maxtbl)
!  logical reglar
!  logical reglsv(mxstge)
!  real relerr
!  real result
!  logical right
!  real rn(4)
!  real rnderr
!  real sing
!  real singnx
!  real slope
!  real stage
!  real step
!  real stepmn
!  real sum1
!  real sumabs
!  real t(maxtbl,maxtbl)
!  real tabs
!  real tabtlm
!  real, parameter :: tljump = 0.01E+00
!  real ts(2049)
!  real vint
!!
!  if ( a == b ) then
!    result = 0.0E+00
!    return
!  end if
! 
!  begin(1:mxstge) = 0.0E+00
!  est(1:mxstge) = 0.0E+00
!  finis(1:mxstge) = 0.0E+00
!  ibegs(1:mxstge) = 0
!  reglsv(1:mxstge) = .false.
! 
!  vint = 0.0E+00
! 
!  rn(1:4) = (/ 0.7142005E+00, 0.3466282E+00, 0.8437510E+00, 0.1263305E+00 /)
! 
!  rnderr = epsilon ( rnderr )
!  result = 0.0E+00
!  error = 0.0E+00
!  ind = 1
!  bma = abs ( b - a )
!  errr = min ( 0.1E+00, max ( abs ( relerr ), 10.0E+00*rnderr) )
!  erra = abs ( abserr )
!  stepmn = max ( bma / 2**mxstge, max ( bma, abs ( a ), abs ( b ) ) * rnderr )
!  stage = 0.5
!  istage = 1
!  curest = 0.0E+00
!  fnsize = 0.0E+00
!  prever = 0.0E+00
!  reglar = .false.
!  beg = a
!  fbeg = func(beg) / 2.0E+00
!  ts(1) = fbeg
!  ibeg = 1
!  end = b
!  fend = func(end) / 2.0E+00
!  ts(2) = fend
!  iend = 2
! 
!10 continue
! 
!  right = .false.
! 
!20 continue
!
!  step = end - beg
!  astep = abs ( step )
! 
!  if ( astep < stepmn ) then
!    ind = 5
!    result = curest + vint
!    return
!  end if
! 
!  t(1,1) = fbeg+fend
!  tabs = abs ( fbeg ) + abs ( fend )
!  l = 1
!  n = 1
!  h2conv = .false.
!  aitken = .false.
!  go to 40
! 
!30 continue
! 
!40 continue
! 
!  lm1 = l
!  l = l+1
!  n2 = n*2
!  fn = n2
!  istep = (iend-ibeg)/n
!
!  if ( istep > 1 ) then
!    go to 60
!  end if
!
!  ii = iend
!  iend = iend + n
!
!  if ( iend > maxts ) then
!    go to 440
!  end if
!
!  hovn = step / fn
! 
!  iii = iend
!  do i = 1, n2, 2
!    ts(iii) = ts(ii)
!    ts(iii-1) = func(end-i*hovn)
!    iii = iii-2
!    ii = ii-1
!  end do
! 
!  istep = 2
! 
!60 continue
! 
!  istep2 = ibeg+istep/2
! 
!  sum1 = 0.0E+00
!  sumabs = 0.0E+00
!  do i = istep2, iend, istep
!    sum1 = sum1 + ts(i)
!    sumabs = sumabs + abs ( ts(i) )
!  end do
! 
!  t(l,1) = t(l-1,1) / 2.0E+00 + sum1 / fn
!  tabs = tabs / 2.0E+00 + sumabs / fn
! 
!  n = n2
!  it = 1
!  vint = step * t(l,1)
!  tabtlm = tabs * rnderr
!  fnsize = max ( fnsize, abs ( t(l,1) ) )
!  ergoal = max ( astep * rnderr * fnsize, &
!    stage * max ( erra , errr * abs ( curest+vint ) ) )
!  fextrp = 1.0E+00
!  do i = 1, lm1
!    fextrp = fextrp * 4.0E+00
!    t(i,l) = t(l,i) - t(l-1,i)
!    t(l,i+1) = t(l,i) + t(i,l) / ( fextrp - 1.0 )
!  end do
! 
!  errer = astep * abs ( t(1,l) )
!  if ( l > 2 ) go to 90
!  if ( abs ( t(1,2) ) <= tabtlm ) go to 290
!  go to 40
! 
!90 continue
! 
!  do i = 2, lm1
!
!    if ( abs ( t(i-1,l) ) > tabtlm ) then
!      diff = t(i-1,lm1) / t(i-1,l)
!    else
!      diff = 0.0E+00
!    end if
!
!    t(i-1,lm1) = diff
!
!  end do
! 
!  if ( abs ( 4.0 - t(1,lm1) ) <= h2tol ) go to 130
!  if ( t(1,lm1) == 0.0 ) go to 120
!  if ( abs ( 2.0 - abs ( t(1,lm1) ) ) < tljump ) go to 280
!  if (l==3) go to 30
!  h2conv = .false.
!  if ( abs ( ( t(1,lm1) - t(1,l-2) ) / t(1,lm1) ) <= aittol ) go to 160
! 
!  if ( .not. reglar .and. l == 4 ) go to 30
! 
!120 continue
! 
!  if ( errer <= ergoal ) go to 310
!  go to 380
!
!130 continue
!
!  if ( .not. h2conv ) then
!    aitken = .false.
!    h2conv = .true.
!  end if
!
!140 continue
!
!  fextrp = 4.0E+00
!
!150 continue
!
!  it = it+1
!  vint = step * t(l,it)
!  errer = abs ( step / (fextrp-1.0) * t(it-1,l))
!  if ( errer <= ergoal ) go to 340
!  if ( it == lm1 ) go to 270
!  if ( t(it,lm1) == 0.0 ) go to 150
!  if ( t(it,lm1) <= fextrp ) go to 270
!
!  if ( abs ( t(it,lm1) / 4.0 - fextrp ) / fextrp < aittol ) then
!    fextrp = fextrp*4.0E+00
!  end if
!
!  go to 150
! 
!160 continue
!
!  if ( t(1,lm1) < aitlow ) then
!    go to 380
!  end if
! 
!  if ( .not. aitken ) then
!    h2conv = .false.
!    aitken = .true.
!  end if
! 
!170 continue
!
!  fextrp = t(l-2,lm1)
!  if ( fextrp > 4.5 ) go to 140
!  if ( fextrp < aitlow ) go to 380
!
!  if ( abs ( fextrp - t(l-3,lm1) ) / t(1,lm1) > h2tol ) then
!    go to 380
!  end if
!
!  sing = fextrp
!  fextm1 = fextrp - 1.0E+00
!
!  ait(1) = 0.0E+00
!  do i = 2, l
!    ait(i) = t(i,1) + (t(i,1)-t(i-1,1)) / fextm1
!    r(i) = t(1,i-1)
!    dif(i) = ait(i) - ait(i-1)
!  end do
!
!  it = 2
!
!190 continue
!
!  vint = step*ait(l)
!
!200 continue
!
!  errer = errer / fextm1
! 
!  if ( errer <= ergoal ) then
!    ind = max ( ind, 2 )
!    go to 340
!  end if
! 
!210 continue
!
!  it = it+1
!  if ( it == lm1 ) go to 270
!
!  if ( it <= 3 ) then
!    h2next = 4.0E+00
!    singnx = 2.0E+00 * sing
!  end if
!
!  if ( h2next < singnx ) go to 230
!  fextrp = singnx
!  singnx = 2.0E+00 * singnx
!  go to 240
!
!230 continue
!
!  fextrp = h2next
!  h2next = 4.0E+00 * h2next
!
!240 continue
! 
!  do i = it, lm1
!    if ( abs ( dif(i+1) ) > tabtlm ) then
!      r(i+1) = dif(i) / dif(i+1)
!    else
!      r(i+1) = 0.0E+00
!    end if
!  end do
! 
!  h2tfex = -h2tol*fextrp
!  if ( r(l) - fextrp < h2tfex ) go to 270
!  if ( r(l-1) - fextrp < h2tfex ) go to 270
!  errer = astep * abs ( dif(l) )
!  fextm1 = fextrp - 1.0E+00
!  do i = it, l
!    ait(i) = ait(i)+dif(i) / fextm1
!    dif(i) = ait(i)-ait(i-1)
!  end do
! 
!  go to 190
! 
!270 continue
!
!  fextrp = max(prever/errer,aitlow)
!  prever = errer
!  if (l<5) go to 40
!  if (l-it>2.and.istage<mxstge) go to 370
!  if (errer/fextrp**(maxtbl-l)<ergoal) go to 40
!  go to 370
! 
!280 continue
!
!  if ( errer > ergoal ) go to 370
!  diff = abs ( t(1,l) ) * 2.0E+00 * fn
!  go to 340
! 
!290 continue
!
!  slope = (fend-fbeg) * 2.0E+00
!  fbeg2 = fbeg * 2.0E+00
! 
!  do i = 1, 4
!    diff = abs ( func ( beg + rn(i) * step ) - fbeg2 - rn(i) * slope )
!    if ( diff > tabtlm ) go to 330
!  end do
! 
!  go to 340
! 
!310 continue
!
!  slope = (fend-fbeg)*2.0E+00
!  fbeg2 = fbeg*2.0E+00
!  i = 1
! 
!320 continue
!
!  diff = abs ( func(beg+rn(i)*step) - fbeg2 - rn(i) * slope )
! 
!330 continue
!
!  errer = max ( errer, astep * diff )
!  if (errer > ergoal) go to 380
!  i = i+1
!  if ( i <= 4 ) go to 320
!  ind = 3
! 
!340 continue
!
!  result = result+vint
!  error = error+errer
! 
!350 continue
!
!  if (right) go to 360
!  istage = istage-1
!  if (istage==0) return
!  reglar = reglsv(istage)
!  beg = begin(istage)
!  end = finis(istage)
!  curest = curest-est(istage+1)+vint
!  iend = ibeg-1
!  fend = ts(iend)
!  ibeg = ibegs(istage)
!  go to 400
! 
!360 continue
!
!  curest = curest+vint
!  stage = stage*2.0E+00
!  iend = ibeg
!  ibeg = ibegs(istage)
!  end = beg
!  beg = begin(istage)
!  fend = fbeg
!  fbeg = ts(ibeg)
!  go to 10
! 
!370 continue
!
!  reglar = .true.
! 
!380 continue
! 
!  if ( istage == mxstge ) then
!    ind = 5
!    result = curest+vint
!    return
!  end if
! 
!390 continue
!
!  if (right) go to 410
!  reglsv(istage+1) = reglar
!  begin(istage) = beg
!  ibegs(istage) = ibeg
!  stage = stage/2.0E+00
!
!400 continue
!
!  right = .true.
!  beg = (beg+end)/2.0E+00
!  ibeg = (ibeg+iend)/2
!  ts(ibeg) = ts(ibeg) / 2.0E+00
!  fbeg = ts(ibeg)
!  go to 20
!
!410 continue
!
!  nnleft = ibeg-ibegs(istage)
!  if (end+nnleft>=maxts) go to 440
!  iii = ibegs(istage)
!  ii = iend
!  do i = iii, ibeg
!    ii = ii+1
!    ts(ii) = ts(i)
!  end do
! 
!  do i = ibeg, ii
!    ts(iii) = ts(i)
!    iii = iii+1
!  end do
! 
!  iend = iend+1
!  ibeg = iend-nnleft
!  fend = fbeg
!  fbeg = ts(ibeg)
!  finis(istage) = end
!  end = beg
!  beg = begin(istage)
!  begin(istage) = end
!  reglsv(istage) = reglar
!  istage = istage+1
!  reglar = reglsv(istage)
!  est(istage) = vint
!  curest = curest+est(istage)
!  go to 10
!
!440 continue
!
!  ind = 4
!
!460 continue
!
!  result = curest+vint
!
!  return
!end
!subroutine chinsp ( func, a, b, epsin, epsout, result )
!!
!!***********************************************************************
!!
!!! CHINSP estimates an integral using a modified Clenshaw-Curtis scheme.
!!
!!
!!  Discussion:
!!
!!    The integral is approximated by Chebyshev polyonomials over each
!!    subinterval.  These are integrated to give the approximate integral.
!!    If the error estimate is unsatisfactory, the integration is repeated
!!    with smaller intervals.
!!
!!    The internal parameter NUPPER is currently set to 9,
!!    corresponding to 1024 subintervals for the unfolded integral,
!!    and 1025 function evaluations.  This parameter may be changed
!!    if necessary.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!    T Havie
!!    BIT 9 (1969), pages 338-350.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real, external FUNC, the name of the function to be
!!    integrated.  The user must declare the name an external
!!    parameter in the calling program, pass the name of the
!!    function in FUNC, and write a function of the form
!!
!!      FUNCTION FUNC(X)
!!
!!    which evaluates the function at the point X.
!!
!!    Input, real A, the lower limit of integration.
!!
!!    Input, real B, the upper limit of integration.
!!
!!    Input, real EPSIN, the relative error tolerance.
!!
!!    Output, real EPSOUT, estimated integration error.
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!  implicit none
!!
!  integer, parameter :: nupper = 9
!!
!  real a
!  real a0
!  real a1
!  real a2
!  real acof(257)
!  real alf
!  real alfnj
!  real alfno
!  real b
!  real bcof(257)
!  real bet
!  real betnj
!  real betno
!  real bounds
!  real ccof(513)
!  real cof
!  real cofmax
!  real const1
!  real const2
!  real deln
!  real deltan
!  real epsin
!  real epsout
!  real error
!  real etank
!  real, external :: func
!  real gamman
!  real hnstep
!  integer i
!  integer index
!  integer j
!  integer k
!  integer ksign
!  integer n
!  integer ncof
!  integer nhalf
!  integer nn
!  real, parameter :: one = 1.0E+00
!  real r1
!  real r2
!  real result
!  real rk
!  real rn
!  real rnderr
!  real rounde
!  real tend
!  real tnew
!  real triarg
!  real umid
!  real wmean
!  real xmin
!  real xplus
!  real xsink
!!
!  if ( a == b ) then
!    result = 0.0E+00
!    return
!  end if
!!
!!  ROUNDE = RNDERR*(R1+R2*N), where R1, R2 are two empirical constants.
!!
!!  Set coefficients in formula for accumulated roundoff error.
!!  N is the current number of function values used.
!!
!  rnderr = epsilon ( 1.0E+00 )
! 
!  r1 = 1.0E+00
!  r2 = 2.0E+00
!  error = epsin
!!
!!  Integration interval parameters.
!!
!  alf = 0.5E+00 * ( b - a )
!  bet = 0.5E+00 * ( b + a )
!!
!!  Parameters for trigonometric recurrence relations.
!!
!  triarg = atan ( 1.0E+00 )
!  alfno = -1.0E+00
!!
!!  Parameters for integration stepsize and loops.
!!
!  rn = 2.0E+00
!  n = 2
!  nhalf = 1
!  hnstep = 1.0E+00
!!
!!  Initial calculation for the end-point approximation.
!!
!  const1 = 0.5E+00 * ( func(a) + func(b) )
!  const2 = func(bet)
!  acof(1) = 0.5E+00 * (const1+const2)
!  acof(2) = 0.5E+00 * (const1-const2)
!  bcof(2) = acof(2)
!  tend = 2.0E+00 * ( acof(1) - acof(2) / 3.0E+00 )
!!
!!  Start actual calculations.
!!
!  do i = 1, nupper
!!
!!  Compute function values.
!!
!    const1 = -sin(triarg)
!    const2 = 0.5E+00 * alfno / const1
!    alfno = const1
!    betno = const2
!    gamman = 1.0E+00 - 2.0E+00 * alfno**2
!    deltan = -2.0E+00 * alfno * betno
!    bcof(1) = 0.0E+00
! 
!    do j = 1, nhalf
!      alfnj = gamman * const1 + deltan*const2
!      betnj = gamman * const2 - deltan*const1
!      xplus = alf * alfnj+bet
!      xmin = -alf * alfnj+bet
!      ccof(j) = func(xplus) + func(xmin)
!      bcof(1) = bcof(1) + ccof(j)
!      const1 = alfnj
!      const2 = betnj
!    end do
! 
!    bcof(1) = 0.5E+00 * hnstep * bcof(1)
!!
!!  Calculation of first B-coefficient finished compute the higher
!!  coefficients if NHALF greater than one.
!!
!    if ( nhalf <= 1 ) go to 60
!    const1 = one
!    const2 = 0.0E+00
!    ncof = nhalf-1
!    ksign = -1
! 
!    do k = 1, ncof
!!
!!  Compute trigonometric sum for B-coefficient.
!!
!      etank = gamman * const1 - deltan*const2
!      xsink = gamman * const2 + deltan*const1
!      cof = 2.0E+00 * ( 2.0E+00 * etank**2 - 1.0E+00 )
!      a2 = 0.0E+00
!      a1 = 0.0E+00
!      a0 = ccof(nhalf)
! 
!      do j = 1, ncof
!        a2 = a1
!        a1 = a0
!        index = nhalf-j
!        a0 = ccof(index) + cof * a1 - a2
!      end do
! 
!      bcof(k+1) = hnstep * (a0-a1) * etank
!      bcof(k+1) = ksign * bcof(k+1)
!      ksign = -ksign
!      const1 = etank
!      const2 = xsink
! 
!    end do
!!
!!  Calculation of B-coefficients finished.
!!
!!  Compute new modified mid-point approximation when the interval
!!  of integration is divided in N equal sub intervals.
!!
!60  continue
! 
!    umid = 0.0E+00
!    rk = rn
!    nn = nhalf+1
!    do k = 1, nn
!      index = nn+1-k
!      umid = umid+bcof(index)/(rk**2-one)
!      rk = rk-2.0E+00
!    end do
! 
!    umid = -2.0E+00 * umid
!!
!!  Compute new C-coefficients for end-point approximation and largest
!!  absolute value of coefficients.
!!
!    nn = n+2
!    cofmax = 0.0E+00
! 
!    do j = 1, nhalf
!      index = nn-j
!      ccof(j) = 0.5E+00 * (acof(j)+bcof(j))
!      ccof(index) = 0.5E+00 * (acof(j)-bcof(j))
!      const1 = abs ( ccof(j) )
!      cofmax = max ( cofmax, const1 )
!      const1 = abs ( ccof(index) )
!      cofmax = max ( cofmax, const1 )
!    end do
! 
!    ccof(nhalf+1) = acof(nhalf+1)
!!
!!  Calculation of new coefficients finished.
!!
!!  Compute new end-point approximation when the interval of
!!  integration is divided in 2N equal sub intervals.
!!
!    wmean = 0.5E+00 * (tend+umid)
!    bounds = 0.5E+00 * (tend-umid)
!    deln = 0.0E+00
!    rk = 2.0E+00 * rn
!    do j = 1, nhalf
!      index = n+2-j
!      deln = deln+ccof(index) / (rk**2-one)
!      rk = rk-2.0E+00
!    end do
! 
!    deln = -2.0E+00 * deln
!    tnew = wmean+deln
!    epsout = abs ( bounds / tnew )
!
!    if ( cofmax < rnderr ) then
!      go to 160
!    end if
!
!    rounde = rnderr*(r1+r2*rn)
!    if ( epsout < rounde ) epsout = rounde
!    if ( error < rounde ) error = rounde
!    if ( epsout > error ) go to 160
!!
!!  Required accuracy obtained or the maximum number of function
!!  values used without obtaining the required accuracy.
!!
!120 continue
! 
!    n = 2*n+1
!    tend = alf*(tend+deln)
!    umid = alf*(umid+deln)
!    deln = alf*deln
!    result = alf*tnew
!    return
!!
!!  If I = NUPPER then the required accuracy is not obtained.
!!
!160 continue
! 
!    if ( i == nupper ) go to 120
! 
!    acof(1:n) = ccof(1:n)
!    acof(n+1) = ccof(n+1)
!    bcof(n+1) = ccof(n+1)
!    tend = tnew
!    nhalf = n
!    n = 2 * n
!    rn = 2.0E+00 * rn
!    hnstep = 0.5E+00 * hnstep
!    triarg = 0.5E+00 * triarg
! 
!  end do
! 
!  return
!end
!subroutine class ( kind, n, alpha, beta, b, a, muzero )
!!
!!***********************************************************************
!!
!!! CLASS sets recurrence coeeficients for various orthogonal polynomials.
!!
!!
!!  Discussion:
!!
!!    CLASS supplies the coefficients A(J), B(J) of the recurrence relation
!!
!!      B(J)*P(J) (X) = (X-A(J))*P(J-1)(X) - B(J-1)*P(J-2)(X)
!!
!!    for the various classical (normalized) orthogonal polynomials,
!!    and the zero-th moment
!!
!!      MUZERO = Integral W(X) DX
!!
!!    of the given polynomial's weight function W(X).  Since the
!!    polynomials are orthonormalized, the tridiagonal matrix is
!!    guaranteed to be symmetric.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real ALPHA, BETA, parameters needed for Laguerre and Jacobi
!!    polynomials.
!!
!!    Input, integer KIND, specifies which polynomial is to be handled:
!!
!!    1: Legendre polynomials P(X) on (-1, +1),
!!    W(X) = 1.
!!
!!    2: Chebyshev polynomials of the first kind T(X) on (-1, +1),
!!    W(X) = 1 / SQRT(1 - X*X)
!!
!!    3: Chebyshev polynomials of the second kind U(X) on (-1, +1),
!!    W(X) = SQRT(1 - X*X)
!!
!!    4: Hermite polynomials H(X) on (-infinity,+infinity),
!!    W(X) = EXP(-X**2)
!!
!!    5: Jacobi polynomials P(ALPHA,BETA)(X) on (-1, +1),
!!    W(X) = (1-X)**ALPHA + (1+X)**BETA,
!!    ALPHA and BETA greater than -1.
!!
!!    6: Laguerre polynomials, L(ALPHA)(X) on (0, +infinity),
!!    W(X) = EXP(-X) * X**ALPHA,
!!    ALPHA greater than -1.
!!
!!    Input, integer N, specifies the number of coefficients to
!!    calculate.
!!
!!    Input, real ALPHA, the value of the ALPHA parameter,
!!    required only for Jacobi or Laguerre polynomials.
!!
!!    Input, real BETA, the value of the BETA parameter,
!!    required only for Jacobi polynomials.
!!
!!    Output, real B(N-1), the offdiagonal coefficients.
!!
!!    Output, real A(N), the diagonal coefficients.
!!
!!    Output, real MUZERO, the zero-th moment, Integral W(X) DX,
!!    of the polynomial's weight function over its interval of
!!    definition.
!!
!  implicit none
!!
!  integer n
!!
!  real a(n)
!  real abi
!  real alpha
!  real b(n-1)
!  real beta
!  real gamma
!  integer i
!  integer kind
!  real muzero
!  real pi
!!
!!  KIND = 1:
!!
!!  Legendre polynomials P(X) on (-1, +1),
!!  W(X) = 1.
!!
!  if ( kind == 1 ) then
! 
!    muzero = 2.0E+00
! 
!    a(1:n) = 0.0E+00
! 
!    do i = 1, n-1
!      b(i) = real(i) / sqrt(4.0*real(i*i) - 1.0)
!    end do
!!
!!  KIND = 2:
!!
!!  Chebyshev polynomials of the first kind T(X) on (-1, +1),
!!  W(X) = 1 / SQRT(1 - X*X)
!!
!  else if ( kind == 2 ) then
! 
!    muzero = pi()
!    a(1:n) = 0.0E+00
!    b(1) = sqrt ( 0.5E+00 )
!    b(1:n-1) = 0.5E+00
!!
!!  KIND = 3:
!!
!!  Chebyshev polynomials of the second kind U(X) on (-1, +1),
!!  W(X) = SQRT(1 - X*X)
!!
!  else if ( kind == 3 ) then
! 
!    muzero = pi() / 2.0E+00
!    a(1:n) = 0.0E+00
!    b(1:n-1) = 0.5E+00
!!
!!  KIND = 4:
!!
!!  Hermite polynomials H(X) on (-infinity,+infinity),
!!  W(X) = EXP(-X**2)
!!
!  else if ( kind == 4 ) then
! 
!    muzero = sqrt ( pi() )
!    a(1:n) = 0.0E+00
!    do i = 1, n-1
!      b(i) = sqrt ( real(i) / 2.0E+00 )
!    end do
!!
!!  KIND = 5:
!!
!!  Jacobi polynomials P(ALPHA,BETA)(X) on (-1, +1),
!!  W(X) = (1-X)**ALPHA + (1+X)**BETA,
!!  ALPHA and BETA greater than -1
!!
!  else if ( kind == 5 ) then
! 
!    muzero = 2.0**( alpha + beta + 1.0E+00 ) * gamma ( alpha + 1.0 ) &
!      * gamma ( beta + 1.0 ) / gamma ( 2.0 + alpha + beta )
! 
!    do i = 1, n
!      a(i) = ( beta**2 - alpha**2 )/ &
!        ((2.0*(i-1)+alpha+beta) * ( 2.0 * real ( i ) + alpha + beta ) )
!    end do
! 
!    abi = 2.0E+00 + alpha+beta
!    b(1) = sqrt ( 4.0*(1.0E+00 + alpha)*(1.0 + beta)/((abi + 1.0)*abi*abi))
! 
!    do i = 2, n-1
!      abi = real ( 2 * i ) + alpha + beta
!      b(i) = sqrt ( 4.0E+00 * real ( i ) * ( real ( i ) + alpha ) &
!        * ( i + beta ) * ( i + alpha + beta) / &
!        ( ( abi*abi - 1.0E+00 ) * abi * abi ) )
!    end do
!!
!!  KIND = 6:
!!
!!  Laguerre polynomials
!!
!!  L(ALPHA)(X) on (0, +infinity),
!!  W(X) = EXP(-X) * X**ALPHA,
!!  ALPHA greater than -1.
!!
!  else if ( kind == 6 ) then
! 
!    muzero = gamma ( alpha + 1.0E+00 )
! 
!    do i = 1, n
!      a(i) = 2.0E+00 * real ( i ) - 1.0E+00 + alpha
!    end do
! 
!    do i = 1, n-1
!      b(i) = sqrt ( real ( i ) * ( real ( i ) + alpha ) )
!    end do
! 
!  end if
! 
!  return
!end
!subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result )
!!
!!***********************************************************************
!!
!!! CSPINT estimates the integral of a tabulated function.
!!
!!
!!  Discussion:
!!
!!    The routine is given the value of a function F(X) at a set of 
!!    nodes XTAB, and estimates
!!
!!      INTEGRAL (A to B) F(X) DX
!!
!!    by computing the cubic natural spline S(X) that interpolates
!!    F(X) at the nodes, and then computing
!!
!!      INTEGRAL (A to B) S(X) DX
!!
!!    exactly.
!!
!!    Other output from the program includes the definite integral
!!    from X(1) to X(I) of S(X), and the coefficients necessary for
!!    the user to evaluate the spline S(X) at any point.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real FTAB(NTAB), contains the tabulated values of
!!    the function, FTAB(I) = F(XTAB(I)).
!!
!!    Input, real XTAB(NTAB), contains the points at which the
!!    function was evaluated.  The XTAB's must be distinct and
!!    in ascending order.
!!
!!    Input, integer NTAB, the number of entries in FTAB and
!!    XTAB.  NTAB must be at least 3.
!!
!!    Input, real A, lower limit of integration.
!!
!!    Input, real B, upper limit of integration.
!!
!!    Output, real Y(3,NTAB), will contain the coefficients
!!    of the interpolating natural spline over each subinterval.
!!
!!    For XTAB(I) <= X <= XTAB(I+1),
!!
!!      S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I))
!!                   + Y(2,I)*(X-XTAB(I))**2
!!                   + Y(3,I)*(X-XTAB(I))**3
!!
!!    Output, real E(NTAB), E(I) = the definite integral from
!!    XTAB(1) to XTAB(I) of S(X).
!!
!!    Workspace, real WORK(NTAB).
!!
!!    Output, real RESULT, the estimated value of the integral.
!!
!  implicit none
!!
!  integer ntab
!!
!  real a
!  real b
!  real e(ntab)
!  real ftab(ntab)
!  integer i
!  integer j
!  real r
!  real result
!  real s
!  real term
!  real u
!  real work(ntab)
!  real xtab(ntab)
!  real y(3,ntab)
!!
!  if ( ntab < 3 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'CSPINT - Fatal error!'
!    write ( *, '(a,i6)' ) '  NTAB must be at least 3, but input NTAB = ',ntab
!    stop
!  end if
! 
!  do i = 1, ntab-1
! 
!    if ( xtab(i+1) <= xtab(i) ) then
!      write ( *, '(a)' ) ' '
!      write ( *, '(a)' ) 'CSPINT - Fatal error!'
!      write ( *, '(a)' ) '  Nodes not in strict increasing order.'
!      write ( *, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
!      write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
!      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
!      stop
!    end if
! 
!  end do
! 
!  s = 0.0E+00
!  do i = 1, ntab-1
!    r = ( ftab(i+1) - ftab(i) ) / ( xtab(i+1) - xtab(i) )
!    y(2,i) = r - s
!    s = r
!  end do
! 
!  result = 0.0E+00
!  s = 0.0E+00
!  r = 0.0E+00
!  y(2,1) = 0.0E+00
!  y(2,ntab) = 0.0E+00
! 
!  do i = 2, ntab-1
!    y(2,i) = y(2,i)+r*y(2,i-1)
!    work(i) = 2.0E+00 * ( xtab(i-1) - xtab(i+1) ) - r * s
!    s = xtab(i+1) - xtab(i)
!    r = s / work(i)
!  end do
! 
!  do j = 2, ntab-1
!    i = ntab+1-j
!    y(2,i) = ((xtab(i+1)-xtab(i))*y(2,i+1)-y(2,i)) / work(i)
!  end do
! 
!  do i = 1, ntab-1
!    s = xtab(i+1)-xtab(i)
!    r = y(2,i+1)-y(2,i)
!    y(3,i) = r / s
!    y(2,i) = 3.0E+00 * y(2,i)
!    y(1,i) = (ftab(i+1)-ftab(i)) / s-(y(2,i)+r)*s
!  end do
! 
!  e(1) = 0.0E+00
!  do i = 1, ntab-1
!    s = xtab(i+1)-xtab(i)
!    term = (((y(3,i)* 0.25E+00 *s+y(2,i) / 3.0 ) *s+y(1,i)* 0.5 )*s+ftab(i))*s
!    e(i+1) = e(i) + term
!  end do
!!
!!  Determine where the endpoints A and B lie in the mesh of XTAB's.
!!
!  r = a
!  u = 1.0E+00
! 
!  do j = 1, 2
!!
!!  The endpoint is less than or equal to XTAB(1).
!!
!    if ( r <= xtab(1) ) then
!      result = result-u*((r-xtab(1))*y(1,1)*0.5E+00 +ftab(1))*(r-xtab(1))
!!
!!  The endpoint is greater than or equal to XTAB(NTAB).
!!
!    else if ( r >= xtab(ntab) ) then
!
!      result = result-u*(e(ntab)+(r-xtab(ntab))*(ftab(ntab)+ &
!        0.5E+00 *(ftab(ntab-1)+(xtab(ntab)-xtab(ntab-1))*y(1,ntab-1)) &
!        *(r-xtab(ntab))))
!!
!!  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
!!
!    else
!      do i = 1,ntab-1
! 
!        if ( r <= xtab(i+1) ) then
!          r = r-xtab(i)
!          result = result-u*(e(i)+(((y(3,i)*0.25*r+y(2,i)/3.0)*r &
!            +y(1,i)*0.5E+00 )*r+ftab(i))*r)
!          go to 120
!        end if
! 
!      end do
! 
!    end if
! 
!  120   continue
! 
!    u = -1.0E+00
!    r = b
! 
!  end do
! 
!  return
!end
!subroutine cubint ( ftab, xtab, ntab, ia, ib, result, error )
!!
!!***********************************************************************
!!
!!! CUBINT approximates an integral using cubic interpolation of data.
!!
!!
!!  Discussion:
!!
!!    The integral to be approximated is
!! 
!!      INTEGRAL (XTAB(IB) to XTAB(IA)) F(X) DX
!!
!!    The routine estimates the error in integration.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!    P E Gill and G F Miller
!!    An algorithm for the integration of unequally spaced data,
!!    Comput J, Number 15, 1972, pages 80-83.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real FTAB(NTAB), contains the tabulated function
!!    values, FTAB(I) = F(XTAB(I)).
!!
!!    Input, real XTAB(NTAB), contains the points at which the
!!    function was tabulated.  XTAB should contain distinct
!!    values, given in ascending order.
!!
!!    Input, integer NTAB, the number of tabulated points.
!!    NTAB must be at least 4.
!!
!!    Input, integer IA, the entry of XTAB at which integration
!!    is to begin.  IA must be no less than 1 and no greater
!!    than NTAB.
!!
!!    Input, integer IB, the entry of XTAB at which integration
!!    is to end.  IB must be no less than 1 and no greater than
!!    NTAB.
!!
!!    Output, real RESULT, the approximate value of the
!!    integral from XTAB(IA) to XTAB(IB) of the function.
!!
!!    Output, real ERROR, an estimate of the error in
!!    integration.
!!
!  implicit none
!!
!  integer ntab
!!
!  real c
!  real d1
!  real d2
!  real d3
!  real error
!  real ftab(ntab)
!  real h1
!  real h2
!  real h3
!  real h4
!  integer i
!  integer ia
!  integer ib
!  integer ind
!  integer it
!  integer j
!  integer k
!  real r1
!  real r2
!  real r3
!  real r4
!  real result
!  real s
!  real term
!  real xtab(ntab)
!!
!  result = 0.0E+00
!  error = 0.0E+00
! 
!  if ( ia == ib ) then
!    return
!  end if
! 
!  if ( ntab < 4 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!    write ( *, '(a,i6)' ) '  NTAB must be at least 4, but input NTAB = ',ntab
!    stop
!  end if
! 
!  if ( ia < 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!    write ( *, '(a,i6)' ) '  IA must be at least 1, but input IA = ',ia
!    stop
!  end if
! 
!  if ( ia > ntab ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!    write ( *, '(a,i6)' ) '  IA must be <= NTAB, but input IA=',ia
!    stop
!  end if
! 
!  if ( ib < 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!    write ( *, '(a,i6)' ) '  IB must be at least 1, but input IB = ',ib
!    stop
!  end if
! 
!  if ( ib > ntab ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'CUBINT - Fatal error!'
!    write ( *, '(a,i6)' ) '  IB must be <= NTAB, but input IB=',ib
!    stop
!  end if
!!
!!  Temporarily switch IA and IB, and store minus sign in IND
!!  so that, while integration is carried out from low X's
!!  to high ones, the sense of the integral is preserved.
!!
!  if ( ia > ib ) then
!    ind = -1
!    it = ib
!    ib = ia
!    ia = it
!  else
!    ind = 1
!  end if
! 
!  s = 0.0E+00
!  c = 0.0E+00
!  r4 = 0.0E+00
!  j = ntab-2
!  if ( ia < ntab-1 .or. ntab == 4 ) then
!    j=max(3,ia)
!  end if
!
!  k = 4
!  if ( ib > 2 .or. ntab == 4 ) then
!    k=min(ntab,ib+2)-1
!  end if
! 
!  do i = j, k
! 
!    if ( i <= j ) then
! 
!      h2 = xtab(j-1)-xtab(j-2)
!      d3 = (ftab(j-1)-ftab(j-2)) / h2
!      h3 = xtab(j)-xtab(j-1)
!      d1 = (ftab(j)-ftab(j-1)) / h3
!      h1 = h2+h3
!      d2 = (d1-d3)/h1
!      h4 = xtab(j+1)-xtab(j)
!      r1 = (ftab(j+1)-ftab(j)) / h4
!      r2 = (r1-d1) / (h4+h3)
!      h1 = h1+h4
!      r3 = (r2-d2) / h1
! 
!      if ( ia <= 1 ) then
!        result = h2 * (ftab(1)+h2*(0.5*d3-h2*(d2/6.0-(h2+h3+h3)*r3/12.)))
!        s = -h2**3 * (h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0
!      end if
! 
!    else
! 
!      h4 = xtab(i+1)-xtab(i)
!      r1 = (ftab(i+1)-ftab(i))/h4
!      r4 = h4+h3
!      r2 = (r1-d1)/r4
!      r4 = r4+h2
!      r3 = (r2-d2)/r4
!      r4 = (r3-d3)/(r4+h1)
! 
!    end if
! 
!    if ( i > ia .and. i <= ib ) then
! 
!      term = h3*((ftab(i)+ftab(i-1))*0.5-h3*h3*(d2+r2+(h2-h4)*r3) / 12.0 )
!      result = result+term
!      c = h3**3*(2.0E+00 *h3*h3+5.*(h3*(h4+h2) + 2.0 * h2 * h4 ) ) / 120.0E+00
!      error = error+(c+s)*r4
! 
!      if ( i /= j ) then
!        s = c
!      else
!        s = s+c+c
!      end if
! 
!    else
! 
!      error = error+r4*s
! 
!    end if
! 
!    if ( i >= k ) then
! 
!      if ( ib >= ntab ) then
!        term = h4*(ftab(ntab) - h4*(0.5*r1+h4*(r2/6.0 +(h3+h3+h4)*r3/12.)))
!        result = result + term
!        error = error - h4**3 * r4 * &
!          ( h4 * ( 3.0 * h4 + 5.0 * h2 ) &
!          + 10.0 * h3 * ( h2 + h3 + h4 ) ) / 60.0E+00
!      end if
! 
!      if ( ib >= ntab-1 ) error=error+s*r4
!    else
!      h1 = h2
!      h2 = h3
!      h3 = h4
!      d1 = r1
!      d2 = r2
!      d3 = r3
!    end if
! 
!  end do
!!
!!  Restore original values of IA and IB, reverse signs
!!  of RESULT and ERROR, to account for integration
!!  that proceeded from high X to low X.
!!
!  if ( ind /= 1 ) then
!    it = ib
!    ib = ia
!    ia = it
!    result = -result
!    error = -error
!  end if
! 
!  return
!end
!subroutine filon_cos ( ftab, a, b, ntab, t, result )
!!
!!***********************************************************************
!!
!!! FILON_COS uses Filon's method on integrals with a cosine factor.
!!
!!
!!  Discussion:
!!
!!    The integral to be approximated has the form:
!!
!!      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
!!
!!    where T is user specified.
!!
!!    The function is interpolated over each subinterval by
!!    a parabolic arc.
!!
!!  Reference:
!!
!!    Abramowitz and Stegun,
!!    Handbook of Mathematical Functions,
!!    pages 890-891.
!!
!!    S M Chase and L D Fosdick,
!!    Algorithm 353, Filon Quadrature,
!!    Communications of the Association for Computing Machinery,
!!    Volume 12, 1969, pages 457-458.
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Modified:
!!
!!    19 February 2002
!!
!!  Parameters:
!!
!!    Input, real FTAB(NTAB), contains the value of the function
!!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!!
!!    Input, real A, B, the limits of integration.
!!
!!    Input, integer NTAB, the number of data points.
!!    NTAB must be odd, and greater than 1.
!!
!!    Input, real T, the multiplier of the X argument of the cosine.
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!  implicit none
!!
!  integer ntab
!!
!  real a
!  real alpha
!  real b
!  real beta
!  real c2n
!  real c2nm1
!  real cost
!  real ftab(ntab)
!  real gamma
!  real h
!  real result
!  real sint
!  real t
!  real theta
!  real xtab(ntab)
!!
!  if ( a == b ) then
!    result = 0.0E+00
!    return
!  end if
! 
!  if ( ntab <= 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'FILON_COS - Fatal error!'
!    write ( *, '(a)' ) '  NTAB < 2'
!    write ( *, '(a,i6)' ) '  NTAB = ', ntab
!    stop
!  end if
! 
!  if ( mod ( ntab, 2 ) /= 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'FILON_COS - Fatal error!'
!    write ( *, '(a)' ) '  NTAB must be odd.'
!    write ( *, '(a,i6)' ) '  NTAB = ', ntab
!    stop
!  end if
!!
!!  Set up a vector of the NTAB X values.
!! 
!  call rvec_even ( a, b, ntab, xtab )
!
!  h = ( b - a ) / real ( ntab - 1 )
!  theta = t * h
!  sint = sin ( theta )
!  cost = cos ( theta )
!
!  alpha = ( theta**2 + theta * sint * cost &
!    - 2.0E+00 * sint**2 ) / theta**3
!
!  beta = ( 2.0E+00 * theta + 2.0E+00 * theta * cost**2 &
!    - 4.0E+00 * sint * cost ) / theta**3
!
!  gamma = 4.0E+00 * ( sint - theta * cost ) / theta**3
!  
!  c2n = sum ( ftab(1:ntab:2) * cos ( t * xtab(1:ntab:2) ) ) &
!    - 0.5E+00 * ( ftab(ntab) * cos ( t * xtab(ntab) ) &
!                + ftab(1) * cos ( t * xtab(1) ) )
!
!  c2nm1 = sum ( ftab(2:ntab-1:2) * cos ( t * xtab(2:ntab-1:2) ) )
! 
!  result = h * ( &
!      alpha * ( ftab(ntab) * sin ( t * xtab(ntab) ) & 
!              - ftab(1) * sin ( t * xtab(1) ) ) &
!    + beta * c2n &
!    + gamma * c2nm1 )
!
!  return
!end
!subroutine filon_sin ( ftab, a, b, ntab, t, result )
!!
!!***********************************************************************
!!
!!! FILON_SIN uses Filon's method on integrals with a sine factor.
!!
!!
!!  Discussion:
!!
!!    The integral to be approximated has the form
!!
!!      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
!!
!!    where T is user specified.
!!
!!    The function is interpolated over each subinterval by
!!    a parabolic arc.
!!
!!  Reference:
!!
!!    Abramowitz and Stegun,
!!    Handbook of Mathematical Functions,
!!    pages 890-891.
!!
!!    S M Chase and L D Fosdick,
!!    Algorithm 353, Filon Quadrature,
!!    Communications of the Association for Computing Machinery,
!!    Volume 12, 1969, pages 457-458.
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Modified:
!!
!!    19 February 2002
!!
!!  Parameters:
!!
!!    Input, real FTAB(NTAB), contains the value of the function
!!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
!!
!!    Input, real A, B, the limits of integration.
!!
!!    Input, integer NTAB, the number of data points, including the
!!    endpoints.  NTAB must be odd, and greater than 1.
!!
!!    Input, real T, multiplier of the X argument of the sine.
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!  implicit none
!!
!  integer ntab
!!
!  real a
!  real alpha
!  real b
!  real beta
!  real cost
!  real ftab(ntab)
!  real gamma
!  real h
!  real result
!  real s2n
!  real s2nm1
!  real sint
!  real t
!  real theta
!  real xtab(ntab)
!!
!  if ( a == b ) then
!    result = 0.0E+00
!    return
!  end if
! 
!  if ( ntab <= 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
!    write ( *, '(a)' ) '  NTAB < 2'
!    write ( *, '(a,i6)' ) '  NTAB = ',ntab
!    stop
!  end if
! 
!  if ( mod ( ntab, 2 ) /= 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
!    write ( *, '(a)' ) '  NTAB must be odd.'
!    write ( *, '(a,i6)' ) '  NTAB = ',ntab
!    stop
!  end if
!!
!!  Set up a vector of the NTAB X values.
!! 
!  call rvec_even ( a, b, ntab, xtab )
!
!  h = ( b - a ) / real ( ntab - 1 )
!  theta = t * h
!  sint = sin ( theta )
!  cost = cos ( theta )
! 
!  alpha = ( theta**2 + theta * sint * cost &
!    - 2.0E+00 * sint**2 ) / theta**3
!
!  beta = ( 2.0E+00 * theta + 2.0E+00 * theta * cost**2 &
!    - 4.0E+00 * sint * cost ) / theta**3
!
!  gamma = 4.0E+00 * ( sint - theta * cost ) / theta**3
!   
!  s2n = sum ( ftab(1:ntab:2) * sin ( t * xtab(1:ntab:2) ) ) &
!    - 0.5E+00 * ( ftab(ntab) * sin ( t * xtab(ntab) ) &
!                + ftab(1) * sin ( t * xtab(1) ) )
!
!  s2nm1 = sum ( ftab(2:ntab-1:2) * sin ( t * xtab(2:ntab-1:2) ) )
!
!  result = h * ( &
!      alpha * ( ftab(1) * cos ( t * xtab(1) ) &
!              - ftab(ntab) * cos ( t * xtab(ntab) ) ) &
!    + beta * s2n &
!    + gamma * s2nm1 )
! 
!  return
!end
!function gamma ( x )
!!
!!*******************************************************************************
!!
!!! GAMMA calculates the Gamma function for a real argument X.
!!
!!
!!  Definition:
!!
!!    GAMMA(X) = Integral ( 0 <= T <= Infinity ) T**(X-1) EXP(-T) DT
!!
!!  Recursion:
!!
!!    GAMMA(X+1) = X * GAMMA(X)
!!
!!  Special values:
!!
!!    GAMMA(0.5) = SQRT(PI)
!!    If N is a positive integer, GAMMA(N+1) = N!, the standard factorial.
!!
!!  Discussion:
!!
!!    Computation is based on an algorithm outlined in reference 1.
!!    The program uses rational functions that approximate the GAMMA
!!    function to at least 20 significant decimal digits.  Coefficients
!!    for the approximation over the interval (1,2) are unpublished.
!!    Those for the approximation for X .GE. 12 are from reference 2.
!!    The accuracy achieved depends on the arithmetic system, the
!!    compiler, the intrinsic functions, and proper selection of the
!!    machine-dependent constants.
!!
!!  Machine-dependent constants:
!!
!!    BETA: radix for the floating-point representation.
!!    MAXEXP: the smallest positive power of BETA that overflows.
!!    XBIG: the largest argument for which GAMMA(X) is representable
!!      in the machine, i.e., the solution to the equation
!!      GAMMA(XBIG) = BETA**MAXEXP.
!!    XINF: the largest machine representable floating-point number;
!!      approximately BETA**MAXEXP.
!!    EPS: the smallest positive floating-point number such that
!!      1.0+EPS .GT. 1.0.
!!    XMININ: the smallest positive floating-point number such that
!!      1/XMININ is machine representable.
!!
!!    Approximate values for some important machines are:
!!
!!                               BETA       MAXEXP        XBIG
!!
!!    CRAY-1         (S.P.)        2         8191        966.961
!!    Cyber 180/855
!!      under NOS    (S.P.)        2         1070        177.803
!!    IEEE (IBM/XT,
!!      SUN, etc.)   (S.P.)        2          128        35.040
!!    IEEE (IBM/XT,
!!      SUN, etc.)   (D.P.)        2         1024        171.624
!!    IBM 3033       (D.P.)       16           63        57.574
!!    VAX D-Format   (D.P.)        2          127        34.844
!!    VAX G-Format   (D.P.)        2         1023        171.489
!!
!!                               XINF         EPS        XMININ
!!
!!    CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
!!    Cyber 180/855
!!      under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
!!    IEEE (IBM/XT,
!!      SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
!!    IEEE (IBM/XT,
!!      SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
!!    IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
!!    VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
!!    VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!!
!!  Reference:
!!
!!    W J Cody,
!!    "An Overview of Software Development for Special Functions",
!!    Lecture Notes in Mathematics, 506,
!!    Numerical Analysis Dundee, 1975,
!!    G. A. Watson (ed.),
!!    Springer Verlag, Berlin, 1976.
!!
!!    Hart et al,
!!    Computer Approximations,
!!    Wiley and sons, New York, 1968.
!!
!!  Author:
!!
!!    W. J. Cody and L. Stoltz,
!!    Applied Mathematics Division,
!!    Argonne National Laboratory,
!!    Argonne, Illinois, 60439.
!!
!!  Parameters:
!!
!!    Input, real X, the argument of the function.
!!
!!    Output, real GAMMA, the value of the function.  The program
!!    returns the value XINF for singularities or when overflow would occur.
!!    The computation is believed to be free of underflow and overflow.
!!
!  implicit none
!!
!  real, parameter, dimension ( 7 ) :: c = (/ &
!    -1.910444077728E-03, &
!     8.4171387781295E-04, &
!    -5.952379913043012E-04, &
!     7.93650793500350248E-04, &
!    -2.777777777777681622553E-03, &
!     8.333333333333333331554247E-02, &
!     5.7083835261E-03 /)
!  real, parameter :: EPS = 1.19E-07
!  real fact
!  real gamma
!  integer i
!  integer n
!  real, parameter, dimension ( 8 ) :: p = (/ &
!    -1.71618513886549492533811E+00, &
!     2.47656508055759199108314E+01, &
!    -3.79804256470945635097577E+02, &
!     6.29331155312818442661052E+02, &
!     8.66966202790413211295064E+02, &
!    -3.14512729688483675254357E+04, &
!    -3.61444134186911729807069E+04, &
!     6.64561438202405440627855E+04 /)
!  logical parity
!  real, parameter :: PI = &
!    3.14159265358979323846264338327950288419716939937510E+00
!  real, parameter, dimension ( 8 ) :: q = (/ &
!    -3.08402300119738975254353e+01, &
!     3.15350626979604161529144e+02, &
!    -1.01515636749021914166146e+03, &
!    -3.10777167157231109440444e+03, &
!     2.25381184209801510330112e+04, &
!     4.75584627752788110767815e+03, &
!    -1.34659959864969306392456e+05, &
!    -1.15132259675553483497211e+05 /)
!  real, parameter :: SQRTPI = 0.9189385332046727417803297E+00
!  real sum1
!  real x
!  real, parameter :: XBIG = 35.040E+00
!  real xden
!  real, parameter :: XINF = 3.4E+38
!  real, parameter :: XMININ = 1.18E-38
!  real xnum
!  real y
!  real y1
!  real ysq
!  real z
!!
!  parity = .false.
!  fact = 1.0E+00
!  n = 0
!  y = x
!!
!!  Argument is negative.
!!
!  if ( y <= 0.0E+00 ) then
!
!    y = - x
!    y1 = aint ( y )
!    gamma = y - y1
!
!    if ( gamma /= 0.0E+00 ) then
!
!      if ( y1 /= aint ( y1 * 0.5E+00 ) * 2.0E+00 ) then
!        parity = .true.
!      end if
!
!      fact = - PI / sin ( PI * gamma )
!      y = y + 1.0E+00
!
!    else
!
!      gamma = XINF
!      return
!
!    end if
!
!  end if
!!
!!  Argument < EPS
!!
!  if ( y < EPS ) then
!
!    if (y >= XMININ ) then
!      gamma = 1.0E+00 / y
!    else
!      gamma = XINF
!      return
!    end if
!
!  else if ( y < 12.0E+00 ) then
!
!    y1 = y
!!
!!  0.0E+00 < argument < 1.0E+00
!!
!    if ( y < 1.0E+00 ) then
!      z = y
!      y = y + 1.0E+00
!!
!!  1.0E+00 < argument < 12.0, reduce argument if necessary.
!!
!    else
!      n = int ( y ) - 1
!      y = y - real ( n )
!      z = y - 1.0E+00
!    end if
!!
!!  Evaluate approximation for 1.0E+00 < argument < 2.0.
!!
!    xnum = 0.0E+00
!    xden = 1.0E+00
!    do i = 1, 8
!      xnum = ( xnum + p(i) ) * z
!      xden = xden * z + q(i)
!    end do
!
!    gamma = xnum / xden + 1.0E+00
!!
!!  Adjust result for case  0.0E+00 < argument < 1.0.
!!
!    if ( y1 < y ) then
!      gamma = gamma / y1
!!
!!  Adjust result for case  2.0E+00 < argument < 12.0.
!!
!    else if ( y1 > y ) then
!
!      do i = 1, n
!        gamma = gamma * y
!        y = y + 1.0E+00
!      end do
!
!    end if
!!
!!  Evaluate for 12 <= argument.
!!
!  else
!
!    if ( y <= XBIG ) then
!
!      ysq = y**2
!      sum1 = c(7)
!      do i = 1, 6
!        sum1 = sum1 / ysq + c(i)
!      end do
!      sum1 = sum1 / y - y + SQRTPI
!      sum1 = sum1 + ( y - 0.5E+00 ) * log ( y )
!      gamma = exp ( sum1 )
!
!    else
!
!      gamma = XINF
!      return
!
!    end if
!
!  end if
!!
!!  Final adjustments and return.
!!
!  if ( parity ) then
!    gamma = - gamma
!  end if
!
!  if ( fact /= 1.0E+00 ) then
!    gamma = fact / gamma
!  end if
!
!  return
!end
!subroutine gaus8 ( func, a, b, err, result, ierr )
!!
!!***********************************************************************
!!
!!! GAUS8 estimates the integral of a function.
!!
!!
!!  Discussion:
!!
!!    GAUS8 integrates real functions of one variable over finite
!!    intervals using an adaptive 8-point Legendre-Gauss
!!    algorithm.
!!
!!    GAUS8 is intended primarily for high accuracy integration or
!!    integration of smooth functions.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Author:
!!
!!    R E Jones,
!!    Sandia National Laboratory,
!!    Los Alamos, New Mexico
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real, external FUNC, name of external function to be integrated.
!!    This name must be in an external statement in the calling program.
!!    FUNC must be a function of one real argument.  The value
!!    of the argument to FUNC is the variable of integration
!!    which ranges from A to B.
!!
!!    Input, real A, the lower limit of integration.
!!
!!    Input, real B, the upper limit of integration.
!!
!!    Input/output, real ERR.
!!    On input, ERR is a requested pseudorelative error
!!    tolerance.  Normally pick a value of ABS ( ERR ) so that
!!    STOL < ABS ( ERR ) <= 1.0E-3 where STOL is the single precision
!!    unit roundoff  = R1MACH(4).
!!
!!    RESULT will normally have no more error than
!!    ABS ( ERR ) times the integral of the absolute value of
!!    FUN(X).  Usually, smaller values for ERR yield more
!!    accuracy and require more function evaluations.
!!
!!    A negative value for ERR causes an estimate of the
!!    absolute error in RESULT to be returned in ERR.  Note that
!!    ERR must be a variable (not a constant) in this case.
!!    Note also that the user must reset the value of ERR
!!    before making any more calls that use the variable ERR.
!!
!!    On output, ERR will be an estimate of the absolute error
!!    in RESULT if the input value of ERR was negative.  ERR is
!!    unchanged if the input value of ERR was non-negative.
!!
!!    The estimated error is solely for information to the user
!!    and should not be used as a correction to the computed integral.
!!
!!    Output, real RESULT, the computed value of the integral.
!!
!!    Output, integer IERR, a status code.
!!
!!    Normal Codes:
!!
!!     1 RESULT most likely meets requested error tolerance, or A = B.
!!    -1 A and B are too nearly equal to allow normal
!!        integration.  RESULT is set to zero.
!!
!!     Abnormal Code:
!!
!!     2 RESULT probably does not meet requested error tolerance.
!!
!  implicit none
!!
!  real a
!  real aa(30)
!  real ae
!  real anib
!  real area
!  real b
!  real c
!  real ce
!  real ee
!  real ef
!  real eps
!  real err
!  real est
!  real, external :: func
!  real g8
!  real gl
!  real glr
!  real gr(30)
!  real h
!  real hh(30)
!  integer i1mach
!  integer, save :: icall = 0
!  integer ierr
!  integer k
!  integer, save :: kml = 6
!  integer, save :: kmx = 5000
!  integer l
!  integer lmn
!  integer lmx
!  integer lr(30)
!  integer mxl
!  integer nbits
!  integer nib
!  integer, save :: nlmn = 1
!  integer nlmx
!  real r1mach
!  real result
!  real tol
!  real vl(30)
!  real vr
!  real, save :: w1
!  real, save :: w2
!  real, save :: w3
!  real, save :: w4
!  real x
!  real, save :: x1
!  real, save :: x2
!  real, save :: x3
!  real, save :: x4
!!
!  data x1, x2, x3, x4/ &
!          1.83434642495649805e-01,     5.25532409916328986e-01, &
!          7.96666477413626740e-01,     9.60289856497536232e-01/
!!
!  data w1, w2, w3, w4/ &
!          3.62683783378361983e-01,     3.13706645877887287e-01, &
!          2.22381034453374471e-01,     1.01228536290376259e-01/
!!
!!  Warning!  Statement function!
!!
!  g8(x,h) = h*((w1*(func(x-x1*h) + func(x+x1*h)) &
!             +w2*(func(x-x2*h) + func(x+x2*h))) &
!            +(w3*(func(x-x3*h) + func(x+x3*h)) &
!             +w4*(func(x-x4*h) + func(x+x4*h))))
!!
!  if ( a == b ) then
!    err = 0.0E+00
!    result = 0.0E+00
!    return
!  end if
! 
!  if ( icall /= 0 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'GAUS8 - Fatal error!'
!    write ( *, '(a)' ) '  GAUS8 was called recursively.'
!    stop
!  end if
!
!  icall = 1
!  k = i1mach(11)
!  anib = r1mach(5) * real(i1mach(11)) / 0.30102000
!  nbits = int(anib)
!  nlmx = min ( 30, (nbits*5)/8 )
!  result = 0.0E+00
!  ierr = 1
!  ce = 0.0E+00
!  result = 0.0E+00
!  lmx = nlmx
!  lmn = nlmn
! 
!  if ( b /= 0.0 ) then
!
!    if ( sign ( 1.0, b ) * a <= 0.0E+00 ) then
!      go to 10
!    end if
!
!    c = abs ( 1.0 - a / b )
!    if ( c > 0.1 ) go to 10
!    if ( c <= 0.0E+00 ) go to 140
!    anib = 0.5 - log(c) / 0.69314718E+00
!    nib = int(anib)
!    lmx = min(nlmx,nbits-nib-7)
!    if ( lmx < 1 ) go to 130
!    lmn = min ( lmn, lmx )
!
!  end if
! 
!10    continue
! 
!  tol = max ( abs ( err ), 2.0**(5-nbits)) / 2.0E+00
!  if ( err == 0.0 ) then
!    tol=sqrt( epsilon ( 1.0E+00 ) )
!  end if
!
!  eps = tol
!  hh(1) = (b-a) / 4.0E+00
!  aa(1) = a
!  lr(1) = 1
!  l = 1
!  est = g8 ( aa(l) + 2.0*hh(l),2.0*hh(l) )
!  k = 8
!  area = abs ( est )
!  ef = 0.5
!  mxl = 0
!!
!!  Compute refined estimates, estimate the error, etc.
!!
!20 continue
! 
!  gl = g8 ( aa(l)+hh(l),hh(l) )
!  gr(l) = g8(aa(l)+3.0*hh(l),hh(l))
!  k = k + 16
!  area = area + ( abs ( gl ) + abs ( gr(l) ) - abs ( est ) )
! 
!  glr = gl + gr(l)
!  ee = abs ( est - glr ) * ef
!  ae = max ( eps * area, tol * abs ( glr ) )
!
!  if ( ee - ae <= 0.0E+00 ) then
!    go to 40
!  else
!    go to 50
!  end if
! 
!30 continue
! 
!  mxl = 1
! 
!40 continue
! 
!  ce = ce + (est-glr)
! 
!  if ( lr(l) <= 0 ) then
!    go to 60
!  else
!    go to 80
!  end if
!!
!!  Consider the left half of this level
!!
!50 continue
!
!  if ( k > kmx ) lmx = kml
!  if ( l >= lmx ) go to 30
!  l = l + 1
!  eps = eps * 0.5E+00
!  ef = ef / sqrt ( 2.0E+00 )
!  hh(l) = hh(l-1) * 0.5E+00
!  lr(l) = -1
!  aa(l) = aa(l-1)
!  est = gl
!  go to 20
!!
!!  Proceed to right half at this level
!!
!60 continue
!
!  vl(l) = glr
!
!70 continue
!
!  est = gr(l-1)
!  lr(l) = 1
!  aa(l) = aa(l) + 4.0E+00 * hh(l)
!  go to 20
!!
!!  Return one level
!!
!80 continue
!
!  vr = glr
!
!90 continue
!
!  if ( l <= 1 ) go to 120
!  l = l - 1
!  eps = eps * 2.0E+00
!  ef = ef * sqrt ( 2.0E+00 )
! 
!  if ( lr(l) <= 0 ) then
!    vl(l) = vl(l+1) + vr
!    go to 70
!  else
!    vr = vl(l+1) + vr
!    go to 90
!  end if
!!
!!  Exit
!!
!120   continue
! 
!  result = vr
!  if ( mxl == 0 .or. abs ( ce ) <= 2.0E+00 * tol * area ) go to 140
!  ierr = 2
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) 'GAUS8 - Warning!'
!  write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
!  icall = 0
!
!  if ( err < 0.0E+00 ) then
!    err = ce
!  end if
!
!  return
! 
!130   continue
! 
!  ierr = -1
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) 'GAUS8 - Warning!'
!  write ( *, '(a)' ) '  A and B are too close to carry out integration.'
! 
!140   continue
! 
!  icall = 0
! 
!  if ( err < 0.0E+00 ) then
!    err = ce
!  end if
! 
!  return
!end
!subroutine gausq2 ( n, d, e, z, ierr )
!!
!!***********************************************************************
!!
!!! GAUSQ2 finds the eigenvalues of a symmetric tridiagonal matrix.
!!
!!
!!  Discussion:
!!
!!    GAUSQ2 finds the eigenvalues and first components of the
!!    eigenvectors of a symmetric tridiagonal matrix by the implicit QL
!!    method.
!!
!!    GAUSQ2 is a translation of an ALGOL procedure,
!!
!!      Martin and Wilkinson,
!!      Numerische Mathematik,
!!      Volume 12, pages 377-383, 1968,
!!
!!    as modified by
!!
!!      Dubrulle,
!!      Numerische Mathematik,
!!      Volume 15, page 450, 1970.
!!
!!    Handbook for Automatic Computation,
!!    vol.ii-linear algebra,
!!    pages 241-248, 1971.
!!
!!    GAUSQ2 is a modified version of the EISPACK routine IMTQL2.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, integer N, is the order of the matrix.
!!
!!    Input/output, real D(N).
!!    On input, D contains the diagonal elements of the matrix.
!!    On output, D contains the eigenvalues in ascending order.
!!    If an error exit is made, the eigenvalues are correct but
!!    unordered for indices 1, 2, ..., IERR-1;
!!
!!    Input/output, real E(N).
!!    On input, E contains the subdiagonal elements of the input matrix
!!    in its first N-1 positions.  E(N) is arbitrary.
!!    On output, E has been destroyed.
!!
!!    Input/output, real Z(N).
!!    On input, Z contains the first row of the identity matrix.
!!    On output, Z contains the first components of the orthonormal
!!    eigenvectors of the symmetric tridiagonal matrix.  If an error exit is
!!    made, Z contains the eigenvectors associated with the stored
!!    eigenvalues.
!!
!!    Output, integer IERR.
!!    0, for normal return,
!!    J, if the j-th eigenvalue has not been determined after 30 iterations.
!!
!  implicit none
!!
!  integer n
!!
!  real b
!  real c
!  real d(n)
!  real e(n)
!  real epmach
!  real f
!  real g
!  integer i
!  integer ierr
!  integer ii
!  integer j
!  integer k
!  integer l
!  integer m
!  integer mml
!  real p
!  real r
!  real s
!  real z(n)
!!
!  epmach = epsilon ( epmach )
! 
!  ierr = 0
!
!  if ( n == 1 ) then
!    return
!  end if
! 
!  e(n) = 0.0E+00
! 
!  do l = 1, n
! 
!    j = 0
!!
!!  Look for a small sub-diagonal element
!!
!    do m = l, n
!
!      if ( m == n ) then
!        exit
!      end if
!
!      if ( abs ( e(m) ) <= epmach * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
!        exit
!      end if
!
!    end do
! 
!10  continue
! 
!    p = d(l)
!    if ( m == l ) go to 20
! 
!    if ( j == 30 ) then
!      ierr = l
!      return
!    end if
! 
!    j = j + 1
!!
!!  Form shift
!!
!    g = (d(l+1) - p) / (2.0E+00 * e(l))
!    r = sqrt ( g*g + 1.0E+00 )
!    g = d(m) - p + e(l) / (g + sign(r, g))
!    s = 1.0E+00
!    c = 1.0E+00
!    p = 0.0E+00
!    mml = m - l
! 
!    do ii = 1, mml
! 
!      i = m - ii
!      f = s * e(i)
!      b = c * e(i)
! 
!      if ( abs ( f ) >= abs ( g ) ) then
! 
!        c = g / f
!        r = sqrt ( c * c + 1.0E+00 )
!        e(i+1) = f * r
!        s = 1.0E+00 / r
!        c = c * s
! 
!      else
! 
!        s = f / g
!        r = sqrt ( s*s + 1.0E+00 )
!        e(i+1) = g * r
!        c = 1.0E+00 / r
!        s = s * c
! 
!      end if
! 
!      g = d(i+1) - p
!      r = (d(i) - g) * s + 2.0E+00 * c * b
!      p = s * r
!      d(i+1) = g + p
!      g = c * r - b
!!
!!  Form the first component of the vector.
!!
!      f = z(i+1)
!      z(i+1) = s * z(i) + c * f
!      z(i) = f * z(i) - s * f
!    end do
! 
!    d(l) = d(l) - p
!    e(l) = g
!    e(m) = 0.0E+00
!    go to 10
! 
!20  continue
! 
!  end do
!!
!!  Order the eigenvalues and eigenvectors.
!!
!  do ii = 2, n
! 
!    i = ii - 1
!    k = i
!    p = d(i)
! 
!    do j = ii, n
!      if ( d(j) < p ) then
!        k = j
!        p = d(j)
!      end if
!    end do
! 
!    if ( k /= i ) then
!      d(k) = d(i)
!      d(i) = p
!      call r_swap ( z(i), z(k) )
!    end if
! 
!  end do
! 
!  return
!end
!subroutine gaussq ( kind, norder, alpha, beta, kpts, endpts, b, xtab, &
!  weight )
!!
!!***********************************************************************
!!
!!! GAUSSQ computes a Gauss quadrature rule.
!!
!!
!!  Discussion:
!!
!!    GAUSSQ computes the nodes and weights for Gaussian-type quadrature
!!    rules with pre-assigned nodes.
!!
!!    These are used when one wishes to approximate
!!
!!      Integral (from A to B)  F(X) W(X) DX
!!
!!    by
!!
!!      Sum (J = 1 to NORDER) WEIGHT(I)*F(XTAB(I))
!!
!!
!!    GAUSSQ includes six integration rules that are applicable
!!    to this problem, for particular weight functions and particular
!!    intervals, including infinite and semi-infinite intervals.
!!
!!    Associated with each weight function W(X) is a set of
!!    orthogonal polynomials.  The nodes XTAB are just the zeroes
!!    of the proper NORDER-th degree polynomial.
!!
!!    GAUSSQ allows the user to modify the rule to require that
!!    one or both of the endpoints of the interval are to be
!!    included as quadrature nodes.
!!
!!  References:
!!
!!    Golub, G H, and Welsch, J H,
!!    Calculation of Gaussian Quadrature Rules,
!!    Mathematics of Computation
!!    Volume 23, April 1969, pages 221-230.
!!
!!    Golub, G H,
!!    Some Modified Matrix Eigenvalue Problems,
!!    SIAM Review
!!    Volume 15, April 1973, pages 318-334, section 7.
!!
!!    Stroud and Secrest,
!!    Gaussian Quadrature Formulas,
!!    Prentice-Hall,
!!    Englewood Cliffs, New Jersey, 1966.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, integer KIND, chooses the quadrature rule to be calculated.
!!    1:  Legendre quadrature,
!!         W(X) = 1
!!         on (-1, 1)
!!
!!    2:  Chebyshev quadrature of the first kind
!!         W(X) = 1/SQRT(1 - X*X)
!!         on (-1, +1)
!!
!!    3:  Chebyshev quadrature of the second kind
!!         W(X) = SQRT(1 - X*X)
!!         on (-1, 1)
!!
!!    4:  Hermite quadrature,
!!         W(X) = EXP(-X*X)
!!         on (-infinity, +infinity)
!!
!!    5:  Jacobi quadrature,
!!         W(X) = (1-X)**ALPHA * (1+X)**BETA
!!         on (-1, 1),
!!         ALPHA > -1, BETA > -1.
!!
!!         Note that KIND = 2 and 3 are a special case of this.
!!
!!    6:  Generalized Laguerre quadrature,
!!         W(X) = EXP(-X)*X**ALPHA
!!         on (0, +infinity),
!!         ALPHA > -1
!!
!!    Input, integer NORDER, the number of points used for the quadrature rule.
!!
!!    Input, real ALPHA.
!!    ALPHA is only required for Gauss-Jacobi and Gauss-Laguerre
!!    quadrature.  Its value is ignored in other cases.
!!
!!    Input, real BETA.
!!    BETA is only required for Gauss-Jacobi quadrature.
!!    Its value is ignored in other cases.
!!
!!    Input, integer KPTS.
!!    KPTS is normally zero.
!!
!!    If KPTS is nonzero, it signals that one or both of the
!!    endpoints of the interval is required to be a node.
!!    This is called Gauss-Radau or Gauss-Lobatto quadrature.
!!    Then KPTS is the number of endpoints that must be
!!    included, either 1 or 2.
!!
!!    Input, real ENDPTS(2).
!!    If KPTS is 1 or 2, ENDPTS contains the locations of the
!!    endpoints to be fixed.
!!
!!    Workspace, real B(NORDER).
!!
!!    Output, real XTAB(NORDER).
!!    XTAB contains the nodes for the quadrature rule.
!!
!!    Output, real WEIGHT(NORDER).
!!    WEIGHT contains the weights for the quadrature rule.
!!
!  implicit none
!!
!  integer norder
!!
!  real alpha
!  real b(norder)
!  real beta
!  real endpts(2)
!  real gam
!  integer i
!  integer ierr
!  integer kind
!  integer kpts
!  real muzero
!  real solve
!  real t1
!  real weight(norder)
!  real xtab(norder)
!!
!!  Get the diagonal coefficients XTAB(1:NORDER) and off-diagonal
!!  coefficients B(1:NORDER-1) and MUZERO.
!!
!  call class ( kind, norder, alpha, beta, b, xtab, muzero )
!!
!!  The matrix of coefficients is assumed to be symmetric.
!!  The array XTAB contains the diagonal elements, the array
!!  B the off-diagonal elements.
!!  Make appropriate changes in the lower right 2 by 2 submatrix.
!!
!!  If KPTS = 1, only XTAB(NORDER) must be changed.
!!
!  if ( kpts == 1 ) then
! 
!    xtab(norder) = endpts(1) + solve(endpts(1),norder,xtab,b) * b(norder-1)**2
!!
!!  If KPTS = 2, XTAB(NORDER) and B(NORDER-1) must be recomputed.
!!
!  else if ( kpts == 2 ) then
! 
!    gam = solve ( endpts(1), norder, xtab, b )
!    t1 = ((endpts(1) - endpts(2)) / (solve(endpts(2),norder,xtab,b) - gam ) )
!    b(norder-1) = sqrt(t1)
!    xtab(norder) = endpts(1) + gam * t1
! 
!  end if
!!
!!  The indices of the elements of B run from 1 to NORDER-1.
!!  The value of B(NORDER) is of no importance.
!!
!!  Now compute the eigenvalues of the symmetric tridiagonal
!!  matrix, which has been modified as necessary.
!!
!!  The method used is a QL-type method with origin shifting.
!!
!  weight(1) = 1.0E+00
!  weight(2:norder) = 0.0E+00
! 
!  call gausq2 ( norder, xtab, b, weight, ierr )
! 
!  do i = 1, norder
!    weight(i) = muzero * weight(i)**2
!  end do
! 
!  return
!end
!subroutine hiordq ( n, y, delt, work, result )
!!
!!***********************************************************************
!!
!!! HIORDQ approximates the integral of a function using equally spaced data.
!!
!!
!!  Discussion:
!!
!!    The method applies the trapezoidal rule to various subsets of the
!!    data, and then applies Richardson extrapolation.
!!
!!  Author:
!!
!!    Alan Kaylor Cline,
!!    Department of Computer Science,
!!    University of Texas at Austin.
!!
!!  Modified:
!!
!!    19 February 2002
!!
!!  Parameters:
!!
!!    Input, integer N, number of data points.
!!
!!    Input, real Y(N), the Y values of the data.
!!
!!    Input, real DELT, the spacing between the X values of the
!!    data.  The actual X values are not needed!
!!
!!    Work array, real WORK(2*(N-1)).  The actual minimum amount
!!    of workspace required is two times the number of integer
!!    divisors of N-1.
!!
!!    Output, real RESULT, the approximation to the integral.
!!
!  implicit none
!!
!  integer n
!!
!  real delt
!  real fac
!  integer i
!  integer j
!  integer jbak
!  integer jj
!  integer k
!  real result
!  real sum2
!  real sum1
!  real work(2*(n-1))
!  real y(n)
!!
!!  Determine initial trapezoidal rule
!!
!  sum1 = ( y(1) + y(n) ) / 2.0E+00
!  j = -1
! 
!  do k = 1, n-1
!!
!!  Check if K divides N-1
!!
!    if ( ((n-1)/k)*k == n-1 ) then
!!
!!  Determine the K-point trapezoidal rule.
!!
!      sum2 = -sum1
!      do i = 1, n, (n-1)/k
!        sum2 = sum2 + y(i)
!      end do
! 
!      j = j + 2
!      work(j) = delt * sum2 * real ( ( n - 1 ) / k )
!      work(j+1) = real ( ((n-1)/k)**2 )
!!
!!  Apply Richardson extrapolation.
!!
!      if ( k /= 1 ) then
! 
!        do jj = 3, j, 2
!          jbak = j+1-jj
!          fac = work(j+1) / ( work(j+1) - work(jbak+1) )
!          work(jbak) = work(jbak+2) + fac * ( work(jbak) - work(jbak+2) )
!        end do
! 
!      end if
! 
!    end if
! 
!  end do
! 
!  result = work(1)
! 
!  return
!end
!function i1mach ( i )
!!
!!*******************************************************************************
!!
!!! I1MACH returns integer machine constants.
!!
!!
!!  I/O unit numbers.
!!
!!    I1MACH(1) = the standard input unit.
!!    I1MACH(2) = the standard output unit.
!!    I1MACH(3) = the standard punch unit.
!!    I1MACH(4) = the standard error message unit.
!!
!!  Words.
!!
!!    I1MACH(5) = the number of bits per integer storage unit.
!!    I1MACH(6) = the number of characters per integer storage unit.
!!
!!  Integers.
!!
!!  Assume integers are represented in the S digit base A form:
!!
!!  Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
!!  where 0 <= X(I)<A for I=0 to S-1.
!!
!!    I1MACH(7) = A, the base.
!!    I1MACH(8) = S, the number of base A digits.
!!    I1MACH(9) = A**S-1, the largest integer.
!!
!!  Floating point numbers
!!
!!  Assume floating point numbers are represented in the T digit base B form:
!!
!!    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
!!
!!  where 0 <= X(I)<B for I=1 to T, 0<X(1) and EMIN <= E <= EMAX
!!
!!    I1MACH(10) = B, the base.
!!
!!  Single precision
!!
!!    I1MACH(11) = T, the number of base B digits.
!!    I1MACH(12) = EMIN, the smallest exponent E.
!!    I1MACH(13) = EMAX, the largest exponent E.
!!
!!  Double precision
!!
!!    I1MACH(14) = T, the number of base B digits.
!!    I1MACH(15) = EMIN, the smallest exponent E.
!!    I1MACH(16) = EMAX, the largest exponent E.
!!
!!  To alter this function for a particular environment, the desired set of DATA
!!  statements should be activated by removing the C from column 1.  On rare
!!  machines, a STATIC statement may need to be added, but probably more systems
!!  prohibit than require it.
!!
!!  Also, the values of I1MACH(1) through I1MACH(4) should be checked for
!!  consistency with the local operating system.  For FORTRAN 77, you may wish
!!  to adjust the data statement so imach(6) is set to 1, and then to comment
!!  out the executable test on I.EQ.6 below.
!!
!!  For IEEE-arithmetic machines (binary standard), the first set of constants
!!  below should be appropriate, except perhaps for IMACH(1) - IMACH(4).
!!
!  integer i
!  integer i1mach
!  integer imach(16)
!  integer output
!!
!  equivalence (imach(4),output)
!!
!!  IEEE arithmetic machines, such as the ATT 3B series, Motorola
!!  68000 based machines such as the SUN 3 and ATT PC 7300, and
!!  8087 based micros such asthe IBM PC and ATT 6300.
!!
!   data imach( 1) /    5 /
!   data imach( 2) /    6 /
!   data imach( 3) /    7 /
!   data imach( 4) /    6 /
!   data imach( 5) /   32 /
!   data imach( 6) /    4 /
!   data imach( 7) /    2 /
!   data imach( 8) /   31 /
!   data imach( 9) / 2147483647 /
!   data imach(10) /    2 /
!   data imach(11) /   24 /
!   data imach(12) / -125 /
!   data imach(13) /  128 /
!   data imach(14) /   53 /
!   data imach(15) / -1021 /
!   data imach(16) /  1024 /
!!
!!  ALLIANT FX/8 UNIX FORTRAN compiler.
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /     6 /
!!      data imach( 4) /     0 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     4 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    32 /
!!      data imach( 9) /2147483647/
!!      data imach(10) /     2 /
!!      data imach(11) /    24 /
!!      data imach(12) /  -126 /
!!      data imach(13) /   128 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1022 /
!!      data imach(16) /  1024 /
!!
!!  AMDAHL machines.
!!
!!      data imach( 1) /   5 /
!!      data imach( 2) /   6 /
!!      data imach( 3) /   7 /
!!      data imach( 4) /   6 /
!!      data imach( 5) /  32 /
!!      data imach( 6) /   4 /
!!      data imach( 7) /   2 /
!!      data imach( 8) /  31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /  16 /
!!      data imach(11) /   6 /
!!      data imach(12) / -64 /
!!      data imach(13) /  63 /
!!      data imach(14) /  14 /
!!      data imach(15) / -64 /
!!      data imach(16) /  63 /
!!
!!  BURROUGHS 1700 system.
!!
!!      data imach( 1) /    7 /
!!      data imach( 2) /    2 /
!!      data imach( 3) /    2 /
!!      data imach( 4) /    2 /
!!      data imach( 5) /   36 /
!!      data imach( 6) /    4 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   33 /
!!      data imach( 9) / Z1FFFFFFFF /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -256 /
!!      data imach(13) /  255 /
!!      data imach(14) /   60 /
!!      data imach(15) / -256 /
!!      data imach(16) /  255 /
!!
!!  BURROUGHS 5700 system.
!!
!!      data imach( 1) /   5 /
!!      data imach( 2) /   6 /
!!      data imach( 3) /   7 /
!!      data imach( 4) /   6 /
!!      data imach( 5) /  48 /
!!      data imach( 6) /   6 /
!!      data imach( 7) /   2 /
!!      data imach( 8) /  39 /
!!      data imach( 9) / O0007777777777777 /
!!      data imach(10) /   8 /
!!      data imach(11) /  13 /
!!      data imach(12) / -50 /
!!      data imach(13) /  76 /
!!      data imach(14) /  26 /
!!      data imach(15) / -50 /
!!      data imach(16) /  76 /
!!
!!  BURROUGHS 6700/7700 systems.
!!
!!      data imach( 1) /   5 /
!!      data imach( 2) /   6 /
!!      data imach( 3) /   7 /
!!      data imach( 4) /   6 /
!!      data imach( 5) /  48 /
!!      data imach( 6) /   6 /
!!      data imach( 7) /   2 /
!!      data imach( 8) /  39 /
!!      data imach( 9) / O0007777777777777 /
!!      data imach(10) /   8 /
!!      data imach(11) /  13 /
!!      data imach(12) / -50 /
!!      data imach(13) /  76 /
!!      data imach(14) /  26 /
!!      data imach(15) / -32754 /
!!      data imach(16) /  32780 /
!!
!!  CDC CYBER 170/180 series using NOS
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   60 /
!!      data imach( 6) /   10 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   48 /
!!      data imach( 9) / O"00007777777777777777" /
!!      data imach(10) /    2 /
!!      data imach(11) /   48 /
!!      data imach(12) / -974 /
!!      data imach(13) / 1070 /
!!      data imach(14) /   96 /
!!      data imach(15) / -927 /
!!      data imach(16) / 1070 /
!!
!!  CDC CYBER 170/180 series using NOS/VE
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /     7 /
!!      data imach( 4) /     6 /
!!      data imach( 5) /    64 /
!!      data imach( 6) /     8 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    63 /
!!      data imach( 9) / 9223372036854775807 /
!!      data imach(10) /     2 /
!!      data imach(11) /    47 /
!!      data imach(12) / -4095 /
!!      data imach(13) /  4094 /
!!      data imach(14) /    94 /
!!      data imach(15) / -4095 /
!!      data imach(16) /  4094 /
!!
!!  CDC CYBER 200 series
!!
!!      data imach( 1) /      5 /
!!      data imach( 2) /      6 /
!!      data imach( 3) /      7 /
!!      data imach( 4) /      6 /
!!      data imach( 5) /     64 /
!!      data imach( 6) /      8 /
!!      data imach( 7) /      2 /
!!      data imach( 8) /     47 /
!!      data imach( 9) / X'00007FFFFFFFFFFF' /
!!      data imach(10) /      2 /
!!      data imach(11) /     47 /
!!      data imach(12) / -28625 /
!!      data imach(13) /  28718 /
!!      data imach(14) /     94 /
!!      data imach(15) / -28625 /
!!      data imach(16) /  28718 /
!!
!!  CDC 6000/7000 series using FTN4.
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   60 /
!!      data imach( 6) /   10 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   48 /
!!      data imach( 9) / 00007777777777777777B /
!!      data imach(10) /    2 /
!!      data imach(11) /   47 /
!!      data imach(12) / -929 /
!!      data imach(13) / 1070 /
!!      data imach(14) /   94 /
!!      data imach(15) / -929 /
!!      data imach(16) / 1069 /
!!
!!  CDC 6000/7000 series using FTN5.
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   60 /
!!      data imach( 6) /   10 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   48 /
!!      data imach( 9) / O"00007777777777777777" /
!!      data imach(10) /    2 /
!!      data imach(11) /   47 /
!!      data imach(12) / -929 /
!!      data imach(13) / 1070 /
!!      data imach(14) /   94 /
!!      data imach(15) / -929 /
!!      data imach(16) / 1069 /
!!
!!  CONVEX C-1.
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   32 /
!!      data imach( 6) /    4 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -128 /
!!      data imach(13) /  127 /
!!      data imach(14) /   53 /
!!      data imach(15) /-1024 /
!!      data imach(16) / 1023 /
!!
!!  CONVEX C-120 (native mode) without -R8 option
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    0 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   32 /
!!      data imach( 6) /    4 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -127 /
!!      data imach(13) /  127 /
!!      data imach(14) /   53 /
!!      data imach(15) / -1023 /
!!      data imach(16) /  1023 /
!!
!!  CONVEX C-120 (native mode) with -R8 option
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /     0 /
!!      data imach( 4) /     6 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     4 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /     2 /
!!      data imach(11) /    53 /
!!      data imach(12) / -1023 /
!!      data imach(13) /  1023 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1023 /
!!      data imach(16) /  1023 /
!!
!!  CONVEX C-120 (IEEE mode) without -R8 option
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    0 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   32 /
!!      data imach( 6) /    4 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -125 /
!!      data imach(13) /  128 /
!!      data imach(14) /   53 /
!!      data imach(15) / -1021 /
!!      data imach(16) /  1024 /
!!
!!  CONVEX C-120 (IEEE mode) with -R8 option
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /     0 /
!!      data imach( 4) /     6 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     4 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /     2 /
!!      data imach(11) /    53 /
!!      data imach(12) / -1021 /
!!      data imach(13) /  1024 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1021 /
!!      data imach(16) /  1024 /
!!
!!  CRAY 1, 2, XMP and YMP.
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /   102 /
!!      data imach( 4) /     6 /
!!      data imach( 5) /    64 /
!!      data imach( 6) /     8 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    63 /
!!      data imach( 9) /  777777777777777777777B /
!!      data imach(10) /     2 /
!!      data imach(11) /    47 /
!!      data imach(12) / -8189 /
!!      data imach(13) /  8190 /
!!      data imach(14) /    94 /
!!      data imach(15) / -8099 /
!!      data imach(16) /  8190 /
!!
!!  DATA GENERAL ECLIPSE S/200.
!!
!!      data imach( 1) /   11 /
!!      data imach( 2) /   12 /
!!      data imach( 3) /    8 /
!!      data imach( 4) /   10 /
!!      data imach( 5) /   16 /
!!      data imach( 6) /    2 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   15 /
!!      data imach( 9) /32767 /
!!      data imach(10) /   16 /
!!      data imach(11) /    6 /
!!      data imach(12) /  -64 /
!!      data imach(13) /   63 /
!!      data imach(14) /   14 /
!!      data imach(15) /  -64 /
!!      data imach(16) /   63 /
!!
!!  ELXSI 6400
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /     6 /
!!      data imach( 4) /     6 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     4 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    32 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /     2 /
!!      data imach(11) /    24 /
!!      data imach(12) /  -126 /
!!      data imach(13) /   127 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1022 /
!!      data imach(16) /  1023 /
!!
!!  HARRIS 220
!!
!!      data imach( 1) /       5 /
!!      data imach( 2) /       6 /
!!      data imach( 3) /       0 /
!!      data imach( 4) /       6 /
!!      data imach( 5) /      24 /
!!      data imach( 6) /       3 /
!!      data imach( 7) /       2 /
!!      data imach( 8) /      23 /
!!      data imach( 9) / 8388607 /
!!      data imach(10) /       2 /
!!      data imach(11) /      23 /
!!      data imach(12) /    -127 /
!!      data imach(13) /     127 /
!!      data imach(14) /      38 /
!!      data imach(15) /    -127 /
!!      data imach(16) /     127 /
!!
!!  HARRIS SLASH 6 and SLASH 7.
!!
!!      data imach( 1) /       5 /
!!      data imach( 2) /       6 /
!!      data imach( 3) /       0 /
!!      data imach( 4) /       6 /
!!      data imach( 5) /      24 /
!!      data imach( 6) /       3 /
!!      data imach( 7) /       2 /
!!      data imach( 8) /      23 /
!!      data imach( 9) / 8388607 /
!!      data imach(10) /       2 /
!!      data imach(11) /      23 /
!!      data imach(12) /    -127 /
!!      data imach(13) /     127 /
!!      data imach(14) /      38 /
!!      data imach(15) /    -127 /
!!      data imach(16) /     127 /
!!
!!  HONEYWELL DPS 8/70 and 600/6000 series.
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /   43 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   36 /
!!      data imach( 6) /    4 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   35 /
!!      data imach( 9) / O377777777777 /
!!      data imach(10) /    2 /
!!      data imach(11) /   27 /
!!      data imach(12) / -127 /
!!      data imach(13) /  127 /
!!      data imach(14) /   63 /
!!      data imach(15) / -127 /
!!      data imach(16) /  127 /
!!
!!  HP 2100, 3 word double precision option with FTN4
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    4 /
!!      data imach( 4) /    1 /
!!      data imach( 5) /   16 /
!!      data imach( 6) /    2 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   15 /
!!      data imach( 9) / 32767 /
!!      data imach(10) /    2 /
!!      data imach(11) /   23 /
!!      data imach(12) / -128 /
!!      data imach(13) /  127 /
!!      data imach(14) /   39 /
!!      data imach(15) / -128 /
!!      data imach(16) /  127 /
!!
!!  HP 2100, 4 word double precision option with FTN4
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    4 /
!!      data imach( 4) /    1 /
!!      data imach( 5) /   16 /
!!      data imach( 6) /    2 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   15 /
!!      data imach( 9) / 32767 /
!!      data imach(10) /    2 /
!!      data imach(11) /   23 /
!!      data imach(12) / -128 /
!!      data imach(13) /  127 /
!!      data imach(14) /   55 /
!!      data imach(15) / -128 /
!!      data imach(16) /  127 /
!!
!!  HP 9000
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /     6 /
!!      data imach( 4) /     7 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     4 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    32 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /     2 /
!!      data imach(11) /    24 /
!!      data imach(12) /  -126 /
!!      data imach(13) /   127 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1015 /
!!      data imach(16) /  1017 /
!!
!!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!!  and PERKIN ELMER (INTERDATA) 3230.
!!
!!      data imach( 1) /   5 /
!!      data imach( 2) /   6 /
!!      data imach( 3) /   7 /
!!      data imach( 4) /   6 /
!!      data imach( 5) /  32 /
!!      data imach( 6) /   4 /
!!      data imach( 7) /   2 /
!!      data imach( 8) /  31 /
!!      data imach( 9) / Z7FFFFFFF /
!!      data imach(10) /  16 /
!!      data imach(11) /   6 /
!!      data imach(12) / -64 /
!!      data imach(13) /  63 /
!!      data imach(14) /  14 /
!!      data imach(15) / -64 /
!!      data imach(16) /  63 /
!!
!!  IBM PC - Microsoft FORTRAN
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /     6 /
!!      data imach( 4) /     0 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     4 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /     2 /
!!      data imach(11) /    24 /
!!      data imach(12) /  -126 /
!!      data imach(13) /   127 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1022 /
!!      data imach(16) /  1023 /
!!
!!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!!
!!      data imach( 1) /     4 /
!!      data imach( 2) /     7 /
!!      data imach( 3) /     7 /
!!      data imach( 4) /     0 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     4 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /     2 /
!!      data imach(11) /    24 /
!!      data imach(12) /  -126 /
!!      data imach(13) /   127 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1022 /
!!      data imach(16) /  1023 /
!!
!!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!!  For the INTERDATA FORTRAN VII compiler, replace the Z's specifying hex
!!  constants with Y's.
!!
!!      data imach( 1) /   5 /
!!      data imach( 2) /   6 /
!!      data imach( 3) /   6 /
!!      data imach( 4) /   6 /
!!      data imach( 5) /  32 /
!!      data imach( 6) /   4 /
!!      data imach( 7) /   2 /
!!      data imach( 8) /  31 /
!!      data imach( 9) / Z'7FFFFFFF' /
!!      data imach(10) /  16 /
!!      data imach(11) /   6 /
!!      data imach(12) / -64 /
!!      data imach(13) /  62 /
!!      data imach(14) /  14 /
!!      data imach(15) / -64 /
!!      data imach(16) /  62 /
!!
!!  PDP-10 (KA processor).
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   36 /
!!      data imach( 6) /    5 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   35 /
!!      data imach( 9) / "377777777777 /
!!      data imach(10) /    2 /
!!      data imach(11) /   27 /
!!      data imach(12) / -128 /
!!      data imach(13) /  127 /
!!      data imach(14) /   54 /
!!      data imach(15) / -101 /
!!      data imach(16) /  127 /
!!
!!  PDP-10 (KI processor).
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   36 /
!!      data imach( 6) /    5 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   35 /
!!      data imach( 9) / "377777777777 /
!!      data imach(10) /    2 /
!!      data imach(11) /   27 /
!!      data imach(12) / -128 /
!!      data imach(13) /  127 /
!!      data imach(14) /   62 /
!!      data imach(15) / -128 /
!!      data imach(16) /  127 /
!!
!!  PDP-11 FORTRANS supporting 32-bit integer arithmetic.
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   32 /
!!      data imach( 6) /    4 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -127 /
!!      data imach(13) /  127 /
!!      data imach(14) /   56 /
!!      data imach(15) / -127 /
!!      data imach(16) /  127 /
!!
!!  PDP-11 FORTRANS supporting 16-bit integer arithmetic.
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   16 /
!!      data imach( 6) /    2 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   15 /
!!      data imach( 9) / 32767 /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -127 /
!!      data imach(13) /  127 /
!!      data imach(14) /   56 /
!!      data imach(15) / -127 /
!!      data imach(16) /  127 /
!!
!!  PRIME 50 series systems with 32-bit integers and 64V MODE instructions,
!!  supplied by Igor Bray.
!!
!!      data imach( 1) /            1 /
!!      data imach( 2) /            1 /
!!      data imach( 3) /            2 /
!!      data imach( 4) /            1 /
!!      data imach( 5) /           32 /
!!      data imach( 6) /            4 /
!!      data imach( 7) /            2 /
!!      data imach( 8) /           31 /
!!      data imach( 9) / :17777777777 /
!!      data imach(10) /            2 /
!!      data imach(11) /           23 /
!!      data imach(12) /         -127 /
!!      data imach(13) /         +127 /
!!      data imach(14) /           47 /
!!      data imach(15) /       -32895 /
!!      data imach(16) /       +32637 /
!!
!!  SEQUENT BALANCE 8000.
!!
!!      data imach( 1) /     0 /
!!      data imach( 2) /     0 /
!!      data imach( 3) /     7 /
!!      data imach( 4) /     0 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     1 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    31 /
!!      data imach( 9) /  2147483647 /
!!      data imach(10) /     2 /
!!      data imach(11) /    24 /
!!      data imach(12) /  -125 /
!!      data imach(13) /   128 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1021 /
!!      data imach(16) /  1024 /
!!
!!  SUN Microsystems UNIX F77 compiler.
!!
!!      data imach( 1) /     5 /
!!      data imach( 2) /     6 /
!!      data imach( 3) /     6 /
!!      data imach( 4) /     0 /
!!      data imach( 5) /    32 /
!!      data imach( 6) /     4 /
!!      data imach( 7) /     2 /
!!      data imach( 8) /    32 /
!!      data imach( 9) /2147483647/
!!      data imach(10) /     2 /
!!      data imach(11) /    24 /
!!      data imach(12) /  -126 /
!!      data imach(13) /   128 /
!!      data imach(14) /    53 /
!!      data imach(15) / -1022 /
!!      data imach(16) /  1024 /
!!
!!  SUN 3 (68881 or FPA)
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    6 /
!!      data imach( 4) /    0 /
!!      data imach( 5) /   32 /
!!      data imach( 6) /    4 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -125 /
!!      data imach(13) /  128 /
!!      data imach(14) /   53 /
!!      data imach(15) / -1021 /
!!      data imach(16) /  1024 /
!!
!!  UNIVAC 1100 series.
!!  Note that the punch unit, I1MACH(3), has been set to 7, which is appropriate
!!  for the UNIVAC-FOR system.  If you have the UNIVAC-FTN system, set it to 1
!!  instead.
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   36 /
!!      data imach( 6) /    6 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   35 /
!!      data imach( 9) / O377777777777 /
!!      data imach(10) /    2 /
!!      data imach(11) /   27 /
!!      data imach(12) / -128 /
!!      data imach(13) /  127 /
!!      data imach(14) /   60 /
!!      data imach(15) /-1024 /
!!      data imach(16) / 1023 /
!!
!!  VAX.
!!
!!      data imach( 1) /    5 /
!!      data imach( 2) /    6 /
!!      data imach( 3) /    7 /
!!      data imach( 4) /    6 /
!!      data imach( 5) /   32 /
!!      data imach( 6) /    4 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   31 /
!!      data imach( 9) / 2147483647 /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -127 /
!!      data imach(13) /  127 /
!!      data imach(14) /   56 /
!!      data imach(15) / -127 /
!!      data imach(16) /  127 /
!!
!!  Z80 microprocessor.
!!
!!      data imach( 1) /    1 /
!!      data imach( 2) /    1 /
!!      data imach( 3) /    0 /
!!      data imach( 4) /    1 /
!!      data imach( 5) /   16 /
!!      data imach( 6) /    2 /
!!      data imach( 7) /    2 /
!!      data imach( 8) /   15 /
!!      data imach( 9) / 32767 /
!!      data imach(10) /    2 /
!!      data imach(11) /   24 /
!!      data imach(12) / -127 /
!!      data imach(13) /  127 /
!!      data imach(14) /   56 /
!!      data imach(15) / -127 /
!!      data imach(16) /  127 /
!!
!  if ( i < 1 .or. i > 16 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'I1MACH - Fatal error!'
!    write ( *, '(a,i6)' ) 'I is out of bounds:',i
!    i1mach = 0
!    stop
!  else
!    i1mach = imach(i)
!  end if
!
!  return
!end
!subroutine iratex ( func, a, b, epsin, epsout, result, ind )
!!
!!***********************************************************************
!!
!!! IRATEX estimates the integral of a function.
!!
!!
!!  Discussion:
!!
!!    IRATEX estimates the integral from A to B of F(X), using the
!!    trapezoidal rule for several stepsizes H.
!!
!!    Then a rational function of H*H is fitted to these results, and 
!!    an estimate of the integral is computed by extrapolating the 
!!    results to a stepsize of zero.
!!
!!  Reference:
!!
!!    R Bulirsch and J Stoer,
!!    Fehlerabschaetzungen und Extrapolation mit rationaled Funktionen
!!      bei Verfahren vom Richardson-Typus,
!!    (Error estimates and extrapolation with rational functions
!!      in processes of Richardson type),
!!    Numerische Mathematik,
!!    Volume 6 (1964), pages 413-427.
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real, external FUNC.
!!    FUNC is the name of the function to be
!!    integrated.  The user must declare the name an external
!!    parameter in the calling program, pass the name of the
!!    function in FUNC, and write a function of the form
!!
!!      FUNCTION FUNC(X)
!!
!!    which evaluates the function at the point X.
!!
!!    Input, real A, the lower limit of integration.
!!
!!    Input, real B, the upper limit of integration.
!!
!!    Input, real EPSIN, the requested relative error tolerance.
!!
!!    Output, real EPSOUT, an estimate of the error in the integration.
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!!    Output, integer IND, error return flag.
!!    IND = 0 if the requested accuracy was not achieved,
!!    IND = 1 if the accuracy was achieved.
!!
!  implicit none
!!
!  real a
!  real arg
!  real b
!  real ba
!  logical bo
!  logical bu
!  real c
!  real d(6)
!  real d1
!  real ddt
!  real den
!  real dt(7)
!  real e
!  real ent
!  real epsin
!  real epsout
!  real, external :: func
!  real gr
!  real hm
!  integer i
!  integer ind
!  integer m
!  integer mr
!  integer n
!  integer np1
!  logical odd
!  real rnderr
!  real result
!  real sm
!  real t
!  real t1
!  real t2
!  real t2a
!  real ta
!  real tab
!  real tb
!  real tnt
!  real v
!  real w
!!
!  if ( a == b ) then
!    result = 0.0E+00
!    return
!  end if
! 
!  rnderr = epsilon ( 1.0E+00 )
!  epsin = max ( epsin, 8.0 * rnderr )
!  ind = 0
!  n = 2
!  np1 = 3
!  ba = b-a
!  t1 = 0.0E+00
!  gr = 0.0E+00
!  sm = 0.0E+00
!  t2a = 0.5E+00 * ( func ( a ) + func ( b ) )
!  t2 = t2a
!  tb = abs ( t2a )
!  c = t2 * ba
!
!  dt(1) = c
!  dt(2:7) = 0.0E+00
! 
!  odd = .true.
!  bu = .false.
! 
!  do m = 1, 15
! 
!    bo = m>=7
!    hm = ba / real(n)
!!
!!  N+1 is odd.
!!
!    if ( odd ) then
! 
!      do i = 1, n, 2
!        arg = a + real(i) * hm
!        w = func(arg)
!        t2 = t2 + w
!        tb = abs ( w ) + tb
!      end do
! 
!      ent = t2
!      tab = tb * abs ( hm )
!      d(1) = 16.0E+00 / 9.0E+00
!      d(3) = 4.0E+00 * d(1)
!      d(5) = 4.0E+00 * d(3)
!!
!!  N+1 is even.
!!
!    else
! 
!      do i = 1, n, 6
!        w = real(i) * hm
!        t1 = t1 + func(a+w)+func(b-w)
!      end do
! 
!      ent = t1+t2a
!      t2a = t2
!      d(1) = 2.25E+00
!      d(3) = 9.0E+00
!      d(5) = 36.0E+00
! 
!    end if
! 
!    ddt = dt(1)
!    t = ent*hm
!    dt(1) = t
!    ent = t
! 
!    if ( m < 7 ) then
!      mr = m
!      w = real ( n * n )
!      d(m) = w
!    else
!      mr = 6
!      d(6) = 64.0E+00
!      w = 144.0E+00
!    end if
! 
!    do i = 1, mr
! 
!      d1 = d(i) * ddt
!      den = d1 - ent
!      e = ent - ddt
!      tnt = ent
!      v = 0.0E+00
!      ent = 0.0E+00
! 
!      if ( abs ( den ) >= epsin ) then
!        e = e / den
!        v = tnt * e
!        ent = d1 * e
!        t = v + t
!        ddt = dt(i+1)
!      end if
! 
!      dt(i+1) = v
! 
!    end do
! 
!    ta = c
!    c = t
!    result = c
!    if ( .not. bo ) then
!      t = t-v
!    end if
!
!    v = t-ta
!    t = v+t
!    epsout = abs ( v )
! 
!    if ( ta < t ) then
!      d1 = ta
!      ta = t
!      t = d1
!    end if
! 
!    bo = bo .or. ( ta < gr .and. t > sm )
! 
!    if ( bu .and. bo .and. epsout < epsin * w * tab ) then
!      go to 140
!    end if
! 
!    gr = ta + epsin
!    sm = t - epsin
!    odd = .not. odd
!    n = np1
!    np1 = n+1
!    bu = bo
!    d(2) = 4.0E+00
!    d(4) = 16.0E+00
! 
!  end do
! 
!  bo = .false.
! 
!  140 continue
! 
!  v = rnderr*tab
! 
!  epsout = max ( epsout, v )
!
!  if (bo) ind = 1
!
!  return
!end
!function pi ( )
!!
!!*******************************************************************************
!!
!!! PI returns the value of pi.
!!
!!
!!  Modified:
!!
!!    04 December 1998
!!
!!  Author:
!!
!!    John Burkardt
!!
!!  Parameters:
!!
!!    Output, real PI, the value of pi.
!!
!  implicit none
!!
!  real pi
!!
!  pi = 3.14159265358979323846264338327950288419716939937510E+00
!
!  return
!end
!subroutine plint ( ftab, xtab, ntab, a, b, result )
!!
!!***********************************************************************
!!
!!! PLINT approximates the integral of unequally spaced data.
!!
!!
!!  Discussion:
!!
!!    The method uses piecewise linear interpolation.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real FTAB(NTAB), the function values, FTAB(I) = F(XTAB(I)).
!!
!!    Input, real XTAB(NTAB), the abscissas at which the
!!    function values are given.  The XTAB's must be distinct
!!    and in ascending order.
!!
!!    Input, integer NTAB, the number of entries in FTAB and
!!    XTAB.  NTAB must be at least 2.
!!
!!    Input, real A, the lower limit of integration.  A should
!!    be, but need not be, near one endpoint of the interval
!!    (X(1), X(NTAB)).
!!
!!    Input, real B, the upper limit of integration.  B should
!!    be, but need not be, near one endpoint of the interval
!!    (X(1), X(NTAB)).
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!  implicit none
!!
!  integer ntab
!!
!  real a
!  real b
!  real fa
!  real fb
!  real ftab(ntab)
!  integer i
!  integer ihi
!  integer ilo
!  integer ind
!  real result
!  real slope
!  real syl
!  real xtab(ntab)
!!
!  if ( a == b ) then
!    result = 0.0E+00
!    return
!  end if
! 
!  if ( ntab < 2 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'PLINT - Fatal error!'
!    write ( *, '(a,i6)' ) '  NTAB < 2, NTAB = ',ntab
!    stop
!  end if
! 
!  do i = 2, ntab
!    if ( xtab(i) <= xtab(i-1) ) then
!      write ( *, '(a)' ) ' '
!      write ( *, '(a)' ) 'PLINT - Fatal error!'
!      write ( *, '(a)' ) '  Nodes not in strict increasing order.'
!      write ( *, '(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
!      write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
!      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
!      stop
!    end if
!  end do
!!
!!  If A > B, temporarily switch A and B, and store sign.
!!
!  if ( a > b ) then
!    syl = b
!    b = a
!    a = syl
!    ind = -1
!  else
!    syl = a
!    ind = 1
!  end if
!!
!!  Find ILO and IHI so that A <= XTAB(ILO) <= XTAB(IHI) <= B
!!  with the possible exception that A and B may be in the same
!!  interval, or completely to the right or left of the XTAB's.
!!
!  ilo = 1
!  ihi = ntab
!  do i = 1, ntab
!    if ( a <= xtab(i) ) then
!      exit
!    end if
!    ilo = ilo+1
!  end do
!  
!  do i = 1, ntab
!    if ( b >= xtab(i) ) then
!      exit
!    end if
!    ihi = ihi-1
!  end do
!!
!!  Treat special cases where A, B lie both to left or both to right
!!  of XTAB interval, or inbetween same pair of XTAB's.
!!
!  if ( ihi == 0 ) then
!    slope = (ftab(2)-ftab(1))/(xtab(2)-xtab(1))
!    fa = ftab(1) + slope*(a-xtab(1))
!    fb = ftab(1) + slope*(b-xtab(1))
!    result = 0.5 * (b-a) * (fa+fb)
!    go to 110
!  else if ( ilo == ntab+1 ) then
!    slope = (ftab(ntab)-ftab(ntab-1))/(xtab(ntab)-xtab(ntab-1))
!    fa = ftab(ntab-1)+slope*(a-xtab(ntab-1))
!    fb = ftab(ntab-1)+slope*(b-xtab(ntab-1))
!    result = 0.5 * (b-a) * (fa+fb)
!    go to 110
!  else if ( ihi+1 == ilo ) then
!    slope = (ftab(ilo)-ftab(ihi))/(xtab(ilo)-xtab(ihi))
!    fa = ftab(ihi)+slope*(a-xtab(ihi))
!    fb = ftab(ihi)+slope*(b-xtab(ihi))
!    result = 0.5 * (b-a) * (fa+fb)
!    go to 110
!  end if
!!
!!  Carry out approximate integration.  We know that ILO is no greater
!!  than IHI-1, but equality is possible; A and B may be on either side
!!  of a single XTAB(I).  That's OK, then the loop below won't be executed
!!  at all.
!!
!  result = 0.0E+00
!  do i = ilo, ihi-1
!    result = result + 0.5 * (xtab(i+1)-xtab(i))*(ftab(i)+ftab(i+1))
!  end do
!!
!!  Add contribution from A-ILO and IHI-B.
!!  Still have to watch out if ILO = 1 or IHI=NTAB...
!!
!  if ( ilo == 1 ) then
!    slope = (ftab(2)-ftab(1)) / (xtab(2)-xtab(1))
!    fa = ftab(1) + slope*(a-xtab(1))
!    result = result + 0.5 * (xtab(ilo)-a)*(fa+ftab(ilo))
!  else
!    slope = (ftab(ilo)-ftab(ilo-1)) / (xtab(ilo)-xtab(ilo-1))
!    fa = ftab(ilo-1) + slope*(a-xtab(ilo-1))
!    result = result + 0.5 * (xtab(ilo)-a)*(fa+ftab(ilo))
!  end if
! 
!  if ( ihi == ntab ) then
!    slope = (ftab(ntab)-ftab(ntab-1)) / (xtab(ntab)-xtab(ntab-1))
!    fb = ftab(ntab-1) + slope*(b-xtab(ntab-1))
!    result = result + 0.5*(b-xtab(ntab))*(fb+ftab(ntab))
!  else
!    slope = (ftab(ihi+1)-ftab(ihi)) / (xtab(ihi+1)-xtab(ihi))
!    fb = ftab(ihi) + slope*(b-xtab(ihi))
!    result = result + 0.5*(b-xtab(ihi))*(fb+ftab(ihi))
!  end if
!!
!!  Restore original values of A and B, reverse sign of integral
!!  because of earlier switch.
!!
!110   continue
! 
!  if ( ind /= 1 ) then
!    ind = 1
!    syl = b
!    b = a
!    a = syl
!    result = -result
!  end if
! 
!  return
!end
!subroutine qnc79 ( func, a, b, err, result, ierr, k )
!!
!!***********************************************************************
!!
!!! QNC79 approximates the integral of F(X) using Newton-Cotes quadrature.
!!
!!
!!  Discussion:
!!
!!    QNC79 is a general purpose program for evaluation of one
!!    dimensional integrals  of user defined functions.  QNC79 will
!!    pick its own points for evaluation of the integrand and these
!!    will vary from problem to problem.
!!
!!    Thus QNC79 is not designed to integrate over data sets.
!!
!!    Moderately smooth integrands will be integrated efficiently
!!    and reliably.  For problems with strong singularities,
!!    oscillations etc., the user may wish to use more sophisticated
!!    routines such as those in QUADPACK.
!!
!!    One measure of the reliability of QNC79 is the output parameter
!!    K, giving the number of integrand evaluations that were needed.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real, external FUNC, the name of the function to be
!!    integrated.  The user must declare the name an external
!!    parameter in the calling program, pass the name of the
!!    function in FUNC, and write a function of the form
!!
!!      FUNCTION FUNC(X)
!!
!!    which evaluates the function at the point X.
!!
!!    Input, real A, lower limit of integral.
!!
!!    Input, real B, upper limit of integral.
!!
!!    Input, real ERR, is a requested error tolerance.
!!    Normally pick a value, 0 .LT. ERR .LT. 1.E-3.
!!
!!    Output, real RESULT, computed value of the integral.
!!    Hopefully RESULT is accurate to within ERR times the
!!    integral of ABS ( FUNC(X) ).
!!
!!    Output, integer IERR, a status code
!!     1 RESULT most likely meets requested error tolerance.
!!    -1 A and B are too nearly equal to allow normal integration.
!!     2 RESULT probably does not meet requested error tolerance.
!!
!!    Output, integer K, the number of function evaluations
!!    actually used to do the integration.
!!    A value of K .GT. 1000 indicates a difficult problem.
!!    Other programs may be more efficient.
!!    QNC79 will gracefully give up if K exceeds 2000.
!!
!  implicit none
!!
!  integer, parameter :: kmx = 2000
!!
!  real a
!  real aa(40)
!  real ae
!  real area
!  real b
!  real bank
!  real blocal
!  real c
!  real ce
!  real ee
!  real ef
!  real eps
!  real err
!  real f(13)
!  real f1(40)
!  real f2(40)
!  real f3(40)
!  real f4(40)
!  real f5(40)
!  real f6(40)
!  real f7(40)
!  real, external :: func
!  real hh(40)
!  integer i
!  integer i1mach
!  integer, save :: icall = 0
!  integer ierr
!  integer k
!  integer, save :: kml = 7
!  integer l
!  integer lmn
!  integer lmx
!  integer lr(40)
!  integer nbits
!  integer nib
!  integer, save :: nlmn = 2
!  integer nlmx
!  real q13
!  real q7
!  real q7l
!  real q7r(40)
!  real r1mach
!  real result
!  real test
!  real tol
!  real vl(40)
!  real vr
!  real w1
!  real w2
!  real w3
!  real w4
!!
!  if ( a == b ) then
!    result = 0.0E+00
!    return
!  end if
! 
!  if ( icall /= 0 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'QNC79 - Fatal error!'
!    write ( *, '(a)' ) '  QNC79 was called recursively!'
!    stop
!  end if
! 
!  icall = 1
!  w1 = 41.0E+00 / 140.0E+00
!  w2  = 216.0 / 140.0E+00
!  w3 = 27.0 / 140.0E+00
!  w4  = 272.0 / 140.0E+00
!  nbits = int ( r1mach(5) * real(i1mach(11)) / 0.30102000 )
!  nlmx = min ( 40, (nbits*4)/5 )
!  result = 0.0E+00
!  ierr = 1
!  ce = 0.0E+00
!  lmx = nlmx
!  lmn = nlmn
!  if ( b == 0.0E+00 ) go to 3
!  if ( sign ( 1.0E+00, b ) * a <= 0.0E+00 ) go to 3
!  c = abs ( 1.0E+00 - a / b )
!
!  if ( c > 0.1E+00 ) then
!    go to 3
!  end if
! 
!  if ( c <= 0.0E+00 ) then
!    ierr  = -1
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'QNC79 - Fatal error!'
!    write ( *, '(a)' ) '  A and B are too close.'
!    stop
!  end if
! 
!  nib = int ( 0.5E+00 - log(c) / log(2.0E+00) )
!  lmx = min ( nlmx, nbits-nib-4)
!
!  if ( lmx < 2 ) then
!    go to 32
!  end if
!
!  lmn = min(lmn,lmx)
! 
!3 continue
! 
!  tol = max ( abs ( err ), 2.0E+00**(5-nbits) )
!  if ( err == 0.0E+00 ) tol = sqrt ( epsilon ( tol ) )
!  eps = tol
!  hh(1) = (b-a) / 12.0E+00
!  aa(1) = a
!  lr(1) = 1
! 
!  do i = 1, 11, 2
!    f(i) = func(a+real(i-1)*hh(1))
!  end do
! 
!  blocal = b
!  f(13) = func(blocal)
!  k = 7
!  l = 1
!  area = 0.0E+00
!  q7 = 0.0E+00
!  ef = 256.0E+00 / 255.0E+00
!  bank = 0.0E+00
!!
!!  Compute refined estimates, estimate the error, etc.
!!
!5 continue
! 
!  do i = 2, 12, 2
!    f(i) = func ( aa(l)+real(i-1)*hh(l) )
!  end do
! 
!  k = k+6
!!
!!  Compute left and right half estimates.
!!
!  q7l = hh(l)*( ( w1*(f(1)+f(7)) + w2*(f(2)+f(6)) ) &
!              + ( w3*(f(3)+f(5)) + w4*f(4) ) )
!
!  q7r(l) = hh(l)*( ( w1*(f(7)+f(13)) + w2*(f(8)+f(12)) ) &
!                + ( w3*(f(9)+f(11)) + w4*f(10) ) )
!!
!!  Update estimate of integral of absolute value.
!!
!  area = area + ( abs ( q7l ) + abs ( q7r(l) ) - abs ( q7 ) )
!!
!!  Do not bother to test convergence before minimum refinement level.
!!
!  if (l<lmn) go to 11
!!
!!  Estimate the error in new value for whole interval, Q13.
!!
!  q13 = q7l+q7r(l)
!  ee = abs ( q7 - q13 ) * ef
!!
!!  Compute nominal allowed error.
!!
!  ae = eps*area
!!
!!  Borrow from bank account, but not too much.
!!
!  test = min(ae+0.8E+00*bank, 10.0E+00*ae)
!!
!!  Don't ask for excessive accuracy.
!!
!  test = max ( test, tol * abs ( q13 ), 0.00003E+00 * tol * area )
!!
!!  Now, did this interval pass or not?
!!
!  if ( ee <= test )go to 8
!  go to 10
!!
!!  Have hit max refinement level - penalize the cumulative error.
!!
!    7 continue
! 
!  ce = ce+(q7-q13)
!  go to 9
!!
!!  On good intervals accumulate the theoretical estimate.
!!
!    8 ce = ce+(q7-q13) / 255.0E+00
!!
!!  Update the bank account.  Don't go into debt.
!!
!9 continue
!
!  bank = bank + (ae-ee)
!  if ( bank < 0.0E+00 ) bank = 0.0E+00
!!
!!  Did we just finish a left half or a right half?
!!
!  if ( lr(l) <= 0 ) go to 15
!  go to 20
!!
!!  Consider the left half of the next deeper level.
!!
!10 continue
! 
!  if ( k > kmx ) lmx = min(kml,lmx)
!  if ( l >= lmx ) go to 7
!
!11 continue
!
!  l = l+1
!  eps = eps * 0.5E+00
!  if ( l <= 17 ) ef = ef/sqrt(2.0E+00)
!  hh(l) = hh(l-1)*0.5E+00
!  lr(l) = -1
!  aa(l) = aa(l-1)
!  q7 = q7l
!  f1(l) = f(7)
!  f2(l) = f(8)
!  f3(l) = f(9)
!  f4(l) = f(10)
!  f5(l) = f(11)
!  f6(l) = f(12)
!  f7(l) = f(13)
!  f(13) = f(7)
!  f(11) = f(6)
!  f(9)  = f(5)
!  f(7)  = f(4)
!  f(5)  = f(3)
!  f(3)  = f(2)
!  go to 5
!!
!!  Proceed to right half at this level
!!
!   15 vl(l) = q13
!   16 q7 = q7r(l-1)
!  lr(l) = 1
!  aa(l) = aa(l)+12.0E+00 * hh(l)
!  f(1)  = f1(l)
!  f(3)  = f2(l)
!  f(5)  = f3(l)
!  f(7)  = f4(l)
!  f(9)  = f5(l)
!  f(11) = f6(l)
!  f(13) = f7(l)
!  go to 5
!!
!!  Left and right halves are done, so go back up a level
!!
!   20 vr = q13
!   22 if ( l <= 1 ) go to 30
!
!  if ( l <= 17 ) then
!    ef = ef * sqrt ( 2.0E+00 )
!  end if
!
!  eps = eps * 2.0E+00
!  l = l-1
! 
!  if ( lr(l) <= 0 ) then
!    vl(l) = vl(l+1)+vr
!    go to 16
!  else
!    vr = vl(l+1)+vr
!    go to 22
!  end if
! 
!   30 continue
! 
!  if ( abs ( ce ) > 2.0E+00 * tol * area ) then
!    ierr = 2
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'QNC79 - Warning!'
!    write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
!  end if
! 
!32    continue
! 
!  result = vr
! 
!  icall = 0
! 
!  return
!end
!subroutine quad ( func, a, b, abserr, relerr, nleast, nmost, work, &
!  result )
!!
!!***********************************************************************
!!
!!! QUAD approximates the integral of F(X) by Romberg integration.
!!
!!
!!  Discussion:
!!
!!    The integration is repeated until convergence of the results, 
!!    or the maximum number of steps is taken.  The Romberg method 
!!    uses the trapezoidal rule, subinterval halving, and Richardson 
!!    extrapolation.
!!
!!    Convergence is declared if either of the following occurs:
!!
!!      ABS ( RESULT - INTEGRAL ) < ABSERR
!!
!!    or
!!
!!      RESULT = integral from A to B (1+Y(X))*FUNC(X)DX for some
!!      function Y with ABS ( Y(X) ) < RELERR+RNDERR  with RNDERR the
!!      machine rounding error.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!    C F Dunkl,
!!    Romberg quadrature to prescribed accuracy,
!!    SHARE file number 7090-1481 TYQUAD
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real, external FUNC, the name of the function to be integrated.
!!    The user must declare the name an external parameter in the
!!    calling program, pass the name of the function in FUNC,
!!    and write a function of the form FUNCTION FUNC(X) which
!!    evaluates the function at the point X.
!!
!!    Input, real A, the lower limit of integration.
!!
!!    Input, real B, the upper limit of integration.
!!
!!    Input, real ABSERR, the absolute error tolerance.
!!
!!    Input, real RELERR, the relative error tolerance.
!!
!!    Input, integer NLEAST, the least number of times the integration
!!    is to be carried out before the convergence test is made.
!!    A value 3 <= NLEAST <= 15 is suggested.
!!
!!    Input, integer NMOST, the most number of times the
!!    integration is to be carried out.
!!
!!    Workspace, real WORK(NMOST+1).
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!  implicit none
!!
!  integer nmost
!!
!  real a
!  real abserr
!  real b
!  real f
!  real fcna
!  real fcnb
!  real fcnxi
!  real, external :: func
!  real h
!  integer i
!  integer j
!  integer k
!  integer nleast
!  integer nx
!  real qx1
!  real qx2
!  real relerr
!  real result
!  real rnderr
!  real sum1
!  real sumabs
!  real t
!  real tabs
!  real work(nmost+1)
!  real x
!!
!  if ( a == b ) then
!    result = 0.0E+00
!    return
!  end if
! 
!  rnderr = epsilon ( 1.0E+00 )
!  qx1 = 0
!  qx2 = 0
!  h = b-a
!  fcna = func(a)
!  fcnb = func(b)
!  tabs = abs ( h ) * ( abs ( fcna ) + abs ( fcnb ) ) / 2.0E+00
!  t = h*(fcna+fcnb) / 2.0E+00
!  nx = 1
! 
!  do i = 1, nmost
! 
!    h = 0.5*h
! 
!    sum1 = 0.0E+00
!    sumabs = 0.0E+00
!    do j = 1, nx
!      fcnxi = func(a+h*(2*j-1))
!      sumabs = sumabs + abs ( fcnxi )
!      sum1 = sum1+fcnxi
!    end do
! 
!    t = 0.5E+00 * t + h * sum1
!    tabs = tabs/2.0E+00 + abs ( h ) * sumabs
!    work(i) = 2.0E+00*(t+h*sum1)/3.0E+00
!
!    if ( i > 1 ) then
!!
!!  Construct difference table for Richardson extrapolation.
!!
!      f = 4.0E+00
!      do j = 2, i
!        k = i+1-j
!        f = f*4.0E+00
!        work(k) = work(k+1)+(work(k+1)-work(k) ) / ( f - 1.0E+00 )
!      end do
!!
!!  Perform acceptance check.
!!
!      if ( i >= nleast ) then
!        x = abs ( work(1)-qx2 ) + abs ( qx2 - qx1 )
!        if ( x <= 3.0E+00 * tabs * ( abs ( relerr ) + rnderr ) ) go to 10
!        if ( x <= 3.0E+00 * abs ( abserr ) ) go to 10
!      end if
!!
!!  Save old result, perform bisection, repeat.
!!
!      qx1 = qx2
!    end if
! 
!    qx2 = work(1)
!    nx = nx*2
! 
!  end do
!!
!!  Accept result.
!!
!10    continue
!  result = work(1)
!  return
!end
!subroutine r_swap ( x, y )
!!
!!*******************************************************************************
!!
!!! R_SWAP swaps two real values.
!!
!!
!!  Modified:
!!
!!    01 May 2000
!!
!!  Author:
!!
!!    John Burkardt
!!
!!  Parameters:
!!
!!    Input/output, real X, Y.  On output, the values of X and
!!    Y have been interchanged.
!!
!  implicit none
!!
!  real x
!  real y
!  real z
!!
!  z = x
!  x = y
!  y = z
!
!  return
!end
!function r1mach(i)
!!
!!*******************************************************************************
!!
!!! R1MACH returns single precision machine constants.
!!
!!
!!  Assume that single precision numbers are stored with a mantissa of T digits
!!  in base B, with an exponent whose value must lie between EMIN and EMAX.  Then
!!  for values of I between 1 and 5, R1MACH will return the following values:
!!
!!    R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!!    R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
!!    R1MACH(3) = B**(-T), the smallest relative spacing.
!!    R1MACH(4) = B**(1-T), the largest relative spacing.
!!    R1MACH(5) = log10(B)
!!
!!  To alter this function for a particular environment, the desired set of data
!!  statements should be activated by removing the C from column 1.
!!
!!  On rare machines a STATIC statement may need to be added.  But probably more
!!  systems prohibit it that require it.
!!
!!  For IEEE-arithmetic machines (binary standard), the first set of constants
!!  below should be appropriate.
!!
!!  Where possible, octal or hexadecimal constants have been used to specify the
!!  constants exactly which has in some cases required the use of EQUIVALENCED
!!  integer arrays.  If your compiler uses half-word integers by default
!!  (sometimes called integer*2), you may need to change INTEGER to INTEGER*4 or
!!  otherwise instruct your compiler to use full-word integers in the next 5
!!  declarations.
!!
!  integer diver(2)
!  integer i
!  integer large(2)
!  integer log10(2)
!  real r1mach
!  integer right(2)
!  real rmach(5)
!  integer small(2)
!!
!  equivalence (rmach(1),small(1))
!  equivalence (rmach(2),large(1))
!  equivalence (rmach(3),right(1))
!  equivalence (rmach(4),diver(1))
!  equivalence (rmach(5),log10(1))
!!
!!  IEEE arithmetic machines, such as the ATT 3B series, Motorola 68000 based
!!  machines such as the SUN 3 and ATT PC 7300, and 8087 based micros such as
!!  the IBM PC and ATT 6300.
!!
!   data small(1) /     8388608 /
!   data large(1) /  2139095039 /
!   data right(1) /   864026624 /
!   data diver(1) /   872415232 /
!   data log10(1) /  1050288283 /
!!
!!  ALLIANT FX/8 UNIX Fortran compiler with the -r8 command line option.  This
!!  option causes all variables declared with 'real' to be of type 'REAL*8' or
!!  DOUBLE PRECISION.  This option does not override the 'real*4' declarations.
!!  These R1MACH numbers below and the coresponding I1MACH are simply the DOUBLE
!!  PRECISION or 'real*8' numbers.  If you use the -r8 your whole code (and the
!!  user libraries you link with, the system libraries are taken care of
!!  automagicly) must be compiled with this option.
!!
!!      data rmach(1) / 2.22507385850721D-308 /
!!      data rmach(2) / 1.79769313486231D+308 /
!!      data rmach(3) / 1.1101827117665D-16 /
!!      data rmach(4) / 2.2203654423533D-16 /
!!      data rmach(5) / 3.01029995663981E-1 /
!!
!!  AMDAHL machines.
!!
!!      data small(1) /    1048576 /
!!      data large(1) / 2147483647 /
!!      data right(1) /  990904320 /
!!      data diver(1) / 1007681536 /
!!      data log10(1) / 1091781651 /
!!
!!  BURROUGHS 1700 system.
!!
!!      data rmach(1) / Z400800000 /
!!      data rmach(2) / Z5FFFFFFFF /
!!      data rmach(3) / Z4E9800000 /
!!      data rmach(4) / Z4EA800000 /
!!      data rmach(5) / Z500E730E8 /
!!
!!  BURROUGHS 5700/6700/7700 systems.
!!
!!      data rmach(1) / O1771000000000000 /
!!      data rmach(2) / O0777777777777777 /
!!      data rmach(3) / O1311000000000000 /
!!      data rmach(4) / O1301000000000000 /
!!      data rmach(5) / O1157163034761675 /
!!
!!  CDC CYBER 170/180 series using NOS
!!
!!      data rmach(1) / O"00014000000000000000" /
!!      data rmach(2) / O"37767777777777777777" /
!!      data rmach(3) / O"16404000000000000000" /
!!      data rmach(4) / O"16414000000000000000" /
!!      data rmach(5) / O"17164642023241175720" /
!!
!!  CDC CYBER 170/180 series using NOS/VE
!!
!!      data rmach(1) / Z"3001800000000000" /
!!      data rmach(2) / Z"4FFEFFFFFFFFFFFE" /
!!      data rmach(3) / Z"3FD2800000000000" /
!!      data rmach(4) / Z"3FD3800000000000" /
!!      data rmach(5) / Z"3FFF9A209A84FBCF" /
!!
!!  CDC CYBER 200 series
!!
!!      data rmach(1) / X'9000400000000000' /
!!      data rmach(2) / X'6FFF7FFFFFFFFFFF' /
!!      data rmach(3) / X'FFA3400000000000' /
!!      data rmach(4) / X'FFA4400000000000' /
!!      data rmach(5) / X'FFD04D104D427DE8' /
!!
!!  CDC 6000/7000 series using FTN4.
!!
!!      data rmach(1) / 00564000000000000000B /
!!      data rmach(2) / 37767777777777777776B /
!!      data rmach(3) / 16414000000000000000B /
!!      data rmach(4) / 16424000000000000000B /
!!      data rmach(5) / 17164642023241175720B /
!!
!!  CDC 6000/7000 series using FTN5.
!!
!!      data rmach(1) / O"00564000000000000000" /
!!      data rmach(2) / O"37767777777777777776" /
!!      data rmach(3) / O"16414000000000000000" /
!!      data rmach(4) / O"16424000000000000000" /
!!      data rmach(5) / O"17164642023241175720" /
!!
!!  CONVEX C-1.
!!
!!      data rmach(1) / '00800000'X /
!!      data rmach(2) / '7FFFFFFF'X /
!!      data rmach(3) / '34800000'X /
!!      data rmach(4) / '35000000'X /
!!      data rmach(5) / '3F9A209B'X /
!!
!!  CONVEX C-120 (native mode) without -R8 option
!!
!!      data rmach(1) / 2.9387360E-39 /
!!      data rmach(2) / 1.7014117E+38 /
!!      data rmach(3) / 5.9604645E-08 /
!!      data rmach(4) / 1.1920929E-07 /
!!      data rmach(5) / 3.0102999E-01 /
!!
!!  CONVEX C-120 (native mode) with -R8 option
!!
!!      data rmach(1) / 5.562684646268007D-309 /
!!      data rmach(2) / 8.988465674311577D+307 /
!!      data rmach(3) / 1.110223024625157D-016 /
!!      data rmach(4) / 2.220446049250313D-016 /
!!      data rmach(5) / 3.010299956639812D-001 /
!!
!!  CONVEX C-120 (IEEE mode) without -R8 option
!!
!!      data rmach(1) / 1.1754945E-38 /
!!      data rmach(2) / 3.4028234E+38 /
!!      data rmach(3) / 5.9604645E-08 /
!!      data rmach(4) / 1.1920929E-07 /
!!      data rmach(5) / 3.0102999E-01 /
!!
!!  CONVEX C-120 (IEEE mode) with -R8 option
!!
!!      data rmach(1) / 2.225073858507202D-308 /
!!      data rmach(2) / 1.797693134862315D+308 /
!!      data rmach(3) / 1.110223024625157D-016 /
!!      data rmach(4) / 2.220446049250313D-016 /
!!      data rmach(5) / 3.010299956639812D-001 /
!!
!!  CRAY 1, 2, XMP and YMP.
!!
!!      data rmach(1) / 200034000000000000000B /
!!      data rmach(2) / 577767777777777777776B /
!!      data rmach(3) / 377224000000000000000B /
!!      data rmach(4) / 377234000000000000000B /
!!      data rmach(5) / 377774642023241175720B /
!!
!!  DATA GENERAL ECLIPSE S/200.
!!  Note - It may be appropriate to include the line: STATIC RMACH(5)
!!
!!      data small /20K,0/
!!      data large /77777K,177777K/
!!      data right /35420K,0/
!!      data diver /36020K,0/
!!      data log10 /40423K,42023K/
!!
!!  ELXSI 6400, assuming real*4 is the default real type.
!!
!!      data small(1) / '00800000'X /
!!      data large(1) / '7F7FFFFF'X /
!!      data right(1) / '33800000'X /
!!      data diver(1) / '34000000'X /
!!      data log10(1) / '3E9A209B'X /
!!
!!  HARRIS 220
!!
!!      data small(1),small(2) / '20000000, '00000201 /
!!      data large(1),large(2) / '37777777, '00000177 /
!!      data right(1),right(2) / '20000000, '00000352 /
!!      data diver(1),diver(2) / '20000000, '00000353 /
!!      data log10(1),log10(2) / '23210115, '00000377 /
!!
!!  HARRIS SLASH 6 and SLASH 7.
!!
!!      data small(1),small(2) / '20000000, '00000201 /
!!      data large(1),large(2) / '37777777, '00000177 /
!!      data right(1),right(2) / '20000000, '00000352 /
!!      data diver(1),diver(2) / '20000000, '00000353 /
!!      data log10(1),log10(2) / '23210115, '00000377 /
!!
!!  HONEYWELL DPS 8/70 and 600/6000 series.
!!
!!      data rmach(1) / O402400000000 /
!!      data rmach(2) / O376777777777 /
!!      data rmach(3) / O714400000000 /
!!      data rmach(4) / O716400000000 /
!!      data rmach(5) / O776464202324 /
!!
!!  HP 2100, 3 word double precision with FTN4
!!
!!      data small(1), small(2) / 40000B,       1 /
!!      data large(1), large(2) / 77777B, 177776B /
!!      data right(1), right(2) / 40000B,    325B /
!!      data diver(1), diver(2) / 40000B,    327B /
!!      data log10(1), log10(2) / 46420B,  46777B /
!!
!!  HP 2100, 4 word double precision with FTN4
!!
!!      data small(1), small(2) / 40000B,       1 /
!!      data large91), large(2) / 77777B, 177776B /
!!      data right(1), right(2) / 40000B,    325B /
!!      data diver(1), diver(2) / 40000B,    327B /
!!      data log10(1), log10(2) / 46420B,  46777B /
!!
!!  HP 9000
!!
!!      r1mach(1) = 1.17549435E-38
!!      r1mach(2) = 1.70141163E+38
!!      r1mach(3) = 5.960464478E-8
!!      r1mach(4) = 1.119209290E-7
!!      r1mach(5) = 3.01030010E-1
!!
!!      data small(1) / 00040000000B /
!!      data large(1) / 17677777777B /
!!      data right(1) / 06340000000B /
!!      data diver(1) / 06400000000B /
!!      data log10(1) / 07646420233B /
!!
!!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!!  and PERKIN ELMER (INTERDATA) 3230.
!!
!!      data rmach(1) / Z00100000 /
!!      data rmach(2) / Z7FFFFFFF /
!!      data rmach(3) / Z3B100000 /
!!      data rmach(4) / Z3C100000 /
!!      data rmach(5) / Z41134413 /
!!
!!  IBM PC - Microsoft FORTRAN
!!
!!      data small(1) / #00800000 /
!!      data large(1) / #7F7FFFFF /
!!      data right(1) / #33800000 /
!!      data diver(1) / #34000000 /
!!      data log10(1) / #3E9A209A /
!!
!!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!!
!!      data small(1)/ Z'00800000' /
!!      data large(1)/ Z'7F7FFFFF' /
!!      data right(1)/ Z'33800000' /
!!      data diver(1)/ Z'34000000' /
!!      data log10(1)/ Z'3E9A209A' /
!!
!!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!!  For the INTERDATA FORTRAN VII compiler replace the Z'S specifying HEX
!!  constants with Y'S.
!!
!!      data rmach(1) / Z'00100000' /
!!      data rmach(2) / Z'7EFFFFFF' /
!!      data rmach(3) / Z'3B100000' /
!!      data rmach(4) / Z'3C100000' /
!!      data rmach(5) / Z'41134413' /
!!
!!  PDP-10 (KA or KI processor).
!!
!!      data rmach(1) / "000400000000 /
!!      data rmach(2) / "377777777777 /
!!      data rmach(3) / "146400000000 /
!!      data rmach(4) / "147400000000 /
!!      data rmach(5) / "177464202324 /
!!
!!  PDP-11 FORTRANS supporting 32-bit integers (integer version).
!!
!!      data small(1) /    8388608 /
!!      data large(1) / 2147483647 /
!!      data right(1) /  880803840 /
!!      data diver(1) /  889192448 /
!!      data log10(1) / 1067065499 /
!!
!!  PDP-11 FORTRANS supporting 32-bit integers (octal version).
!!
!!      data rmach(1) / O00040000000 /
!!      data rmach(2) / O17777777777 /
!!      data rmach(3) / O06440000000 /
!!      data rmach(4) / O06500000000 /
!!      data rmach(5) / O07746420233 /
!!
!!  PDP-11 FORTRANS supporting 16-bit integers (integer version).
!!
!!      data small(1),small(2) /   128,     0 /
!!      data large(1),large(2) / 32767,    -1 /
!!      data right(1),right(2) / 13440,     0 /
!!      data diver(1),diver(2) / 13568,     0 /
!!      data log10(1),log10(2) / 16282,  8347 /
!!
!!  PDP-11 FORTRANS supporting 16-bit integers (octal version).
!!
!!      data small(1),small(2) / O000200, O000000 /
!!      data large(1),large(2) / O077777, O177777 /
!!      data right(1),right(2) / O032200, O000000 /
!!      data diver(1),diver(2) / O032400, O000000 /
!!      data log10(1),log10(2) / O037632, O020233 /
!!
!!  SEQUENT BALANCE 8000.
!!
!!      data small(1) / $00800000 /
!!      data large(1) / $7F7FFFFF /
!!      data right(1) / $33800000 /
!!      data diver(1) / $34000000 /
!!      data log10(1) / $3E9A209B /
!!
!!  SUN Microsystems UNIX F77 compiler.
!!
!!      data rmach(1) / 1.17549435E-38 /
!!      data rmach(2) / 3.40282347E+38 /
!!      data rmach(3) / 5.96016605E-08 /
!!      data rmach(4) / 1.19203321E-07 /
!!      data rmach(5) / 3.01030010E-01 /
!!
!!  SUN 3 (68881 or FPA)
!!
!!      data small(1) / X'00800000' /
!!      data large(1) / X'7F7FFFFF' /
!!      data right(1) / X'33800000' /
!!      data diver(1) / X'34000000' /
!!      data log10(1) / X'3E9A209B' /
!!
!!  UNIVAC 1100 series.
!!
!!      data rmach(1) / O000400000000 /
!!      data rmach(2) / O377777777777 /
!!      data rmach(3) / O146400000000 /
!!      data rmach(4) / O147400000000 /
!!      data rmach(5) / O177464202324 /
!!
!!  VAX/ULTRIX F77 compiler.
!!
!!      data small(1) /       128 /
!!      data large(1) /    -32769 /
!!      data right(1) /     13440 /
!!      data diver(1) /     13568 /
!!      data log10(1) / 547045274 /
!!
!!  VAX-11 with FORTRAN IV-PLUS compiler.
!!
!!      data rmach(1) / Z00000080 /
!!      data rmach(2) / ZFFFF7FFF /
!!      data rmach(3) / Z00003480 /
!!      data rmach(4) / Z00003500 /
!!      data rmach(5) / Z209B3F9A /
!!
!!  VAX/VMS version 2.2.
!!
!!      data rmach(1) /       '80'X /
!!      data rmach(2) / 'FFFF7FFF'X /
!!      data rmach(3) /     '3480'X /
!!      data rmach(4) /     '3500'X /
!!      data rmach(5) / '209B3F9A'X /
!!
!!  VAX/VMS 11/780
!!
!!      data small(1) / Z00000080 /
!!      data large(1) / ZFFFF7FFF /
!!      data right(1) / Z00003480 /
!!      data diver(1) / Z00003500 /
!!      data log10(1) / Z209B3F9A /
!!
!!  Z80 microprocessor.
!!
!!      data small(1), small(2) /     0,    256 /
!!      data large(1), large(2) /    -1,   -129 /
!!      data right(1), right(2) /     0,  26880 /
!!      data diver(1), diver(2) /     0,  27136 /
!!      data log10(1), log10(2) /  8347,  32538 /
!!
!  if ( i < 1 .or. i > 5 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'R1MACH - Fatal error!'
!    write ( *, '(a,i6)' ) 'I is out of bounds = ',i
!    r1mach = 0.0E+00
!    stop
!  else
!    r1mach = rmach(i)
!  end if
!
!  return
!end
!subroutine rminsp ( func, a, b, epsin, epsout, iop, result )
!!
!!***********************************************************************
!!
!!! RMINSP approximates the integral of a function using Romberg integration.
!!
!!
!!  Discussion:
!!
!!    Both the midpoint and trapezoidal rule are used,
!!    the intervals are repeatedly bisected, and Richardson
!!    extrapolation is carried out to achieve a high accuracy.
!!
!!    RMINSP can carry out a cosine-transform of the integral.  The
!!    only effect this has is to handle a function F(X) which has
!!    singularities near the endpoints.  This transform is based on
!!    the fact that
!!
!!      Integral from -1 to 1  ( F(    X))       DX
!!
!!    equals
!!
!!      Integral from  0 to PI ( F(COS(Z))*SIN(Z) )  DZ
!!
!!    If suitable accuracy is not achieved, the internal variable
!!    NUPPER might be increased.  Its current value of 9 corresponds
!!    to a maximum of 1024 subintervals and 1025 function evaluations.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!    T Havie,
!!    Algorithm 257,
!!    Communications of the Association for Computing Machinery,
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real, external FUNC, the name of the function to be
!!    integrated.  The user must declare the name an external
!!    parameter in the calling program, pass the name of the
!!    function in FUNC, and write a function of the form
!!
!!       FUNCTION FUNC(X)
!!
!!    which evaluates the function at the point X.
!!
!!    Input, real A, lower limit of integration.
!!
!!    Input, real B, upper limit of integration.
!!
!!    Input, real EPSIN, requested relative error tolerance.
!!
!!    Output, real EPSOUT, estimated achieved relative error.
!!
!!    Input, integer IOP, method switch:
!!    1, Use ordinary algorithm.
!!    2, Use cosine transformation to decrease effect of
!!       singularities near the endpoints.
!!
!!    Output, real RESULT, the approximate value of the integral.
!!
!  implicit none
!!
!  integer, parameter :: nupper = 9
!!
!  real a
!  real acof(11)
!  real alf
!  real alfnj
!  real alfno
!  real ar
!  real b
!  real bcof(nupper+1)
!  real bet
!  real betnj
!  real betno
!  real const1
!  real const2
!  real deltan
!  real endpts
!  real epsin
!  real epsout
!  real error
!  real, parameter :: fac1 = 0.411233516712057E+00
!  real, parameter :: fac2 = 0.822467033441132E+00
!  real factor
!  real, external :: func
!  real gamman
!  real hnstep
!  integer i
!  integer index
!  integer iop
!  integer iout
!  integer j
!  integer n
!  integer nhalf
!  real pi
!  real r1
!  real r2
!  real rn
!  real rnderr
!  real result
!  real rounde
!  real tend
!  real triarg
!  real umid
!  real xmin
!  real xplus
!!
!  result = 0.0E+00
! 
!  if ( a == b ) then
!    return
!  end if
!!
!!  Set coefficients in formula for accumulated roundoff error,
!!  rounde = rnderr*(r1+r2*n), where r1, r2 are two empirical constants
!!  and n is the current number of function values used.
!!
!  rnderr = epsilon ( 1.0E+00 )
! 
!  r1 = 1.0E+00
!  r2 = 2.0E+00
!  if ( iop==2 ) r1 = 50.0E+00
!  if ( iop==1 ) r2 = 0.01E+00 * r2
!  error = epsin
!!
!!  Initial calculations.
!!
!  alf = 0.5E+00 * (b-a)
!  bet = 0.5E+00 * (b+a)
!  acof(1) = func(a)+func(b)
!  bcof(1) = func(bet)
!!
!!  Modified Romberg algorithm, ordinary case.
!!
!  if ( iop /= 2 ) then
!
!    hnstep = 2.0E+00
!    bcof(1) = hnstep*bcof(1)
!    factor = 1.0E+00
!!
!!  Modified Romberg, cosine transformed case.
!!
!  else
!    hnstep = pi()
!    ar = fac1
!    endpts = acof(1)
!    acof(1) = fac2*acof(1)
!    bcof(1) = hnstep*bcof(1)-ar*endpts
!    factor = 4.0E+00
!    ar = ar/4.0E+00
!    triarg = pi() / 4.0E+00
!    alfno = -1.0E+00
!  end if
! 
!  hnstep = 0.5E+00 * hnstep
!  nhalf = 1
!  n = 2
!  rn = 2.0E+00
!  acof(1) = 0.5E+00 * (acof(1)+bcof(1))
!  acof(2) = acof(1)-(acof(1)-bcof(1))/(4.0E+00*factor-1.0E+00)
!!
!!  End of initial calculation.
!!
!!  Start actual calculations.
!!
!  do i = 1, nupper
! 
!    umid = 0.0E+00
!!
!!  Modified Romberg algorithm, ordinary case.
!!  compute first element in mid-point formula for ordinary case
!!
!    if ( iop == 1 ) then
! 
!      alfnj = 0.5E+00*hnstep
! 
!      do j = 1, nhalf
!        xplus = alf*alfnj+bet
!        xmin = -alf*alfnj+bet
!        umid = umid + func(xplus)+func(xmin)
!        alfnj = alfnj+hnstep
!      end do
! 
!      umid = hnstep*umid
!!
!!  Modified Romberg algorithm, cosine transformed case
!!  compute first element in mid-point formula for cosine transformed
!!  Romberg algorithm
!!
!    else if ( iop == 2 ) then
! 
!      const1 = -sin(triarg)
!      const2 = 0.5E+00 * alfno/const1
! 
!      alfno = const1
!      betno = const2
!      gamman = 1.0E+00 - 2.0E+00 * alfno**2
!      deltan = -2.0E+00 * alfno*betno
! 
!      do j = 1, nhalf
!        alfnj = gamman*const1+deltan*const2
!        betnj = gamman*const2-deltan*const1
!        xplus = alf*alfnj+bet
!        xmin = -alf*alfnj+bet
!        umid = umid+betnj*(func(xplus)+func(xmin))
!        const1 = alfnj
!        const2 = betnj
!      end do
! 
!      umid = hnstep*umid-ar*endpts
!      ar = ar / 4.0E+00
! 
!    end if
!!
!!  Modified Romberg algorithm, calculate (i+1)-th row in the U table
!!
!    const1 = 4.0E+00 * factor
!    index = i+1
! 
!    do j = 2, i+1
!      tend = umid + ( umid - bcof(j-1) ) / ( const1 - 1.0E+00 )
!      bcof(j-1) = umid
!      umid = tend
!      const1 = 4.0E+00 * const1
!    end do
! 
!    bcof(i+1) = tend
!    xplus = const1
!!
!!  Calculation of (i+1)-th row in the U table is finished
!!
!!  Test to see if the required accuracy has been obtained.
!!
!    epsout = 1.0E+00
!    iout = 1
! 
!    do j = 1, index
! 
!      const1 = 0.5E+00 * ( acof(j) + bcof(j) )
!      const2 = 0.5E+00 * abs ( ( acof(j) - bcof(j) ) / const1 )
! 
!      if ( const2 <= epsout ) then
!        epsout = const2
!        iout = j
!      end if
! 
!      acof(j) = const1
! 
!    end do
!!
!!  Testing on accuracy finished
!!
!    if (iout==index) iout=iout+1
!    acof(index+1) = acof(index)-(acof(index)-bcof(index))/(xplus-1.0)
!    rounde = rnderr*(r1+r2*rn)
!
!    epsout = max ( epsout, rounde )
!    error = max ( error, rounde )
! 
!    if ( epsout <= error ) go to 10
! 
!    nhalf = n
!    n = 2 * n
!    rn = 2.0E+00 * rn
!    hnstep = 0.5E+00 * hnstep
!
!    if ( iop > 1 ) then
!      triarg = 0.5 * triarg
!    end if
! 
!  end do
!!
!!  Accuracy not reached with maximum number of subdivisions
!!
!  n = nhalf
!!
!!  Calculation for modified Romberg algorithm finished
!!
!  10  continue
! 
!  n = 2*n
!  index = index-1
!  n = n+1
!  j = iout
!  if ((j-1)>=index) j = index
!  tend = alf * (2.0E+00*acof(j)-bcof(j))
!  umid = alf * bcof(j)
!  result = alf * acof(iout)
!
!  return
!end
!subroutine rvec_even ( alo, ahi, n, a )
!!
!!*******************************************************************************
!!
!!! RVEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!!
!!
!!  Modified:
!!
!!    31 October 2000
!!
!!  Author:
!!
!!    John Burkardt
!!
!!  Parameters:
!!
!!    Input, real ALO, AHI, the low and high values.
!!
!!    Input, integer N, the number of values.
!!
!!    Output, real A(N), N evenly spaced values.
!!    Normally, A(1) = ALO and A(N) = AHI.
!!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!!
!!
!  implicit none
!!
!  integer n
!!
!  real a(n)
!  real ahi
!  real alo
!  integer i
!!
!  if ( n == 1 ) then
!
!    a(1) = 0.5E+00 * ( alo + ahi )
!
!  else
!
!    do i = 1, n
!      a(i) = ( real ( n - i ) * alo + real ( i - 1 ) * ahi ) / real ( n - 1 )
!    end do
!
!  end if
!
!  return
!end
subroutine simp ( func, a, b, eps, result )
!
!***********************************************************************
!
!! SIMP approximates the integral of a function using an adaptive Simpson's rule.
!
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    J N Lyness,
!    Algorithm 379,
!    SQUANK - Simpson quadrature used adaptively, noise killed,
!    Communications of the Association for Computing Machinery,
!    Volume 13 (1970), pages 260-263.
!
!    W M McKeeman and L Tesler,
!    Algorithm 182,
!    Nonrecursive adaptive integration,
!    Communications of the Association for Computing Machinery,
!    Volume 6 (1963), page 315.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!    FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real A, the lower limit of integration.
!
!    Input, real B, the upper limit of integration.
!
!    Input, real EPS, the requested error tolerance.
!
!    Output, real RESULT, the approximation to the integral.
!
  implicit none
!
  integer(kind=4), parameter :: maxlev = 30
!
  real (kind=8) a
  real (kind=8) a1
  real (kind=8) absar
  real (kind=8) b
  real (kind=8) da
  real (kind=8) dx(maxlev)
  real (kind=8) ep
  real (kind=8) eps
  real (kind=8) epsp(maxlev)
  real (kind=8) est
  real (kind=8) est1
  real (kind=8) est2(maxlev)
  real (kind=8) est3(maxlev)
  real (kind=8) f1
  real (kind=8) f2(maxlev)
  real (kind=8) f3(maxlev)
  real (kind=8) f4(maxlev)
  real (kind=8) fa
  real (kind=8) fb
  real (kind=8) fbp(maxlev)
  real (kind=8) fm
  real (kind=8) fmp(maxlev)
  real (kind=8), external :: func
  integer (kind=4) i
  integer (kind=4) j
  integer (kind=4) l
  integer (kind=4) lvl
  integer (kind=4)nrtr(maxlev)
  real (kind=8) pval(maxlev,3)
  real (kind=8) result
  real (kind=8) sum1
  real (kind=8) sx
  real (kind=8) x2(maxlev)
  real (kind=8) x3(maxlev)
!
  result = 0.0E+00
  if ( a == b ) then
    return
  end if
 
  ep = eps
  a1 = a
  nrtr(1:maxlev) = 0
  pval(1:maxlev,1:3) = 0.0E+00
 
  lvl = 0
  absar = 0.0E+00
  est = 0.0E+00
  da = b - a1

  fa = func ( a1 )
  fm = 4.0E+00 * func ( (a1+b) * 0.5E+00 )
  fb = func ( b )
!
!  1 = RECUR
!
   30 continue
 
  lvl = lvl + 1
  dx(lvl) = da / 3.0E+00
  sx = dx(lvl)/6.0E+00
  f1 = 4.0E+00 * func(0.5*dx(lvl)+a1)
  x2(lvl) = a1+dx(lvl)
  f2(lvl) = func(x2(lvl))
  x3(lvl) = x2(lvl)+dx(lvl)
  f3(lvl) = func(x3(lvl))
  epsp(lvl) = ep
  f4(lvl) = 4.0E+00 * func(dx(lvl)*0.5E+00+x3(lvl))
  fmp(lvl) = fm
  est1 = sx*(fa+f1+f2(lvl))
  fbp(lvl) = fb
  est2(lvl) = sx * (f2(lvl)+f3(lvl)+fm)
  est3(lvl) = sx * (f3(lvl)+f4(lvl)+fb)
  sum1 = est1+est2(lvl)+est3(lvl)
  absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) &
    + abs ( est3(lvl) )
  if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) go to 40
  if ( lvl >= maxlev ) go to 50
!
!  2 = UP
!
40 continue
 
  if ( lvl > 1 ) then
    lvl = lvl-1
  end if

  l = nrtr(lvl)

  if ( l == 0 ) then
    go to 50
  end if

  pval(lvl,l) = sum1

  if ( l == 1 ) go to 60
  if ( l == 2 ) go to 70
  if ( l == 3 ) go to 80
 
50 continue
 
  nrtr(lvl) = 1
  est = est1
  fm = f1
  fb = f2(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
60 continue
 
  nrtr(lvl) = 2
  fa = f2(lvl)
  fm = fmp(lvl)
  fb = f3(lvl)
  est = est2(lvl)
  a1 = x2(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
70 continue
 
  nrtr(lvl) = 3
  fa = f3(lvl)
  fm = f4(lvl)
  fb = fbp(lvl)
  est = est3(lvl)
  a1 = x3(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
80 continue

  sum1 = pval(lvl,1)+pval(lvl,2)+pval(lvl,3)

  if ( lvl > 1 ) then
    go to 40
  end if
 
90 continue
 
  result = sum1
 
  return
end

subroutine simp2 ( func, a, b, eps, result )
!
!***********************************************************************
!
!! SIMP approximates the integral of a function using an adaptive Simpson's rule.
!
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    J N Lyness,
!    Algorithm 379,
!    SQUANK - Simpson quadrature used adaptively, noise killed,
!    Communications of the Association for Computing Machinery,
!    Volume 13 (1970), pages 260-263.
!
!    W M McKeeman and L Tesler,
!    Algorithm 182,
!    Nonrecursive adaptive integration,
!    Communications of the Association for Computing Machinery,
!    Volume 6 (1963), page 315.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!    FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real A, the lower limit of integration.
!
!    Input, real B, the upper limit of integration.
!
!    Input, real EPS, the requested error tolerance.
!
!    Output, real RESULT, the approximation to the integral.
!
implicit none
!
  integer(kind=4), parameter :: maxlev = 30
!
  real (kind=8) a
  real (kind=8) a1
  real (kind=8) absar
  real (kind=8) b
  real (kind=8) da
  real (kind=8) dx(maxlev)
  real (kind=8) ep
  real (kind=8) eps
  real (kind=8) epsp(maxlev)
  real (kind=8) est
  real (kind=8) est1
  real (kind=8) est2(maxlev)
  real (kind=8) est3(maxlev)
  real (kind=8) f1
  real (kind=8) f2(maxlev)
  real (kind=8) f3(maxlev)
  real (kind=8) f4(maxlev)
  real (kind=8) fa
  real (kind=8) fb
  real (kind=8) fbp(maxlev)
  real (kind=8) fm
  real (kind=8) fmp(maxlev)
  real (kind=8), external :: func
  integer (kind=4) i
  integer (kind=4) j
  integer (kind=4) l
  integer (kind=4) lvl
  integer (kind=4)nrtr(maxlev)
  real (kind=8) pval(maxlev,3)
  real (kind=8) result
  real (kind=8) sum1
  real (kind=8) sx
  real (kind=8) x2(maxlev)
  real (kind=8) x3(maxlev)
!
  result = 0.0E+00
  if ( a == b ) then
    return
  end if
 !print *, "limits inside intlib", a, b
  ep = eps
  a1 = a
  nrtr(1:maxlev) = 0
  pval(1:maxlev,1:3) = 0.0E+00
 
  lvl = 0
  absar = 0.0E+00
  est = 0.0E+00
  da = b - a1

  fa = func ( a1 )
  fm = 4.0E+00 * func ( (a1+b) * 0.5E+00 )
  fb = func ( b )
!
!  1 = RECUR
!
   30 continue
 
  lvl = lvl + 1
  dx(lvl) = da / 3.0E+00
  sx = dx(lvl)/6.0E+00
  f1 = 4.0E+00 * func(0.5*dx(lvl)+a1)
  x2(lvl) = a1+dx(lvl)
  f2(lvl) = func(x2(lvl))
  x3(lvl) = x2(lvl)+dx(lvl)
  f3(lvl) = func(x3(lvl))
  epsp(lvl) = ep
  f4(lvl) = 4.0E+00 * func(dx(lvl)*0.5E+00+x3(lvl))
  fmp(lvl) = fm
  est1 = sx*(fa+f1+f2(lvl))
  fbp(lvl) = fb
  est2(lvl) = sx * (f2(lvl)+f3(lvl)+fm)
  est3(lvl) = sx * (f3(lvl)+f4(lvl)+fb)
  sum1 = est1+est2(lvl)+est3(lvl)
  absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) &
    + abs ( est3(lvl) )
  if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) go to 40
  if ( lvl >= maxlev ) go to 50
!
!  2 = UP
!
40 continue
 
  if ( lvl > 1 ) then
    lvl = lvl-1
  end if

  l = nrtr(lvl)

  if ( l == 0 ) then
    go to 50
  end if

  pval(lvl,l) = sum1

  if ( l == 1 ) go to 60
  if ( l == 2 ) go to 70
  if ( l == 3 ) go to 80
 
50 continue
 
  nrtr(lvl) = 1
  est = est1
  fm = f1
  fb = f2(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
60 continue
 
  nrtr(lvl) = 2
  fa = f2(lvl)
  fm = fmp(lvl)
  fb = f3(lvl)
  est = est2(lvl)
  a1 = x2(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
70 continue
 
  nrtr(lvl) = 3
  fa = f3(lvl)
  fm = f4(lvl)
  fb = fbp(lvl)
  est = est3(lvl)
  a1 = x3(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
80 continue

  sum1 = pval(lvl,1)+pval(lvl,2)+pval(lvl,3)

  if ( lvl > 1 ) then
    go to 40
  end if
 
90 continue
 
  result = sum1
 
  return
end

!subroutine simpne ( x, y, num, result )
!!
!!***********************************************************************
!!
!!! SIMPNE approximates the integral of unevenly spaced data.
!!
!!
!!  Discussion:
!!
!!    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
!!    to the data and integrates that exactly.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real X(NUM), contains the X values of the data, in order.
!!
!!    Input, real Y(NUM), contains the Y values of the data.
!!
!!    Input, integer NUM, number of data points.  NUM must be at least 3.
!!
!!    Output, real RESULT.
!!    RESULT is the approximate value of the integral.
!!
!  implicit none
!!
!  integer num
!!
!  real del(3)
!  real e
!  real f
!  real feints
!  real g(3)
!  integer i
!  integer n
!  real pi(3)
!  real result
!  real sum1
!  real x(num)
!  real x1
!  real x2
!  real x3
!  real y(num)
!!
!  result = 0.0E+00
! 
!  if ( num <= 2 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'SIMPNE - Fatal error!'
!    write ( *, '(a)' ) '  NUM <= 2.'
!    stop
!  end if
! 
!  n = 1
! 
!  do
! 
!    x1 = x(n)
!    x2 = x(n+1)
!    x3 = x(n+2)
!    e = x3*x3-x1*x1
!    f = x3*x3*x3-x1*x1*x1
!    feints = x3-x1
!    del(1) = x3-x2
!    del(2) = x1-x3
!    del(3) = x2-x1
!    g(1) = x2+x3
!    g(2) = x1+x3
!    g(3) = x1+x2
!    pi(1) = x2*x3
!    pi(2) = x1*x3
!    pi(3) = x1*x2
! 
!    sum1 = 0.0E+00
!    do i = 1, 3
!      sum1 = sum1 + y(n-1+i)*del(i)*(f/3.0E+00-g(i)*0.5E+00*e+pi(i)*feints)
!    end do
!    result = result - sum1 / ( del(1) * del(2) * del(3) )
! 
!    n = n+2
!
!    if ( n + 1 >= num ) then
!      exit
!    end if
!
!  end do
! 
!  if ( mod(num,2) /= 0 ) then
!    return
!  end if
!
!  n = num-2
!  x3 = x(num)
!  x2 = x(num-1)
!  x1 = x(num-2)
!  e = x3*x3-x2*x2
!  f = x3*x3*x3-x2*x2*x2
!  feints = x3-x2
!  del(1) = x3-x2
!  del(2) = x1-x3
!  del(3) = x2-x1
!  g(1) = x2+x3
!  g(2) = x1+x3
!  g(3) = x1+x2
!  pi(1) = x2*x3
!  pi(2) = x1*x3
!  pi(3) = x1*x2
! 
!  sum1 = 0.0E+00
!  do i = 1, 3
!    sum1 = sum1 + y(n-1+i) * del(i) * &
!      ( f / 3.0E+00 - g(i) * 0.5E+00 * e + pi(i) * feints )
!  end do
! 
!  result = result - sum1 / ( del(1) * del(2) * del(3) )
! 
!  return
!end
!subroutine simpsn ( h, y, num, result )
!!
!!***********************************************************************
!!
!!! SIMPSN approximates the integral of evenly spaced data.
!!
!!
!!  Discussion:
!!
!!    Simpson's rule is used.
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real H, specifies the increment between the
!!    X values.  Note that the actual X values are not needed,
!!    just the constant spacing!
!!
!!    Input, real Y(NUM), the data.
!!
!!    Input, integer NUM, the number of data points.  NUM must be at least 3.
!!
!!    Output, real RESULT, the value of the integral
!!    from the first to the last point.
!!
!  implicit none
!!
!  integer num
!!
!  real del(3)
!  real f
!  real g(3)
!  real h
!  integer i
!  integer n
!  real pii(3)
!  real result
!  real sum1
!  real y(num)
!!
!  result = 0.0E+00
! 
!  if ( num <= 2 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'SIMPSN - Fatal error!'
!    write ( *, '(a,i6)' ) '  NUM < 2, NUM = ', num
!    stop
!  end if
! 
!  if ( mod ( num, 2 ) == 0 ) then
!    n = num-1
!  else
!    n = num
!  end if
! 
!  result = y(1) + y(n) + 4.0E+00 * y(n-1)
!  do i = 2, n-2, 2
!    result = result + 4.0E+00 * y(i) + 2.0E+00 * y(i+1)
!  end do
!  result = h * result / 3.0E+00
! 
!  if ( mod(num,2) == 1 ) then
!    return
!  end if
! 
!  f = h*h*h
!  del(1) = h
!  del(2) = -2.0E+00 * h
!  del(3) = h
!  g(1) = h
!  g(2) = 0.0E+00
!  g(3) = -h
!  pii(1) = 0.0E+00
!  pii(2) = -h*h
!  pii(3) = 0.0E+00
!  n = n-1
! 
!  sum1 = 0.0E+00
!  do i = 1, 3
!    sum1 = sum1 + y(n-1+i) * del(i) * &
!      ( f / 3.0E+00 - g(i) * 0.5E+00 * h * h + pii(i) * h )
!  end do
! 
!  result = result + 0.5E+00 * sum1 / h**3
! 
!  return
!end
!function solve ( shift, n, a, b )
!!
!!***********************************************************************
!!
!!! SOLVE solves a special linear system.
!!
!!
!!  Discussion:
!!
!!    SOLVE solves for the N-th component of the solution DELTA to the equation
!!
!!      (Jn - shift*Identity) * DELTA  = En,
!!
!!    En is the vector of all zeroes except for 1 in the N-th position.
!!
!!    The matrix Jn is symmetric tridiagonal, with diagonal
!!    elements A(I), off-diagonal elements B(I).  This equation
!!    must be solved to obtain the appropriate changes in the lower
!!    2 by 2 submatrix of coefficients for orthogonal polynomials.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Parameters:
!!
!!    Input, real SHIFT, the value of the factor that multiplies
!!    the identity matrix, in the definition of the system matrix.
!!
!!    Input, integer N, the index of the desired component.
!!
!  implicit none
!!
!  integer n
!!
!  real a(n)
!  real alpha
!  real b(n)
!  integer i
!  real shift
!  real solve
!!
!  alpha = a(1) - shift
!  do i = 2, n-1
!    alpha = a(i) - shift - b(i-1)**2 / alpha
!  end do
! 
!  solve = 1.0E+00 / alpha
! 
!  return
!end
!!subroutine timestamp ( )
!!!
!!!*******************************************************************************
!!!
!!!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!!!
!!!
!!!  Example:
!!!
!!!    May 31 2001   9:45:54.872 AM
!!!
!!!  Modified:
!!!
!!!    31 May 2001
!!!
!!!  Author:
!!!
!!!    John Burkardt
!!!
!!!  Parameters:
!!!
!!!    None
!!!
!!  implicit none
!!!
!!  character ( len = 8 ) ampm
!!  integer d
!!  character ( len = 8 ) date
!!  integer h
!!  integer m
!!  integer mm
!!  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
!!    'January  ', 'February ', 'March    ', 'April    ', &
!!    'May      ', 'June     ', 'July     ', 'August   ', &
!!    'September', 'October  ', 'November ', 'December ' /)
!!  integer n
!!  integer s
!!  character ( len = 10 )  time
!!  integer values(8)
!!  integer y
!!  character ( len = 5 ) zone
!!!
!!  call date_and_time ( date, time, zone, values )
!!
!!  y = values(1)
!!  m = values(2)
!!  d = values(3)
!!  h = values(5)
!!  n = values(6)
!!  s = values(7)
!!  mm = values(8)
!!
!!  if ( h < 12 ) then
!!    ampm = 'AM'
!!  else if ( h == 12 ) then
!!    if ( n == 0 .and. s == 0 ) then
!!      ampm = 'Noon'
!!    else
!!      ampm = 'PM'
!!    end if
!!  else
!!    h = h - 12
!!    if ( h < 12 ) then
!!      ampm = 'PM'
!!    else if ( h == 12 ) then
!!      if ( n == 0 .and. s == 0 ) then
!!        ampm = 'Midnight'
!!      else
!!        ampm = 'AM'
!!      end if
!!    end if
!!  end if
!!
!!  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
!!    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )
!!
!!  return
!!end
!subroutine wedint ( ftab, h, ntab, result )
!!
!!***********************************************************************
!!
!!! WEDINT uses Weddle's rule to integrate data at equally spaced points.
!!
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Author:
!!
!!    John Burkardt
!!
!!  Parameters:
!!
!!    Input, real FTAB(NTAB), contains the tabulated data values.
!!
!!    Input, real H, is the spacing between the points at which the data
!!    was evaluated.
!!
!!    Input, integer NTAB, is the number of data points.  (NTAB-1) must be
!!    divisible by 6.
!!
!!    Output, real RESULT, is the approximation to the integral.
!!
!  implicit none
!!
!  integer ntab
!!
!  real ftab(ntab)
!  real h
!  integer i
!  real result
!!
!  result = 0.0E+00
! 
!  if ( ntab <= 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'WEDINT - Fatal error!'
!    write ( *, '(a)' ) '  NTAB < 2'
!    write ( *, '(a,i6)' ) '  NTAB = ', ntab
!    stop
!  end if
! 
!  if ( mod ( ntab, 6 ) /= 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'WEDINT - Fatal error!'
!    write ( *, '(a)' ) '  NTAB must equal 6*N+1 for some N!'
!    stop
!  end if
! 
!  do i = 1, ntab-6, 6
!    result = result & 
!      +           ftab(i) &
!      + 5.0E+00 * ftab(i+1) &
!      +           ftab(i+2) &
!      + 6.0E+00 * ftab(i+3) &
!      +           ftab(i+4) &
!      + 5.0E+00 * ftab(i+5) &
!      +           ftab(i+6)
!  end do
! 
!  result = 3.0E+00 * h * result / 10.0E+00
! 
!  return
!end


subroutine specialpositions(myx, myy)
	!determine special positions on the interface
	use globalinfotest
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
		!call interpbridge(5, myy(tempmaxi:tempmaxi+4), myx(tempmaxi:tempmaxi+4),  yend, xflagb)
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
			!call interpbridge(min(flagb+2, XN+1), myy(1:min(flagb+2, XN+1)), myx(1:min(flagb+2, XN+1)), -R, xflagl)
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





subroutine trapz1(f,a,b,h,r)
!===========================================================
! int_trap.f: integration by trapezoid rule of f(x) on [a,b]
!-------------------------------------------------
! f     - Function to integrate (supplied above)
! a	- Lower limit of integration
! b	- Upper limit of integration
! R	- Result of integration (out)
! n	- number of intervals
!=================================================
Real*8 a, b, f, r ,dx, x,h
!double precision a,b,f,r,dx,x
Integer*4 n, i
! interger*4 i


!correct to even

n= (b-a)/h + 1

n=2*ceiling(((b-a)/h + 1)/2)


r = 0.d0

dx = (b-a)/(n-1)


Do i=2,n-1
x = a+(i-1)*dx
r = r + f(x)
End do


r = (r + (f(a)+f(b))/2.d0)*dx

Return
end subroutine trapz1








 
 
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

use globalinfotest
implicit none

Real*8  ff, ss
real (kind=8) xy(XN+1)
Integer*4  ii,aa,bb

 
ss= 0.d0

Do ii=aa,bb-1 
 
ss= ss + 0.5*(ff(ii+1)+ff(ii))*(xy(ii+1)-xy(ii))

  ! print *, "ii", ii
 !  print *, ff(ii+1), ff(ii)
!print *, "xy(ii+1)", xy(ii+1),xy(ii)
!print *, "ss", ss
End do


Return
end subroutine trapznu


