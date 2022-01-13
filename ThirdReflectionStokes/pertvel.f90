subroutine wTN () !(wu, wv)
use globalinfo
implicit none



real (kind = 8)  wbackflow(XN+1), wsidesphere(XN+1), wbelowsphere(XN+1), wabovesphere(XN+1),&
tempbackflow, tempsidesphere, tempbelowsphere, tempabovesphere !, wu(XN+1), wv(XN+1),
real (kind=8), external :: w2IntegrandBackflow, w2IntegrandBelowSphere, w2IntegrandPartialSphere,&
w2IntegrandZetaSphere, w2IntegrandZetaVert, w2IntegrandR, w2IntegrandZ
integer (kind=4) wi
real (kind=8) AEbelow,AEside,AEabove!,AreaEntrain, AreaReflux, AreaSpherePortion
real (kind=8), external :: RefluxAreaFunction,EntrainPartialSphere,EntrainZetaSphere,EntrainZetaVert,EntrainBelowSphere, AreaElementZeta,AreaElementRho

integer  ierr


wbackflow		= real(0.0, kind=8)
wsidesphere 	= real(0.0, kind=8)
wbelowsphere 	= real(0.0, kind=8)
wabovesphere 	= real(0.0, kind=8)
AreaReflux      = real(0.0, kind=8)
AreaEntrain     = real(0.0, kind=8)

call specialpositions(sx, sy)





!========== FIND ENTRAINMENT AND REFLUX VOLUMES =========!

!reflux

if(flagb <XN+1.0 ) then
!print *, "w0"
call simp(RefluxAreaFunction, xflagb, sx(XN+1), integthres, AreaReflux)

endif
if (flagu /= 0.0) then
AreaSpherePortion= (4/3*pi*R**3);


call simp(EntrainPartialSphere,real(0.0,kind=8),xflagl,integthres,AEbelow)


call simp(EntrainZetaSphere, max(-R, sy(1)), R, integthres, AEside)


call simp(EntrainZetaVert, R, yend, integthres, AEabove)

AreaEntrain=AEbelow+AEside+AEabove;

elseif (flagl /= 0.0) then


!find area of the sphere

AreaSpherePortion= (pi*(yend+R)**2/3)*(3*R - (yend+R));


call simp(EntrainPartialSphere, real(0.0, kind=8), xflagl, integthres,AEbelow)

call simp(EntrainZetaSphere, max(-R, sy(1)), yend, integthres, AEside)


AreaEntrain=AEbelow+AEside;


else
!print *, "w6"
AreaSpherePortion= 0.0;

call simp(EntrainBelowSphere, max(sx(1), real(0.0, kind=8)), xflagb, integthres, AEbelow)

AreaEntrain=AEbelow;

endif



!========== FIND perturbation velocity w  =========!

do RorZ = 0, 1
if (minval(sy) < maxval(sy)) then
do wi = 1, XN+1
px = sx(wi)
py = sy(wi)



!========= Initialize Interpolated  ===== !



IF (ALLOCATED(cwinterpx)) DEALLOCATE(cwinterpx, STAT=ierr)
IF (ALLOCATED(cwinterpy)) DEALLOCATE(cwinterpy, STAT=ierr)

!initialize interface interpolated
ALLOCATE(cwinterpx(1000000), STAT=ierr)
IF (ierr /= 0) PRINT*, "cintertpx : Allocation failed"

ALLOCATE(cwinterpy(1000000), STAT=ierr)
IF (ierr /= 0) PRINT*, "cintertpy : Allocation failed"

cwinternalcount=0.0;
cwinterpx =0;
cwinterpy = 0;



tempbackflow	= real(0.0, kind=8)
tempsidesphere	= real(0.0, kind=8)
tempbelowsphere	= real(0.0, kind=8)
tempabovesphere	= real(0.0, kind=8)

if(flagb <XN+1.0 ) then
!print *, "w0"
call simp(w2IntegrandBackflow, xflagb, sx(XN+1), integthres, tempbackflow)

endif
if (flagu /= 0.0) then
!print *, "w1"
call simp(w2IntegrandPartialSphere, real(0.0, kind=8), xflagl, integthres, tempbelowsphere)
!print *, "w2"
call simp(w2IntegrandZetaSphere, max(-R, sy(1)), R, integthres, tempsidesphere)
!print *, "w3"
call simp(w2IntegrandZetaVert, R, yend, integthres, tempabovesphere)
!print *, "w3.5"
elseif (flagl /= 0.0) then
!print *, "w4"
call simp(w2IntegrandPartialSphere, real(0.0, kind=8), xflagl, integthres, tempbelowsphere)
!print *, "w5"
call simp(w2IntegrandZetaSphere, max(-R, sy(1)), yend, integthres, tempsidesphere)
!print *, "w5.5"
else
!print *, "w6"

!print *, "below sphere limits", real(0.0, kind=8), xflagb
call simp(w2IntegrandBelowSphere, max(sx(1), real(0.0, kind=8)), xflagb, integthres, tempbelowsphere)

!print *, "belowsphere", tempbelowsphere

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
!			print *, "wbackflow", wbackflow
!			print *, "wbelowsphere",wbelowsphere(2)
!			print *, "wsidesphere",wsidesphere(2)
!			print *, "wabovesphere",wabovesphere(2)
endif
!print *, "wv", wv


end do

!print *, "flow components", wbelowsphere(3), wsidesphere(3)
if(abs(wbelowsphere(1)) > 1000.0) then
stop
endif
!print *, "w7"
end subroutine wTN



!========== FOR VOLUME TRACK =========!

function AreaElementRho (rho)
use globalinfo
implicit none

real (kind =8) rho, AreaElementRho

AreaElementRho = 2*pi*rho;
return
end


function AreaElementZeta (zeta)
use globalinfo
implicit none

real (kind =8) zeta, AreaElementZeta

AreaElementZeta = 2*pi*myrho;
return
end


function RefluxAreaFunction(rho)
use globalinfo
implicit none

real (kind=8) rho, zcoord, RefluxAreaFunction
real (kind=8), external :: AreaElementZeta

call interpbridge( XN+2-max(flagb-2, 1), sx(max(flagb-2,1):XN+1), sy(max(flagb-2, 1):XN+1), rho, zcoord)

myrho   = rho
call simp2(AreaElementZeta, yend, zcoord, integthres, RefluxAreaFunction)
end


function EntrainZetaSphere(zeta)
use globalinfo
implicit none

real (kind=8) zeta, xupper, xlower, EntrainZetaSphere
real (kind=8), external :: AreaElementRho

myzeta  = zeta
call interpbridge(min(flagb+2, XN+1), sy(1:min(flagb+2, XN+1)), sx(1:min(flagb+2, XN+1)), zeta, xupper)

xlower  = sqrt(R**2 - zeta**2)

call simp2(AreaElementRho, xlower, xupper, integthres, EntrainZetaSphere)
end

function EntrainPartialSphere(rho)
use globalinfo
implicit none

real (kind=8) rho, zcoord, EntrainPartialSphere
real (kind=8), external :: AreaElementZeta

myrho = rho

call interpbridge(min(flagl+3, XN+1), sx(1:min(flagl+3, XN+1)), sy(1:min(flagl+3, XN+1)), rho, zcoord)


call simp2(AreaElementZeta, zcoord, -R, integthres, EntrainPartialSphere)
end

function EntrainZetaVert(zeta)
use globalinfo
implicit none

real (kind=8) zeta, xupper, EntrainZetaVert
real (kind=8), external :: AreaElementRho
integer (kind=4) temp(1), tempmini

myzeta  = zeta

call interpbridge(min(flagb+2, XN+1), sy(1:min(flagb+2, XN+1)), sx(1:min(flagb+2, XN+1)), zeta, xupper)

call simp2(AreaElementRho, real(0.0, kind=8), xupper, integthres,EntrainZetaVert)
end


function EntrainBelowSphere(rho)
use globalinfo
implicit none

real (kind=8) rho, zcoord, EntrainBelowSphere
real (kind=8), external :: AreaElementZeta



call interpbridge(max(flagl+2, XN+1), sx(1:max(flagl+2, XN+1)), sy(1:max(flagl+2, XN+1)), rho, zcoord)

myrho   = rho
call simp2(AreaElementZeta, zcoord, yend, integthres, EntrainBelowSphere)
end


!========== FOR VOLUME TRACK =========!

function w2IntegrandBackflow(rho)
use globalinfo
implicit none

real (kind=8) rho, zcoord, w2IntegrandBackflow
real (kind=8), external :: w2IntegrandZeta

!print *, "w1"
call interpbridge( XN+2-max(flagb-2, 1), sx(max(flagb-2,1):XN+1), sy(max(flagb-2, 1):XN+1), rho, zcoord)
!print *, "w1"
myrho   = rho
call simp2(w2IntegrandZeta, yend, zcoord, integthres, w2IntegrandBackflow)

cwinternalcount = cwinternalcount +1;
cwinterpx(cwinternalcount) = rho;
cwinterpy(cwinternalcount) = zcoord;

end

function w2IntegrandZetaSphere(zeta)
use globalinfo
implicit none

real (kind=8) zeta, xupper, xlower, w2IntegrandZetaSphere
real (kind=8), external :: w2IntegrandRho

myzeta  = zeta
!print *, "w2"
call interpbridge(min(flagb+2, XN+1), sy(1:min(flagb+2, XN+1)), sx(1:min(flagb+2, XN+1)), zeta, xupper)
!print *, "w2"
xlower  = sqrt(R**2 - zeta**2)
call simp2(w2IntegrandRho, xlower, xupper, integthres, w2IntegrandZetaSphere)

cwinternalcount = cwinternalcount +1;
cwinterpx(cwinternalcount) = xupper;
cwinterpy(cwinternalcount) = zeta;
end

function w2IntegrandPartialSphere(rho)
use globalinfo
implicit none

real (kind=8) rho, zcoord, w2IntegrandPartialSphere
real (kind=8), external :: w2IntegrandZeta

myrho = rho
!print *, "w3"
call interpbridge(min(flagl+3, XN+1), sx(1:min(flagl+3, XN+1)), sy(1:min(flagl+3, XN+1)), rho, zcoord)
!print *, "w3"
call simp2(w2IntegrandZeta, zcoord, -R, integthres, w2IntegrandPartialSphere)

cwinternalcount = cwinternalcount +1;
cwinterpx(cwinternalcount) = rho;
cwinterpy(cwinternalcount) = zcoord;
end

function w2IntegrandZetaVert(zeta)
use globalinfo
implicit none

real (kind=8) zeta, xupper, w2IntegrandZetaVert
real (kind=8), external :: w2IntegrandRho
integer (kind=4) temp(1), tempmini

myzeta  = zeta
!temp = minloc( sy(1:min(flagb, XN+1)))
!tempmini = temp(1)
!print *, "w4"
!call interpbridge(min(flagb, XN+1)-tempmini+1, sy(tempmini:min(flagb, XN+1)), sx(tempmini:min(flagb, XN+1)),zeta, xupper)
call interpbridge(min(flagb+2, XN+1), sy(1:min(flagb+2, XN+1)), sx(1:min(flagb+2, XN+1)), zeta, xupper)
!print *, "w4"
!print *, yend
call simp2(w2IntegrandRho, real(0.0, kind=8), xupper, integthres,w2IntegrandZetaVert)


cwinternalcount = cwinternalcount +1;
cwinterpx(cwinternalcount) = xupper;
cwinterpy(cwinternalcount) = zeta;
end


function w2IntegrandBelowSphere(rho)
use globalinfo
implicit none

real (kind=8) rho, zcoord, w2IntegrandBelowSphere
real (kind=8), external :: w2IntegrandZeta


!print *, "w5p2"
call interpbridge(max(flagl+2, XN+1), sx(1:max(flagl+2, XN+1)), sy(1:max(flagl+2, XN+1)), rho, zcoord)
!print *, "w5p2"
myrho   = rho
!print *, zcoord, yend, myrho
call simp2(w2IntegrandZeta, zcoord, yend, integthres, w2IntegrandBelowSphere)

cwinternalcount = cwinternalcount +1;
cwinterpx(cwinternalcount) = rho;
cwinterpy(cwinternalcount) = zcoord;
end




!=====================================================
!Integrands
!=====================================================

function w2IntegrandZeta(zeta)
use globalinfo
implicit none

real (kind=8) zeta, w2IntegrandZeta
real (kind=8), external :: w2IntegrandR, w2IntegrandZ, w2IntegrandRFF, w2IntegrandZFF



Real*8 ellipticE, ellipticK,ellipticE1, ellipticK1
DOUBLE PRECISION k,k1, kbar, tempK, tempE, DRF, DRD, ex, ey, ez
integer ier



myzeta = zeta


if ( ((py-myzeta)**2+(px-myrho)**2)>singthres) then

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



if (RorZ >0.0) then

if (px**2+py**2 >= (FFR)**2) then
w2IntegrandZeta= w2IntegrandRFF()
else
w2IntegrandZeta =w2IntegrandR(ellipticK,ellipticE,ellipticK1,ellipticE1)
endif

else
if (px**2+py**2 >= (FFR)**2) then
w2IntegrandZeta = w2IntegrandZFF()
!print *,'FarField', 'px',px,'py',py, 'zeta',zeta,' rho',myrho
!print *, "integrand value", w2IntegrandZeta
else

w2IntegrandZeta =w2IntegrandZ(ellipticK,ellipticE,ellipticK1,ellipticE1)
!print*,'3D Oseen',  'px',px,'py',py, 'zeta',zeta,' rho',myrho
!print *, "integrand value", w2IntegrandZeta
endif

endif





!print *, "integrand value", w2IntegrandZeta
!print *, "eval pts", zeta, px, py
!print *, "eval pt and result", myzeta,  w2IntegrandZeta, myt
!endif
end



function w2IntegrandRho(rho)
use globalinfo
implicit none

real (kind=8) rho, w2IntegrandRho
real (kind=8), external :: w2IntegrandR, w2IntegrandZ, w2IntegrandRFF, w2IntegrandZFF


Real*8 ellipticE, ellipticK,ellipticE1, ellipticK1
DOUBLE PRECISION k, k1, kbar, tempK, tempE, DRF, DRD, ex, ey, ez
integer ier




myrho = rho
if ( ((py-myzeta)**2+(px-myrho)**2)>singthres) then

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


if (RorZ >0.0) then

if (px**2+py**2 >= (FFR)**2) then
w2IntegrandRho = w2IntegrandRFF()
else
w2IntegrandRho = w2IntegrandR(ellipticK,ellipticE,ellipticK1,ellipticE1)
endif

else
if (px**2+py**2 >= (FFR)**2) then
w2IntegrandRho = w2IntegrandZFF()
else
w2IntegrandRho = w2IntegrandZ(ellipticK,ellipticE,ellipticK1,ellipticE1)
endif

endif
end





!====================== Far Field kernel for px^2+py^2 > 4R^2 =====================!
function w2IntegrandRFF()
use globalinfo
implicit none

real (kind=8) w2IntegrandRFF, ellipticE, ellipticK
DOUBLE PRECISION k, kbar, tempK, tempE, DRF, DRD, ex, ey, ez
integer ier

w2IntegrandRFF = 0.0
if (px > 0.0 .and. ((py-myzeta)**2+(px-myrho)**2)>singthres) then
k       = 4.0*px*myrho/((py-myzeta)**2+(px-myrho)**2)
kbar    = k/(k+1.0)

ex = 0.0
ey = 1.0-kbar
ez = 1.0

tempK 	= DRF(ex, ey, ez, ier)
tempE 	= tempK-1.0/3.0*kbar*DRD(ex, ey, ez, ier)
!tempK=0.0
!tempE=0.0

ellipticE	= sqrt(1.0+k)*tempE
ellipticK	= (1.0/sqrt(1.0+k))*tempK

w2IntegrandRFF = 2.0*(py-myzeta)/(pi*px*sqrt((px-myrho)**2+(py-myzeta)**2)&
*((px+myrho)**2+(py-myzeta)**2))*((px**2-myrho**2-(py-myzeta)**2)*ellipticE&
+((myrho+px)**2+(py-myzeta)**2)*ellipticK)&
-3.0*R*px*py*(2.0*myzeta**2+myrho**2)/(2.0*sqrt(px**2+py**2)**3*sqrt(myzeta**2+myrho**2)**3)&
+(R**3*px*(px**2*(py + 5.0*myzeta)* (2.0*myzeta**2 - myrho**2)&
+ py*(py**2*(2.0*myzeta**2 - myrho**2) + 10.0*py*myzeta*(-2.0*myzeta**2 + myrho**2)&
+3.0*(2.0*myzeta**4 + 3.0*myzeta**2*myrho**2 + myrho**4))))/(2.0*(px**2 + py**2)**(2.5)*(myzeta**2 + myrho**2)**(2.5))&
- (3.0*R**5*px*(8.0*px**4*(2.0*myzeta**3 - 3.0*myzeta*myrho**2) +px**2* (-8.0*py**2*(2.0*myzeta**3&
- 3.0*myzeta*myrho**2) + py*(-136.0*myzeta**4 + 296.0*myzeta**2*myrho**2 - 23.0*myrho**4)&
+ 8.0*myzeta*(2.0*myzeta**4 + myzeta**2*myrho**2 - myrho**4)) -4.0*py**2*(4.0*py**2*(2.0*myzeta**3&
- 3.0*myzeta*myrho**2) + 8.0*myzeta*(2.0*myzeta**4 + myzeta**2*myrho**2 - myrho**4)&
- 3.0*py*(12.0*myzeta**4 - 22.0*myzeta**2* myrho**2 + myrho**4))))/(16.0*(px**2 + py**2)**(3.5)*&
(myzeta**2 + myrho**2)**(3.5))

w2IntegrandRFF = w2IntegrandRFF*myrho
endif
end

function w2IntegrandZFF()
use globalinfo
implicit none

real (kind=8) w2IntegrandZFF, ellipticE, ellipticK
DOUBLE PRECISION k, kbar, tempK, tempE, DRF, DRD, ex, ey, ez
integer ier

w2IntegrandZFF = 0.0
!if (px > 0.0 .and. ((py-myzeta)**2+(px-myrho)**2)>1.0E-006) then
if (((py-myzeta)**2+(px-myrho)**2)>singthres) then
k       = 4.0*px*myrho/((py-myzeta)**2+(px-myrho)**2)
kbar    = k/(k+1.0)

ex = 0.0
ey = 1.0-kbar
ez = 1.0

tempK 	= DRF(ex, ey, ez, ier)
tempE 	= DRF(ex, ey, ez, ier)-1.0/3.0*kbar*DRD(ex, ey, ez, ier)
!tempK=0.0
!tempE=0.0

ellipticE	= sqrt(1.0+k)*tempE
ellipticK	= (1.0/sqrt(1.0+k))*tempK

w2IntegrandZFF = 4.0*((py-myzeta)**2*ellipticE+((px+myrho)**2+(py-myzeta)**2)*ellipticK)/&
(pi*sqrt((px-myrho)**2+(py-myzeta)**2)*((px+myrho)**2+(py-myzeta)**2))&
-3.0*R*(px**2+2.0*py**2)*(2.0*myzeta**2+myrho**2)/(2.0*sqrt(px**2+py**2)**3*sqrt(myzeta**2+myrho**2)**3)&
-R**3/(2.0*sqrt(px**2+py**2)**5*sqrt(myrho**2+myzeta**2)**5)*&
(px**4*(-2.0*myzeta**2 + myrho**2) - 2.0*py**2*(2.0*myzeta**4 + 3.0*myzeta**2*myrho**2&
+ myrho**4 + py**2*(2.0*myzeta**2 - myrho**2) + 5.0*py*myzeta*(-2.0*myzeta**2 + myrho**2))&
+ px**2*(2.0*myzeta**4 + 3.0*myzeta**2*myrho**2 + myrho**4 + 5.0*py*myzeta*(-2.0*myzeta**2 + myrho**2)&
+py**2* (-6.0*myzeta**2 + 3.0*myrho**2)))&
-(3.0*R**5*(px**4*(8.0*myzeta**4 - 24.0*myzeta**2*myrho**2 + 3.0*myrho**4 + 8.0*py*(2.0*myzeta**3&
- 3.0*myzeta*myrho**2)) - 8.0*py**3*(4.0*myzeta**5 + 2.0*myzeta**3*myrho**2 - 2.0*myzeta*myrho**4&
+ py**2*(4.0*myzeta**3 - 6.0*myzeta*myrho**2) - py*(12.0*myzeta**4 - 22.0*myzeta**2* myrho**2&
+ myrho**4)) - 8.0*px**2*py*(-6.0*myzeta**5 - 3.0*myzeta**3* myrho**2 + 3.0*myzeta* myrho**4&
+ py**2* (2.0* myzeta**3 - 3.0* myzeta* myrho**2) + py* (22.0*myzeta**4 - 45.0*myzeta**2* myrho**2&
+ 3.0*myrho**4))))/(16.0*(px**2 + py**2)**(3.5)* (myzeta**2 + myrho**2)**(3.5))

w2IntegrandZFF = w2IntegrandZFF*myrho
endif
end
!====================== Far Field kernel for px^2+py^2 > 4R^2 =====================!





function w2IntegrandR(ellipticK,ellipticE,ellipticK1, ellipticE1)
! this is the horizontal component of the velocity
use globalinfo
implicit none

real (kind=8) w2IntegrandR,ellipticE, ellipticK,ellipticE1, ellipticK1,ILogR
real (kind=8), external :: I1R,I2R,I3R,I4R,I5R,I6R,I7R,I8R,LogTermThetaR

w2IntegrandR=0.0;
ILogR = 0.0;


if (px > 0.0) then



if (abs (myzeta*px + py * myrho ) < logsing) then
!call simp2(LogTermThetaR, real(0.0,kind=8), real(2.0*pi,kind=8),integthres, ILogR)

call trapz1(LogTermThetaR,  real(0.0,kind=8),real(2.0*pi,kind=8), logtrapzbig, ILogR)

else
call trapz1(LogTermThetaR,  real(0.0,kind=8),real(2.0*pi,kind=8), logtrapz, ILogR)
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
use globalinfo
implicit none

real (kind=8) w2IntegrandZ,ellipticE, ellipticK,ellipticE1, ellipticK1,ILogZ
real (kind=8), external :: I1Z,I2Z,I3Z,I4Z,I5Z,I6Z,I7Z,I8Z,LogTermThetaZ


ILogZ = 0.0;

if (abs (myzeta*px + py * myrho ) < logsing .and. abs (myzeta*px + py * myrho ) > 0.0) then


call trapz1(LogTermThetaZ,real(0.0,kind=8),real(2.0*pi,kind=8),logtrapzbig, ILogZ)


!call simp2(LogTermThetaZ, real(0.0,kind=8), real(2.0*pi,kind=8),integthres, ILogZ)
elseif (abs (myzeta*px + py * myrho ) > logsing )then

call trapz1(LogTermThetaZ, real(0.0,kind=8),real(2.0*pi,kind=8),logtrapz, ILogZ)

!elseif (rho=0.0 .and. px =0.0 ) then

!ILogZ=2*pi (-(R**2 - px**2  - py**2)*(R**2 - myzeta**2))/(2.0*(myzeta**2))&
!*((-6*myzeta**4)/ &
!(R*(-(myzeta*py) + R**2)
!(myzeta + abs(myzeta))*&
!abs(-(myzeta*py) + R**2)))

endif




w2IntegrandZ = I1Z(ellipticK,ellipticE)+I2Z(ellipticK1,ellipticE1)&
+I3Z(ellipticK1,ellipticE1)+I4Z(ellipticK1,ellipticE1)&
+I5Z(ellipticK1,ellipticE1)+I6Z(ellipticK1,ellipticE1) &
+I7Z(ellipticK1,ellipticE1)+I8Z(ellipticK1,ellipticE1)+ILogZ

w2IntegrandZ = w2IntegrandZ*myrho/pi

!print *, "w2IntegrandZ", w2IntegrandZ


end



function I1R(ellipticK,ellipticE)
use globalinfo
implicit none

Real*8 I1R, ellipticE, ellipticK

I1R=0
if ( ((py-myzeta)**2+(px-myrho)**2)>singthres) then


I1R= 2.0*px**(-1.0)*(px**2.0+(-2.0)*px*myrho+myrho**2.0+(py+(-1.0)*myzeta)**2.0)**(-1.0/2.0)* &
(px**2.0+2.0*px*myrho+myrho**2.0+(py+(-1.0)*myzeta)**2.0)**(-1.0)*(py+(-1.0)*myzeta)*( &
(px**2.0+(-1.0)*myrho**2.0+(-1.0)*(py+(-1.0)*myzeta)**2.0)*ellipticE+(px**2.0+2.0* &
px*myrho+myrho**2.0+(py+(-1.0)*myzeta)**2.0)*ellipticK)

endif

end




function I1Z(ellipticK,ellipticE)
use globalinfo
implicit none

Real*8 I1Z, ellipticE, ellipticK

I1Z=0

if ( ((py-myzeta)**2+(px-myrho)**2)>singthres) then


I1Z=4.0*(py+(-1.0)*myzeta)**2.0*(px**2.0+(-2.0)*px*myrho+myrho**2.0+((-1.0)*py+myzeta)**2.0) &
**(-1.0/2.0)*(px**2.0+2.0*px*myrho+myrho**2.0+((-1.0)*py+myzeta)**2.0)**(-1.0)* &
ellipticE+4.0*(px**2.0+(-2.0)*px*myrho+myrho**2.0+((-1.0)*py+myzeta)**2.0)**( &
-1.0/2.0)*ellipticK;
endif

end




function I2R (ellipticK1, ellipticE1)
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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
use globalinfo
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

!if (myt>896) then
!print *, "theta",theta, "px", px, "py", py, "rho", myrho, "zeta",myzeta,"integrand",LogTermThetaZ
!endif

return
end










