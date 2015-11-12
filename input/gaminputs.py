# script to convert Xu's specifications to ESLCode inptus
from numpy import *
A=2.0;
rhocoef=sqrt(2.)*1.02
R0=1.71
roverR = 0.2
T=3000.
B0=15.
kperprho = 0.1375
q=3.0
epsilon = roverR
emaxoverT = 6.0

# assume Xu intends kperp to be defined in terms of rho at B=B0
B0gauss = B0*1.e4
rho = rhocoef*sqrt(A*T)/B0gauss
print "rho(m) = ", rho
kperp = kperprho/rho
wavelen = 2.*pi/kperp
print "wavelength = ", wavelen
rbar = roverR*R0
print "rmin = ", rbar-wavelen/2.
print "rmax = ", rbar+wavelen/2.
print "RB = ", B0*R0
#Bpol at R=R0
Bpol0 = epsilon*B0/q
print "Bpol0 = ", Bpol0
dpsidr = R0*Bpol0
print "dpsidr = ", dpsidr

vparmax = sqrt(emaxoverT/(0.5*A))
mumax = 2.0*emaxoverT/B0
print "vparmax, mumax = ", vparmax,mumax
