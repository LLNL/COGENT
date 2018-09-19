read ~/utils/basis/uedge_step.f



#import single-null grid from gridue file
echo=no

# initialize packages --
package flx;package grd;package aph;package bbb;package com

# mesh construction --
# Set some parameters that usually get set by subroutine flxrun
gengrid=0
mhdgeo=1
geometry="snull"
nxpt=1

call gallot("Xpoint_indices",0)
call readgrid("gridue",runid)


nx=nxm
ny=nym

#win on
#read plotmesh
#plot zm(ixpt1,,0) rm(ixpt1,,0) mark=circle color=red
#plot zm(ixpt2,,0) rm(ixpt2,,0) mark=circle color=green
#plot zm(,iysptrx1,0) rm(,iysptrx1,0) mark=cross color=red
#attr legend = no
#frame



allocate
tes=10*ev
tis=10*ev
nis=1e20
ngs=1e16
ups=0

##-solving only the temperature diffusion equation
#isteon=1
#istion=0
#isnion=0
#isupon=0
#isupgon=0
#isngon=0


#2#-reset to flat initial conditions-
double te_offset=1e0
tes=(te_offset+100.)*ev; te=tes
#integer i
#do i=0,nx+1
#  tes(i,)=(10.0-9.0*(yyc-min(yyc))/(max(yyc)-min(yyc)))*ev
#enddo
#te=tes



#2a#-try simple homogeneous boundary conditions
iflcore=0; tcoree = te_offset+1.0
istewc= 1; tedge =  te_offset
istepfc=1; tewalli= te_offset; tewallo=te_offset
ibctepl=0; tepltl = te_offset
ibctepr=0; tepltr = te_offset

#2b#-fix parallel heat conduction at fixed value
concap=1
afix=1e-10 ##1e-1 ###5e0 ###50e0


#3#-make sure the right set of equations is used
isnion = 0
istion = 0        
isteon = 1        
isupon = 0 
isngon = 0        
isphion = 0      


##############-try initial tiny step to get things going
read tecond.bas
dtreal=1e-8; exmain

#-reset initial profile
integer i
real yycn(0:ny+1)=(yyc-min(yyc))/(max(yyc)-min(yyc))


#tes=1.0*ev #-make tanh() profile above xpt
#do i=ixpt1+1,ixpt2
#  ###tes(i,)=(10.0-9.0*(yyc-min(yyc))/(max(yyc)-min(yyc)))*ev
#  tes(i,)=(1.0+0.5*(1.0-tanh((yycn-0.5)*10)))*ev
#enddo

#-set step-function
tes=1.0*ev
tes(,0:iysptrx1)=2.0*ev

te=tes

win on
plot te(ixmp,)/ev yyc
attr legend=no
attr labels=no
#######################################################



#4#-prepare time-evolution run
integer istp
integer nstp=10
double tmax=1e-4
double dtnow=iota(nstp)*0 + tmax/(nstp-1) #-uniform time step
###double dtnow=sqrt(iota(nstp))*(3*tmax)/(2*nstp*sqrt(nstp)) #-growing as square root
###double dtnow=iota(nstp)*(2*tmax)/(nstp*(nstp-1)) #-linearly growing


#-allocate memory for time history
double te_stor(nstp,0:nx+1,0:ny+1)
double tim_stor(nstp)

  #-initial state at t=0
  istp=1
  tim_stor(istp)=0.0
  te_stor(istp,,)=te/ev ###-1d0 ##-for retaining accuracy



#5#-start a time evolution run-
restart=1

do istp=2,nstp

  istp  
  uedge_step(dtnow(istp-1),0)
  
  tim_stor(istp)=tim_stor(istp-1)+dtnow(istp-1)
  << "total time=" << tim_stor(istp)

  te_stor(istp,,)=te/ev ###-1d0 ##-for retaining accuracy
  plot te(ixmp,)/ev yyc
  
enddo

integer nmax=nstp

#6#-show time history at a specific point
win on 2
nf
plot te_stor(1:nmax,nx/2,ny/2) tim_stor(1:nmax)
#attr legend=no
#attr labels=no


#7#-save time history and grid data
double rc=rm(,,0)
double zc=zm(,,0)
create uedge_test1.pdb
write rc, zc, te_stor, ev, tim_stor, nmax
close




#restart=1
#isbcwdt=1
#dtreal=1e-6; exmain
#read rdinitdt
#read rdcontdt
#
#isbcwdt=0
#dtreal=1e20
#exmain
