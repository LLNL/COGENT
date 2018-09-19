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

win on
read plotmesh
plot zm(ixpt1,,0) rm(ixpt1,,0) mark=circle color=red
plot zm(ixpt2,,0) rm(ixpt2,,0) mark=circle color=green
plot zm(,iysptrx1,0) rm(,iysptrx1,0) mark=cross color=red
attr legend = no
frame



allocate
tes=10*ev
tis=10*ev
nis=1e20
ngs=1e16
ups=0


isteon=1
istion=1
isnion=1
isupon=1
isupgon=0
isngon=0

restart=1
isbcwdt=1
dtreal=1e-6; exmain
read rdinitdt
read rdcontdt

isbcwdt=0
dtreal=1e20
exmain
