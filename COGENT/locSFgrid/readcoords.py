from numpy import *
from gist import *
class gridreader:
    def __init__(self):
        pass
    def readcoords(self,blockno):
        fcoord = open("coords"+`blockno`,"r")
        line1 = fcoord.readline()
        line1split = line1.split()
        self.nrad = eval(line1split[0])
        self.npol = eval(line1split[1])
        ndat = self.nrad*self.npol
        self.coords = zeros([ndat,2],"d")
        for irad in range(ndat):
            linesplit = fcoord.readline().split()
            self.coords[irad,0] = eval(linesplit[0])
            self.coords[irad,1] = eval(linesplit[1])
        fcoord.close()

    def plotthetas(self,color="black"):
        # plots the "xi" coordinate lines
        xarr = zeros(self.npol,"d")
        zarr = zeros(self.npol,"d")
        index = 0
        for irad in range(self.nrad):
            for ipol in range(self.npol):
                xarr[ipol] = self.coords[index,0]
                zarr[ipol] = self.coords[index,1]
                index += 1
            plg(zarr,xarr,marks=0,color=color)
    def plotflines(self,iskip=1,color="black"):
        # plots the flux surfaces
        xarr = zeros(self.npol,"d")
        zarr = zeros(self.npol,"d")
        for irad in range(0,self.nrad,iskip):
            index = irad
            for ipol in range(self.npol):
                xarr[ipol] = self.coords[index,0]
                zarr[ipol] = self.coords[index,1]
                index += self.nrad
            plg(zarr,xarr,marks=0,color=color)

    def plotone(self,irad):
        # plots just one of the curves
        xarr = zeros(self.npol,"d")
        zarr = zeros(self.npol,"d")
        index = irad*self.npol
        for ipol in range(self.npol):
            xarr[ipol] = self.coords[index,0]
            zarr[ipol] = self.coords[index,1]
            index += 1
        plg(zarr,xarr)

r = gridreader()
r.readcoords(0)

