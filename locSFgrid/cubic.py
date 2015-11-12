# plot cubic functions for various eta, xi
from numpy import *
from gist import *
zm0 = .2;zmax = 2.0;Rmaj = 10.
def getcoefs(eta,xi):
    global alpha, beta
    alpha = xi/(1.-xi)
    beta = (zm0/zmax)**3
    a = 1.-3.*alpha**2
    b = 6.*alpha**2
    c = -3.*alpha**2
    d = -eta*beta
    return a, b, c, d

def cubicarr(eta,xi,x):
    a,b,c,d = getcoefs(eta,xi)
    return a*x**3 + b*x**2 + c*x + d

def plot (eta, xi):
    xmax = zmax/sqrt(3.)
    xarr = xmax*arange(101)/100.
    disc = cubicarr(eta,xi, xarr)
    plg(disc,xarr)

def guess(eta, xi):
    a,b,c,d = getcoefs(eta,xi)
    if (abs(a) > 1.e-6):
        myguess = (-b + sqrt(b**2-4*a*c))/(2.*a)
    else:
        myguess = -c/b
    return myguess

def cubicprime(eta,xi,x):
    a,b,c,d = getcoefs(eta,xi)
    return 3.*a*x**2 + 2.*b*x + c

def guessprime(eta,xi):
    myguess = guess(eta,xi)
    deriv = cubicprime(eta,xi,myguess)
    print "myguess,deriv = ", myguess, deriv

def findroot(eta,xi,myguess=None):
    a,b,c,d = getcoefs(eta,xi)
    if (myguess == None):
        myguess = guess(eta,xi)
    # start newton with myguess
    y = myguess
    niter = 0
    resid = abs(cubicarr(eta,xi,y))
    while (abs(resid > 1.e-12) and niter < 40):
        f = cubicarr(eta,xi,y)
        fprime = cubicprime(eta,xi,y)
        y = y - f/fprime
        niter += 1
        resid = abs(f)
        # print "y, resid = ", y, resid
    z = y*zmax
    x = alpha*(zmax-z)
    return eta,xi,a,b,c,d,myguess,y,niter
    # return eta,xi,x,z, niter


# UNIT VECTORS
def unvecs(x,z):
    denomi = 1./(Rmaj+x)
    Bx = -(-3.*x**2+3.*z**2)*denomi
    Bz = -6.*x*z*denomi
    bmag = sqrt(Bx**2 + Bz**2)
    bx = Bx/bmag
    bz = Bz/bmag
    return bx, bz

def getB(x,z):
    denomi = 1./(Rmaj+x)
    Bx = -(-3.*x**2+3.*z**2)*denomi
    Bz = -6.*x*z*denomi
    return Bx, Bz
