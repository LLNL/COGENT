from math import *

# Number of radial points defining the physical block
nr = 17

# Number of points defining the radial block extension
nr_extend = 12

# Number of poloidal points defining the physical block
np = 65

# Number of points defining the poloidal block extension
np_extend = 12

# Minor radius minimum
r_min = 0.2

# Minor radius maximum
r_max = 0.6

# Major radius
R0 = 1.6

# Triangularity factor (set to 0 for no triangularity)
#delta = 0.416
delta = 0.

# Elongation factor (set to 1 for circular annulus)
#kappa = 1.66
kappa = 1.

dr = (r_max - r_min) / (nr - 1)
dp = 2.*pi / (np - 1)

#text_file = open("miller_mapping", "w")
text_file = open("miller_mapping_annulus", "w")

text_file.write(format("%d %d %d %d\n"%(nr, nr_extend, np, np_extend)))

for j in range(-np_extend,np+np_extend):
  theta = j*dp

  for i in range(-nr_extend,nr+nr_extend):
    r = i*dr + r_min

    R = R0 + r*cos(theta + asin(delta)*sin(theta))
    Z = kappa*r*sin(theta)

    text_file.write(format("%20.16e %20.16e\n"%(R,Z)))

text_file.close()
