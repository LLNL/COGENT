import glob, os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import tempfile
from plot_blocks import plot_all_blocks, plot_valid_block, plot_extended_block, get_block_names

script_location = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
grid_generator_path = script_location + '/grid_generator'
smooth_ghosts_path = script_location + '/smooth_ghosts'

######### USER-SPECIFIED INPUT ###################################

psiNormInner = 0.89                              # normalized flux at the inner core boundary
psiNormOuter = 1.11                             # normalized flux at the outer sol boundary
n_psiCORE = 22                                  # number of radial cells in the core region
n_thetaCORE = 128                                # number of poloidal cells in the core region
n_thetaInnerLeg = 6                             # number of poloidal cells in the inner leg
n_thetaOuterLeg = 6                             # number of poloidal cells in the outer leg
n_psiPF = 6                                    # number of radial cells in the PF region
n_thetaMCORE = 0                               # number of poloidal cells in the mcore region (set to 0 for original 8 blocks)

n_extrpInnerLeg = 0                             # number of poloidal extrapolation cells in the inner leg
n_extrpOuterLeg = 2                             # number of poloidal extrapolation cells in the outer leg

n_psiGhost = 2                                  # number of radial ghost cells
n_thetaGhost = 2                                # number of poloidal ghost cells

trans_rad = -0.2                                 # transition radius (set to a negative value for flux-aligned grids)
trans_length = 0.6                              # transition length = distance from X point beyond which grid is locally orthogonal
blending_factor = 0.6                           # blending factor weighting the real and block-aligned flux

init_guessXr = 1.4                              # initial (r) guess for Xpt
init_guessXz = -1.05                            # initial (z) guess for Xpt
init_guessOr = 1.6                              # initial (r) guess for Opt
init_guessOz = 0.0                              # initial (z) guess for Opt

compute_critical_points = 1                     # should we compute critical points (if not, rely on the data in DCT file).
output_gridUE = 0                               # should we output data for gridUE generator?

run_GridGenerator = True                        # should we run grid generator (or only do the plotting of existing output)?
smooth_Ghosts = False                           # should we extrapolate smooth ghost cells?
plot_block_valid = True                         # should we plot grid blocks with physical ghosts?
plot_block_ghost = True                         # should we plot grid blocks with all ghosts?

############## ADD SMOOTH GHOST CELLS ######################################

def smoothGhosts(input_file_name, mapping_file_name):

   if (len(glob.glob(smooth_ghosts_path + "/*.ex")) > 1):
      print "Multiple smoothGhosts executables are detected, the user should manually set executable_name"
      sys.exit()
   else:
      executable_name = os.path.basename(glob.glob(smooth_ghosts_path + "/*.ex")[0])

   args = ['cd ' + smooth_ghosts_path + ';./' + executable_name, input_file_name, mapping_file_name]

   str_args = [ str(x) for x in args ]

   subprocess.call(' '.join(str_args), shell=True)

############## RUN GRID-GENERATOR ######################################

if run_GridGenerator:

   if not (len(sys.argv) == 3):
      print 'Usage: python pyRun.py <DCT_file> <output_file>'
      sys.exit()
   else:      
      DCT_file_name = os.getcwd() + '/' + sys.argv[1]
      mapping_file_name = os.getcwd() + '/' + sys.argv[2]

   if (len(glob.glob(grid_generator_path + "/*.ex")) > 1):
      print "Multiple gridGenerator executables are detected, the user should manually set executable_name"
      sys.exit()
   else:
      executable_name = os.path.basename(glob.glob(grid_generator_path + "/*.ex")[0])

   if smooth_Ghosts:
      fd_tmp,grid_file = tempfile.mkstemp(text=True)
   else:
      grid_file = mapping_file_name

   #replace ';./' with '; mpirun -np 2 ./' for parallel usage
   args = ['cd ' + grid_generator_path + ';./' + executable_name, DCT_file_name, grid_file, psiNormInner,
            psiNormOuter, n_psiCORE, n_thetaCORE, n_thetaInnerLeg, n_thetaOuterLeg, n_psiPF, n_thetaMCORE,
            n_extrpInnerLeg, n_extrpOuterLeg, n_psiGhost, n_thetaGhost, trans_rad, trans_length, blending_factor,
            init_guessXr, init_guessXz, init_guessOr, init_guessOz, compute_critical_points, output_gridUE]

   str_args = [ str(x) for x in args ]

   subprocess.call(' '.join(str_args), shell=True)

   if smooth_Ghosts:
      smoothGhosts(grid_file, mapping_file_name)

############  FUNCTION TO READ GRID-BLOCK DATA  #####################

#def Read_Two_Column_File(file_name):
#   with open(file_name, 'r') as data:
#      first_line = data.readline()
#      x = []
#      y = []
#      for line in data:
#         p = line.split()
#         x.append(float(p[0]))
#         y.append(float(p[1]))
#
#   p = first_line.split()
#
#   dim0 = int(p[0])
#   dim1 = int(p[1])
#
#   dataR = np.reshape(x, (dim0,dim1), order="F")
#   dataZ = np.reshape(y, (dim0,dim1), order="F")
#
#   return dataR, dataZ, dim0, dim1

##########  PLOT ENTIRE GRID (INCLUDES PHYSICAL GHOSTS) ##############################

print "PLOTTING ENTIRE GRID"

mapping_file_name = os.getcwd() + '/' + sys.argv[2]

plot_all_blocks(mapping_file_name)

#for file in glob.glob(grid_generator_path + "/output/coords*"):
#   dataR, dataZ, dim0, dim1 = Read_Two_Column_File(file)
#
#   for i in range(dim0):
#      plt.plot(dataR[i,:], dataZ[i,:], color='b', linewidth=0.2)
#
#   for i in range(dim1):
#      plt.plot(dataR[:,i], dataZ[:,i], color='b', linewidth=0.2)
#

plt.axes().set_aspect('equal')
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.savefig(grid_generator_path + '/output/EntireGrid.pdf')
plt.close('all')

#########  PLOT GRID BLOCKS WITH PHYSICAL GHOSTS ONLY  ##############################

print "PLOTTING GRID BLOCKS WITH PHYSICAL GHOSTS ONLY"

if plot_block_valid:
   
   block_names = get_block_names(mapping_file_name)

   if (len(block_names) == 8):
      Ncols = 4
      Nrows = 2
   elif (len(block_names) == 10):
      Ncols = 5
      Nrows = 2
   else:
      print 'Number of blocks is not 8 or 10'
      sys.exit()

   fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols)

   for ax, blockID in zip(axes.flat, block_names):
      ax.set_title(blockID.upper())
      plot_valid_block(mapping_file_name, blockID, ax)

   fig.set_tight_layout(True)
   fig.savefig(grid_generator_path + '/output/GridBlocksWithPhysGhosts')
   plt.close('all')

#if plot_block_valid:
#   Ncols = 1
#   Nrows = 1
#   blocks = []
#   
#   if (len(glob.glob(grid_generator_path + "/output/coords*")) == 8):
#      Ncols = 4
#      Nrows = 2
#      blocks = ['LCORE', 'RCORE', 'LCSOL', 'RCSOL', 'LSOL', 'RSOL', 'LPF', 'RPF']
#   else:
#      Ncols = 5
#      Nrows = 2
#      blocks = ['LCORE', 'RCORE', 'LCSOL', 'RCSOL', 'LSOL', 'RSOL', 'LPF', 'RPF', 'MCORE', 'MCSOL']
#
#   fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols)
#
#   for ax, file, blockID in zip(axes.flat, glob.glob(grid_generator_path + "/output/coords*"), blocks):
#
#      dataR, dataZ, dim0, dim1 = Read_Two_Column_File(file)
#
#      ax.set_title(blockID)
#
#      start = np.min(dataR)
#      end = np.max(dataR)
#      tick_spacing = (end - start)/2.001 #increase over 2.0 to include left-end tick
#      ax.xaxis.set_ticks(np.arange(start, end, tick_spacing))
#      ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
#
#      for i in range(dim0):
#         ax.plot(dataR[i,:], dataZ[i,:], color='b', linewidth=0.2)
#   
#      for i in range(dim1):
#         ax.plot(dataR[:,i], dataZ[:,i], color='b', linewidth=0.2)
#
#   fig.set_tight_layout(True)
#   fig.savefig(grid_generator_path + '/output/GridBlocksWithPhysGhosts')
#   plt.close('all')

#########  PLOT GRID BLOCKS WITH ALL GHOSTS  ##############################

print "PLOTTING GRID BLOCKS WITH ALL GHOSTS "

if plot_block_ghost:
   
   block_names = get_block_names(mapping_file_name)

   if (len(block_names) == 8):
      Ncols = 4
      Nrows = 2
   elif (len(block_names) == 10):
      Ncols = 5
      Nrows = 2
   else:
      print 'Number of blocks is not 8 or 10'
      sys.exit()

   fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols)

   for ax, blockID in zip(axes.flat, block_names):
      ax.set_title(blockID.upper())
      plot_extended_block(mapping_file_name, blockID, ax)

   fig.set_tight_layout(True)
   fig.savefig(grid_generator_path + '/output/GridBlocksWithAllGhosts')
   plt.close('all')

#if plot_block_ghost:
#   Ncols = 1
#   Nrows = 1
#   blocks = []
#   
#   if (len(glob.glob(grid_generator_path + "/output/extended_coords*")) == 8):
#      Ncols = 4
#      Nrows = 2
#      blocks = ['LCORE', 'RCORE', 'LCSOL', 'RCSOL', 'LSOL', 'RSOL', 'LPF', 'RPF']
#   else:
#      Ncols = 5
#      Nrows = 2
#      blocks = ['LCORE', 'RCORE', 'LCSOL', 'RCSOL', 'LSOL', 'RSOL', 'LPF', 'RPF', 'MCORE', 'MCSOL']
#
#   fig, axes = plt.subplots(nrows=Nrows, ncols=Ncols)
#
#   for ax, file, blockID in zip(axes.flat, glob.glob(grid_generator_path + "/output/extended_coords*"), blocks):
#   
#      dataR, dataZ, dim0, dim1 = Read_Two_Column_File(file)
#      
#      ax.set_title(blockID)
#      
#      start = np.min(dataR)
#      end = np.max(dataR)
#      tick_spacing = (end - start)/2.001 #increase over 2.0 to include left-end tick
#      ax.xaxis.set_ticks(np.arange(start, end, tick_spacing))
#      ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
#      
#      for i in range(dim0):
#         ax.plot(dataR[i,:], dataZ[i,:], color='b', linewidth=0.2)
#      
#      for i in range(dim1):
#         ax.plot(dataR[:,i], dataZ[:,i], color='b', linewidth=0.2)
#
#   fig.set_tight_layout(True)
#   fig.savefig(grid_generator_path + '/output/GridBlocksWithAllGhosts')
#   plt.close('all')

