import os
import sys
from plot_blocks import plot_valid_block, plot_extended_block

if not (len(sys.argv) == 2):
   print 'Usage: python pyRun_interactive.py <file_name>'
   sys.exit()
else:      
   file_name = os.getcwd() + '/' + sys.argv[1]

plot_ghosts = raw_input("Do you want to plot ghost cells? ")

include_ghosts = (plot_ghosts == 'yes' or plot_ghosts == 'Yes' or plot_ghosts == 'y' or plot_ghosts == 'Y')

keep_reading = True

while keep_reading:
   block = raw_input("Which block do you want to plot (type 'none' or 'q' to quit)? ")
   if block == 'none' or block == 'q':
      keep_reading = False
   else:
      if include_ghosts:
         plot_extended_block(file_name, block)
      else:
         plot_valid_block(file_name, block)

