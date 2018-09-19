import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoLocator
import numpy as np

mpl.rcParams['lines.linewidth'] = 0.2

def plot_all_blocks(file_name):

   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.set(xlabel='R (m)', ylabel='Z (m)')

   block_names = get_block_names(file_name)

   for block_name in block_names:
      plot_valid_block(file_name, block_name, ax)

   ax.autoscale(enable=True, axis='both', tight=True)
   ax.get_xaxis().set_major_locator(AutoLocator())

   plt.show();


def get_block_names(file_name):

   block_names = [];

   with open(file_name,'r') as f:

      block_start = 0

      for i,line in enumerate(f):
         if i==block_start:

            data = line.split()

            this_block_name = data[0]
            n0              = int(data[1])
            n0_extend       = int(data[2])
            n1              = int(data[3])
            n1_extend       = int(data[4])
            dim0 = n0 + 2*n0_extend;
            dim1 = n1 + 2*n1_extend;
            block_size = dim0 * dim1

            block_start += block_size+1

            block_names.append(this_block_name)

   return block_names

def read_block(file_name, block_name):

   with open(file_name,'r') as f:

      block_start = 0
      block_data = [] 
      found_block = False

      for i,line in enumerate(f):

         if i==block_start:

            data = line.split()

            this_block_name = data[0]
            n0              = int(data[1])
            n0_extend       = int(data[2])
            n1              = int(data[3])
            n1_extend       = int(data[4])
            dim0 = n0 + 2*n0_extend;
            dim1 = n1 + 2*n1_extend;
            block_size = dim0 * dim1

            block_start += block_size+1

            read_block = this_block_name == block_name

            if read_block:
               block_dim = (dim0,dim1)
               block_n = (n0,n1)
               block_n_extend = (n0_extend,n1_extend)
               found_block = True

         elif read_block:

            data = map(float,line.split())
            block_data.append(data)

      if found_block:            
         na = np.asfarray(block_data)
         R = np.reshape(na[:,0],block_dim,order='F')
         Z = np.reshape(na[:,1],block_dim,order='F')

         return R, Z, block_n, block_n_extend
      else:
         print 'The block \'' + block_name + '\' was not found in the file ' + file_name
         raise

def plot_valid_block(file_name, block_name, ax=0):

   block_names = get_block_names(file_name)

   if block_name in block_names:

      R,Z,n,n_extend = read_block(file_name, block_name)

      R_valid = R[n_extend[0]:n_extend[0]+n[0],n_extend[1]:n_extend[1]+n[1]]
      Z_valid = Z[n_extend[0]:n_extend[0]+n[0],n_extend[1]:n_extend[1]+n[1]]

      create_plot = ax==0

      if create_plot:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.set(xlabel='R (m)', ylabel='Z (m)', title=block_name.upper())
         ax.axis('equal')

      start = np.min(R_valid)
      end = np.max(R_valid)
      tick_spacing = (end - start)/2.001 #increase over 2.0 to include left-end tick
      ax.xaxis.set_ticks(np.arange(start, end, tick_spacing))
      ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

      plot_grid(R_valid,Z_valid,ax)

      if create_plot:
         plt.show()
   else:
      print 'Uncognized block name: ', block_name

def plot_extended_block(file_name, block_name, ax=0):

   R,Z,n,n_extend = read_block(file_name, block_name)

   create_plot = ax==0

   if create_plot:
      fig = plt.figure()
      ax = fig.add_subplot(111)
      ax.set(xlabel='R (m)', ylabel='Z (m)', title=block_name.upper())

   start = np.min(R)
   end = np.max(R)
   tick_spacing = (end - start)/2.001 #increase over 2.0 to include left-end tick
   ax.xaxis.set_ticks(np.arange(start, end, tick_spacing))
   ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

   plot_extended_grid(R,Z,n_extend,ax)

   if create_plot:
      plt.show()

def plot_grid(R,Z,ax):
   
   dim = R.shape
   assert Z.shape == dim

   for i in range(dim[0]):
      ax.plot(R[i,:],Z[i,:],'k')

   for j in range(dim[1]):
      ax.plot(R[:,j],Z[:,j],'k')

def plot_extended_grid(R,Z,n_extend,ax):
   
   dim = R.shape
   assert Z.shape == dim

   for i in range(0,n_extend[0]):
      ax.plot(R[i,:],Z[i,:],'r')
   for i in range(n_extend[0],dim[0]-n_extend[0]):
      ax.plot(R[i,0:n_extend[1]+1],Z[i,0:n_extend[1]+1],'r')
      ax.plot(R[i,n_extend[1]:dim[1]-n_extend[1]],Z[i,n_extend[1]:dim[1]-n_extend[1]],'k')
      ax.plot(R[i,dim[1]-n_extend[1]-1:dim[1]],Z[i,dim[1]-n_extend[1]-1:dim[1]],'r')
   for i in range(dim[0]-n_extend[0],dim[0]):
      ax.plot(R[i,:],Z[i,:],'r')

   for j in range(0,n_extend[1]):
      ax.plot(R[:,j],Z[:,j],'r')
   for j in range(n_extend[1],dim[1]-n_extend[1]):
      ax.plot(R[0:n_extend[0]+1,j],Z[0:n_extend[0]+1,j],'r')
      ax.plot(R[n_extend[0]:dim[0]-n_extend[0],j],Z[n_extend[0]:dim[0]-n_extend[0],j],'k')
      ax.plot(R[dim[0]-n_extend[0]-1:dim[0],j],Z[dim[0]-n_extend[0]-1:dim[0],j],'r')
   for j in range(dim[1]-n_extend[1],dim[1]):
      ax.plot(R[:,j],Z[:,j],'r')
