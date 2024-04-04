'''
Python script to create plots from the solution of a
HyPar simulation.
'''

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def plotSolution( a_ndims:  int,
                  a_nvars:  int,
                  a_size:   np.ndarray,
                  a_time:   float,
                  a_x:      np.ndarray,
                  a_U:      np.ndarray,
                  a_fname:  str ):

  plt_dir_name='plots'
  if not os.path.exists(plt_dir_name):
        os.makedirs(plt_dir_name)

  if a_ndims == 1 :

    font = {'size':22}
    matplotlib.rc('font', **font)
    figsize=(15,9)

    for var in range(a_nvars):
      try:
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()
        ax.set( xlim=(np.min(a_x), np.max(a_x)),
                ylim=(np.min(a_U[var::a_nvars]),
                      np.max(a_U[var::a_nvars]) ) )
        ax.plot(a_x, a_U[var::a_nvars], lw=2)
        ax.set_title('var {:}, t={:.3}'.format(var,a_time))
        plt.grid(visible=True, linestyle=':', linewidth=1)
        if a_nvars == 1 :
          plt_fname = plt_dir_name+'/'+a_fname
        else:
          plt_fname = plt_dir_name+'/var'+f'{var:02d}'+'_'+a_fname
        print('Saving %s' % plt_fname)
        plt.savefig(plt_fname)
        plt.close()
      except:
        print('Error in plotBinarySolution: unable to generate plot for var %d' % var)
        print('Error in plotBinarySolution: nvars is %d' % a_nvars)

  elif a_ndims == 2 :

    font = {'size':22}
    matplotlib.rc('font', **font)
    colormap='jet'
    figsize=(12,10)

    x = a_x[:a_size[0]]
    y = a_x[a_size[0]:]
    y2d, x2d = np.meshgrid(y, x)

    for var in range(a_nvars):
      try:
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()
        ax.set( xlim=(np.min(x), np.max(x)),
                ylim=(np.min(y), np.max(y)) )
        sol2d = np.transpose(a_U.reshape(a_size[1],a_size[0],a_nvars))
        plot = ax.pcolor(x2d, y2d, sol2d[var,:,:], cmap=colormap)
        ax.set_title('var {:}, t={:.3}'.format(var,a_time))
        fig.colorbar(plot, ax=ax)
        if a_nvars == 1 :
          plt_fname = plt_dir_name+'/'+a_fname
        else:
          plt_fname = plt_dir_name+'/var'+f'{var:02d}'+'_'+a_fname
        print('Saving %s' % plt_fname)
        plt.savefig(plt_fname)
        plt.close()
      except:
        print('Error in plotBinarySolution: unable to generate plot for var %d' % var)
        print('Error in plotBinarySolution: nvars is %d' % a_nvars)

  else :

    print('plotSolution.py: No plotting implemented for a_ndims=%d', a_ndims)

