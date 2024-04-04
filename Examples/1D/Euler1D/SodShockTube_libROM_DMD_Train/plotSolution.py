'''
Python script to create an animation from the unsteady solution of a
HyPar simulation.
- op_overwrite must be set to "no"
- solution must be in binary format
'''
import os
import time
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import numpy.linalg as LA

hypar_dir = os.environ.get('HYPAR_DIR')
hypar_dir_python = hypar_dir + '/Examples/Python'
sys.path.append(hypar_dir_python)

import modHyParUtils as hyparutils

font = {'size':22}
matplotlib.rc('font', **font)

figsize=(15,9)
plt_dir_name='plots'

fom_plt_style = '-'
rom_plt_style = 'o:'

'''
Set up the simulation parameters
'''
sim_path = '.'
sim_inp_data = hyparutils.readHyParInpFile(sim_path+'/simulation.inp')
solver_inp_data = hyparutils.readHyParInpFile(sim_path+'/solver.inp')

if not sim_inp_data:
    nsims = 1
else:
    nsims = int(sim_inp_data['nsims'][0])

print('Number of simulations: ', nsims)

if solver_inp_data['op_file_format'] != 'binary':
    raise Exception("op_file_format must be 'binary' in solver.inp.")

ndims = int(solver_inp_data['ndims'][0])
nvars = int(solver_inp_data['nvars'][0])
size = np.int32(solver_inp_data['size'])

if not os.path.exists(plt_dir_name):
      os.makedirs(plt_dir_name)

if solver_inp_data['op_overwrite'] == 'no':

  niter = int(solver_inp_data['n_iter'][0])
  dt = float(solver_inp_data['dt'][0])
  t_final = dt*niter

  op_write_iter = int(solver_inp_data['file_op_iter'][0])
  dt_snapshots = op_write_iter*dt
  n_snapshots = int(niter/op_write_iter) + 1

  print('Simulation parameters:')
  print('  ndims = ', ndims)
  print('  nvars = ', nvars)
  print('  grid size = ', size)
  print('  simulation dt = ', dt)
  print('  niter = ', niter)
  print('  final time = ', t_final)
  print('  snapshot dt = ', dt_snapshots)
  print('  number of snapshots = ', n_snapshots)

  '''
  Load simulation data (solution snapshots)
  '''
  x,sol_fom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                          nsims,
                                                          n_snapshots,
                                                          ndims,
                                                          nvars,
                                                          size )
  sol_fom_snapshots = np.float32(sol_fom_snapshots)

  x,sol_rom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                          nsims,
                                                          n_snapshots,
                                                          ndims,
                                                          nvars,
                                                          size,
                                                          op_root='op_rom')
  sol_rom_snapshots = np.float32(sol_rom_snapshots)

  for var in range(nvars):
    for i in range(n_snapshots):
      for s in range(nsims):
        solution_fom_snapshots_sim = sol_fom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
        solution_rom_snapshots_sim = sol_rom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
        diff_fom_rom = solution_fom_snapshots_sim - solution_rom_snapshots_sim
        fom_solution_norm = LA.norm(solution_fom_snapshots_sim[i,var::nvars])
        err_norm = LA.norm(diff_fom_rom[i,var::nvars])
        if fom_solution_norm > 1e-14:
          err_norm = err_norm / fom_solution_norm
        # Plot solutions
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()
        ax.set( xlim=(np.min(x), np.max(x)),
                ylim=(np.min(sol_rom_snapshots[i,var::nvars]),
                      np.max(sol_rom_snapshots[i,var::nvars]) ) )
        ax.plot(x, solution_fom_snapshots_sim[i, var::nvars], fom_plt_style,lw=2, label='FOM')
        ax.plot(x, solution_rom_snapshots_sim[i, var::nvars], rom_plt_style,lw=2, label='ROM')
        ax.legend()
        ax.set_title('var {:}, t={:.3},  diff between ROM and FOM={:.3e}'.format(var,i*dt_snapshots,err_norm))
        plt.grid(visible=True, linestyle=':', linewidth=1)
        if nsims > 1:
          plt_fname = plt_dir_name+'/fig_'+f'{var:02d}'+'_'+f'{s:04d}'+'_'+f'{i:05d}'+'.png'
        else:
          plt_fname = plt_dir_name+'/fig_'+f'{var:02d}'+'_'+f'{i:05d}'+'.png'
        print('Saving %s' % plt_fname)
        plt.savefig(plt_fname)
        # Plot diff
        fig = plt.figure(figsize=figsize)
        ax = plt.axes()
        ax.set( xlim=(np.min(x), np.max(x)),
                ylim=(np.min(diff_fom_rom[:,var::nvars]),
                      np.max(diff_fom_rom[:,var::nvars]) ) )
        ax.plot(x, diff_fom_rom[i, var::nvars], lw=2)
        ax.set_title('var {:}, t={:.3}'.format(var,i*dt_snapshots))
        plt.grid(visible=True, linestyle=':', linewidth=1)
        if nsims > 1:
          plt_fname = plt_dir_name+'/fig_diff_'+f'{var:02d}'+'_'+f'{s:04d}'+'_'+f'{i:05d}'+'.png'
        else:
          plt_fname = plt_dir_name+'/fig_diff_'+f'{var:02d}'+'_'+f'{i:05d}'+'.png'
        print('Saving %s' % plt_fname)
        plt.savefig(plt_fname)

else:

  niter = int(solver_inp_data['n_iter'][0])
  dt = float(solver_inp_data['dt'][0])
  t_final = dt*niter

  n_snapshots = 1

  print('Simulation parameters:')
  print('  ndims = ', ndims)
  print('  nvars = ', nvars)
  print('  grid size = ', size)
  print('  final time = ', t_final)
  print('  number of snapshots = ', n_snapshots)

  '''
  Load simulation data (solution snapshots)
  '''
  x,sol_fom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                          nsims,
                                                          n_snapshots,
                                                          ndims,
                                                          nvars,
                                                          size )
  sol_fom_snapshots = np.float32(sol_fom_snapshots)

  x,sol_rom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                          nsims,
                                                          n_snapshots,
                                                          ndims,
                                                          nvars,
                                                          size,
                                                          op_root='op_rom')
  sol_rom_snapshots = np.float32(sol_rom_snapshots)

  for var in range(nvars):
    for s in range(nsims):
      solution_fom_snapshots_sim = sol_fom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
      solution_rom_snapshots_sim = sol_rom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
      diff_fom_rom = solution_fom_snapshots_sim - solution_rom_snapshots_sim
      fom_solution_norm = LA.norm(solution_fom_snapshots_sim[0,var::nvars])
      err_norm = LA.norm(diff_fom_rom[0,var::nvars])
      if fom_solution_norm > 1e-14:
        err_norm = err_norm / fom_solution_norm
      # Plot solutions
      fig = plt.figure(figsize=figsize)
      ax = plt.axes()
      ax.set( xlim=(np.min(x), np.max(x)),
              ylim=(np.min(sol_rom_snapshots[0,var::nvars]),
                    np.max(sol_rom_snapshots[0,var::nvars]) ) )
      ax.plot(x, solution_fom_snapshots_sim[0, var::nvars], fom_plt_style,lw=2, label='FOM')
      ax.plot(x, solution_rom_snapshots_sim[0, var::nvars], rom_plt_style,lw=2, label='ROM')
      ax.legend()
      ax.set_title('var {:}, t={:.3}, diff between ROM and FOM={:.3}'.format(var,t_final,err_norm))
      plt.grid(visible=True, linestyle=':', linewidth=1)
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_'+f'{var:02d}'+'_'+f'{s:04d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_'+f'{var:02d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      # Plot diff
      fig = plt.figure(figsize=figsize)
      ax = plt.axes()
      ax.set( xlim=(np.min(x), np.max(x)),
              ylim=(np.min(diff_fom_rom[:,var::nvars]),
                    np.max(diff_fom_rom[:,var::nvars]) ) )
      ax.plot(x, diff_fom_rom[0, var::nvars], lw=2)
      ax.set_title('var {:}, t={:.3}'.format(var,t_final))
      plt.grid(visible=True, linestyle=':', linewidth=1)
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_diff_'+f'{var:02d}'+'_'+f'{s:04d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_diff_'+f'{var:02d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)


