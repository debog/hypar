'''
Python script to create streamline plots from the solution
of a HyPar simulation. This particular file is for the
2D Navier-Stokes/Euler physics, where the solution vector
components are (rho, rho*u, rho*v, e), but it can be modified
for other physics by appropriately computing the Cartesian
velocity components from the HyPar solution.

- if op_overwrite is set to "no", a plot is generated
for each variable (solution vector component) and each
simulation time for which the solution is available.
- if op_overwrite is set to "yes", a single plot is
created for for each variable (solution vector component).
- solution must be in binary format

Make sure the environment variable "HYPAR_DIR" is set
and points to the correct location (/path/to/hypar)
'''

import os
import time
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

hypar_dir = os.environ.get('HYPAR_DIR')
hypar_dir_python = hypar_dir + '/Examples/Python'
sys.path.append(hypar_dir_python)

import modHyParUtils as hyparutils

font = {'size':22}
matplotlib.rc('font', **font)
colormap='jet'

plt_dir_name='plots'
fig_nv=1
fig_nh=3
figsize=(1.4*fig_nh*10,fig_nv*10)

density = 4
linewidth = 1

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
  grid, soln_fom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                              nsims,
                                                              n_snapshots,
                                                              ndims,
                                                              nvars,
                                                              size )
  grid, soln_rom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                              nsims,
                                                              n_snapshots,
                                                              ndims,
                                                              nvars,
                                                              size,
                                                              op_root='op_rom')
  soln_fom_snapshots = np.float32(soln_fom_snapshots)
  soln_rom_snapshots = np.float32(soln_rom_snapshots)
  x = grid[:size[0]]
  y = grid[size[0]:]
  print('2D Domain:');
  print(' x: ', np.min(x), np.max(x))
  print(' x.shape: ', x.shape)
  print(' y: ', np.min(y), np.max(y))
  print(' y.shape: ', y.shape)
  y2d, x2d = np.meshgrid(y, x)

  for i in range(n_snapshots):
    for s in range(nsims):
      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      soln_fom_snapshots_sim = soln_fom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
      sol2d_fom = np.transpose(soln_fom_snapshots_sim.reshape(n_snapshots,size[1],size[0],nvars))[:,:,:,i]
      sol2d_fom_u = sol2d_fom[1,:,:] / sol2d_fom[0,:,:]
      sol2d_fom_v = sol2d_fom[2,:,:] / sol2d_fom[0,:,:]
      sol2d_fom_speed = np.sqrt(np.square(sol2d_fom_u)+np.square(sol2d_fom_v))
      soln_rom_snapshots_sim = soln_rom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
      sol2d_rom = np.transpose(soln_rom_snapshots_sim.reshape(n_snapshots,size[1],size[0],nvars))[:,:,:,i]
      sol2d_rom_u = sol2d_rom[1,:,:] / sol2d_rom[0,:,:]
      sol2d_rom_v = sol2d_rom[2,:,:] / sol2d_rom[0,:,:]
      sol2d_rom_speed = np.sqrt(np.square(sol2d_rom_u)+np.square(sol2d_rom_v))
      diff_speed = sol2d_fom_speed - sol2d_rom_speed
      plot_fom = axes[0].streamplot(x2d.T, y2d.T, sol2d_fom_u.T, sol2d_fom_v.T,
                                    color=sol2d_fom_speed.T,
                                    linewidth=linewidth,
                                    cmap=colormap,
                                    density=density)
      plot_rom = axes[1].streamplot(x2d.T, y2d.T, sol2d_rom_u.T, sol2d_rom_v.T,
                                    color=sol2d_rom_speed.T,
                                    linewidth=linewidth,
                                    cmap=colormap,
                                    density=density)
      plot_diff = axes[2].pcolor(x2d, y2d, diff_speed, cmap=colormap)
      axes[0].set_title('FOM Streamline(u,v), t={:.3}'.format(i*dt_snapshots))
      axes[1].set_title('ROM Streamline(u,v), t={:.3}'.format(i*dt_snapshots))
      axes[2].set_title('Diff Velocity Magnitude, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_fom.lines, ax=axes[0])
      fig.colorbar(plot_rom.lines, ax=axes[1])
      fig.colorbar(plot_diff, ax=axes[2])
      if nsims > 1:
        plt_fname = plt_dir_name+'/streamlines_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/streamlines_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

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
  grid,soln_fom_snapshots = hyparutils.getSolutionSnapshots(  sim_path,
                                                              nsims,
                                                              n_snapshots,
                                                              ndims,
                                                              nvars,
                                                              size )
  grid, soln_rom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                              nsims,
                                                              n_snapshots,
                                                              ndims,
                                                              nvars,
                                                              size,
                                                              op_root='op_rom')
  soln_fom_snapshots = np.float32(soln_fom_snapshots)
  soln_rom_snapshots = np.float32(soln_rom_snapshots)
  x = grid[:size[0]]
  y = grid[size[0]:]
  print('2D Domain:');
  print(' x: ', np.min(x), np.max(x))
  print(' x.shape: ', x.shape)
  print(' y: ', np.min(y), np.max(y))
  print(' y.shape: ', y.shape)
  y2d, x2d = np.meshgrid(y, x)

  for s in range(nsims):
    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    soln_fom_snapshots_sim = soln_fom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
    sol2d_fom = np.transpose(soln_fom_snapshots_sim.reshape(size[1],size[0],nvars))
    sol2d_fom_u = sol2d_fom[1,:,:] / sol2d_fom[0,:,:]
    sol2d_fom_v = sol2d_fom[2,:,:] / sol2d_fom[0,:,:]
    sol2d_fom_speed = np.sqrt(np.square(sol2d_fom_u)+np.square(sol2d_fom_v))
    soln_rom_snapshots_sim = soln_rom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
    sol2d_rom = np.transpose(soln_rom_snapshots_sim.reshape(size[1],size[0],nvars))
    sol2d_rom_u = sol2d_rom[1,:,:] / sol2d_rom[0,:,:]
    sol2d_rom_v = sol2d_rom[2,:,:] / sol2d_rom[0,:,:]
    sol2d_rom_speed = np.sqrt(np.square(sol2d_rom_u)+np.square(sol2d_rom_v))
    diff_speed = sol2d_fom_speed - sol2d_rom_speed
    plot_fom = axes[0].streamplot(x2d.T, y2d.T, sol2d_fom_u.T, sol2d_fom_v.T,
                                  color=sol2d_fom_speed.T,
                                  linewidth=linewidth,
                                  cmap=colormap,
                                  density=density)
    plot_rom = axes[1].streamplot(x2d.T, y2d.T, sol2d_rom_u.T, sol2d_rom_v.T,
                                  color=sol2d_rom_speed.T,
                                  linewidth=linewidth,
                                  cmap=colormap,
                                  density=density)
    plot_diff = axes[2].pcolor(x2d, y2d, diff_speed, cmap=colormap)
    axes[0].set_title('FOM Streamline(u,v), t={:.3}'.format(t_final))
    axes[1].set_title('ROM Streamline(u,v), t={:.3}'.format(t_final))
    axes[2].set_title('Diff Velocity Magnitude, t={:.3}'.format(t_final))
    fig.colorbar(plot_fom.lines, ax=axes[0])
    fig.colorbar(plot_rom.lines, ax=axes[1])
    fig.colorbar(plot_diff, ax=axes[2])
    if nsims > 1:
      plt_fname = plt_dir_name+'/streamlines_'+f'{s:02d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/streamlines.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()
