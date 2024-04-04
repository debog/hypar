'''
Python script to plot the solution of flow around a cylinder.

- If op_overwrite is set to "no", a plot is generated
  for each variable  and each simulation time for which the solution
  is available.
- If op_overwrite is set to "yes", a single plot is
  created for for each variable (solution vector component).
- Solution must be in binary format

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

fig_nv=1
fig_nh=1
figsize=(1.6*fig_nh*10,fig_nv*10)
plt_dir_name='plots'

# axis to slice (0 = x, 1 = y, 2 = z)
slice_axis = 2
slice_loc = 0.5
# plot area
plt_xmin = -4
plt_xmax = 10
plt_ymin = -6
plt_ymax = 6

# fluid properties
gamma = 1.4

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
  grid, soln_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                          nsims,
                                                          n_snapshots,
                                                          ndims,
                                                          nvars,
                                                          size )
  soln_snapshots = np.float32(soln_snapshots)
  x = grid[:size[0]]
  y = grid[size[0]:(size[0]+size[1])]
  z = grid[(size[0]+size[1]):]
  print('2D Domain:');
  print(' x: ', np.min(x), np.max(x))
  print(' x.shape: ', x.shape)
  print(' y: ', np.min(y), np.max(y))
  print(' y.shape: ', y.shape)

  for i in range(n_snapshots):
    for s in range(nsims):
      soln_snapshots_sim = soln_snapshots[s*n_snapshots:(s+1)*n_snapshots]
      sol_3d = np.transpose(soln_snapshots_sim.reshape(n_snapshots,size[2],size[1],size[0],nvars))[:,:,:,:,i]
      print('3D solution shape: ', sol_3d.shape)
      if slice_axis == 0:
        print('Slicing along x at %f' % slice_loc)
        sol_2d = np.squeeze(sol_3d[:,int(slice_loc*sol_3d.shape[1]),:,:])
        y2d, x2d = np.meshgrid(z, y)
      elif slice_axis == 1:
        print('Slicing along y at %f' % slice_loc)
        sol_2d = np.squeeze(sol_3d[:,:,int(slice_loc*sol_3d.shape[2]),:])
        y2d, x2d = np.meshgrid(z, x)
      else:
        print('Slicing along z at %f' % slice_loc)
        sol_2d = np.squeeze(sol_3d[:,:,:,int(slice_loc*sol_3d.shape[3])])
        y2d, x2d = np.meshgrid(y, x)
      print('2D slice shape: ', sol_2d.shape)

      sol_2d_rho = sol_2d[0,:,:]
      sol_2d_u = sol_2d[1,:,:] / sol_2d_rho
      sol_2d_v = sol_2d[2,:,:] / sol_2d_rho
      sol_2d_w = sol_2d[3,:,:] / sol_2d_rho
      sol_2d_speed = np.sqrt(np.square(sol_2d_u)+np.square(sol_2d_v)+np.square(sol_2d_w))
      sol_2d_p = (gamma-1)*(sol_2d[4,:,:]-0.5*sol_2d_rho*(np.square(sol_2d_speed)))

      print('Solution:')
      print('  Density range: ', np.min(sol_2d_rho), np.max(sol_2d_rho))
      print('  Pressure range: ', np.min(sol_2d_p), np.max(sol_2d_p))
      print('  u range: ', np.min(sol_2d_u), np.max(sol_2d_u))
      print('  v range: ', np.min(sol_2d_v), np.max(sol_2d_v))
      print('  w range: ', np.min(sol_2d_w), np.max(sol_2d_w))

      iblank = np.ones(sol_2d_rho.shape)
      for ii in range(iblank.shape[0]):
        for jj in range(iblank.shape[1]):
          s = np.sqrt(x[ii]*x[ii]+y[jj]*y[jj])
          if s < 1.0:
            iblank[ii,jj] = np.nan

      sol_2d_rho = iblank * sol_2d_rho
      sol_2d_p = iblank * sol_2d_p
      sol_2d_u = iblank * sol_2d_u
      sol_2d_v = iblank * sol_2d_v
      sol_2d_w = iblank * sol_2d_w
      sol_2d_speed = iblank * sol_2d_speed
      sol_2d_u_wake = np.ma.masked_greater(sol_2d_u, 0)

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes.set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_rho = axes.pcolor(x2d, y2d, sol_2d_rho, cmap=colormap)
      axes.set_title('Density, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_rho, ax=axes)
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_density_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
          plt_fname = plt_dir_name+'/fig_density_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes.set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_p = axes.pcolor(x2d, y2d, sol_2d_p, cmap=colormap)
      axes.set_title('Pressure, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_p, ax=axes)
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_pressure_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
          plt_fname = plt_dir_name+'/fig_pressure_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes.set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_vel_magn = axes.pcolor(x2d, y2d, sol_2d_speed, cmap=colormap)
      axes.set_title('Velocity magnitude, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_vel_magn, ax=axes)
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_velocity_magnitude_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_velocity_magnitude_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes.set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_u = axes.pcolor(x2d, y2d, sol_2d_u, cmap=colormap)
      axes.set_title('x-Velocity, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_u, ax=axes)
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_u_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_u_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes.set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_v = axes.pcolor(x2d, y2d, sol_2d_v, cmap=colormap)
      axes.set_title('y-Velocity, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_v, ax=axes)
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_v_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_v_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes.set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_p = axes.pcolor(x2d, y2d, sol_2d_p, cmap='RdBu')
      plot_u_wake = axes.pcolor(x2d, y2d, sol_2d_u_wake, cmap=colormap)
      axes.set_title('Pressure, x-Velocity in wake (u < 0), t={:.3}'.format(i*dt_snapshots))
      cbar_p = fig.colorbar(plot_p, ax=axes)
      cbar_u_wake = fig.colorbar(plot_u_wake, ax=axes)
      cbar_p.set_label("Pressure")
      cbar_u_wake.set_label("u")
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_u_wake_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_u_wake_'+f'{i:05d}'+'.png'
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
  grid, soln_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                          nsims,
                                                          n_snapshots,
                                                          ndims,
                                                          nvars,
                                                          size )
  soln_snapshots = np.float32(soln_snapshots)
  x = grid[:size[0]]
  y = grid[size[0]:(size[0]+size[1])]
  z = grid[(size[0]+size[1]):]
  print('2D Domain:');
  print(' x: ', np.min(x), np.max(x))
  print(' x.shape: ', x.shape)
  print(' y: ', np.min(y), np.max(y))
  print(' y.shape: ', y.shape)

  for s in range(nsims):
    soln_snapshots_sim = soln_snapshots[s*n_snapshots:(s+1)*n_snapshots]
    sol_3d = np.transpose(soln_snapshots_sim.reshape(n_snapshots,size[2],size[1],size[0],nvars))[:,:,:,:]
    print('3D solution shape: ', sol_3d.shape)
    if slice_axis == 0:
      print('Slicing along x at %f' % slice_loc)
      sol_2d = np.squeeze(sol_3d[:,int(slice_loc*sol_3d.shape[1]),:,:])
      y2d, x2d = np.meshgrid(z, y)
    elif slice_axis == 1:
      print('Slicing along y at %f' % slice_loc)
      sol_2d = np.squeeze(sol_3d[:,:,int(slice_loc*sol_3d.shape[2]),:])
      y2d, x2d = np.meshgrid(z, x)
    else:
      print('Slicing along z at %f' % slice_loc)
      sol_2d = np.squeeze(sol_3d[:,:,:,int(slice_loc*sol_3d.shape[3])])
      y2d, x2d = np.meshgrid(y, x)
    print('2D slice shape: ', sol_2d.shape)

    sol_2d_rho = sol_2d[0,:,:]
    sol_2d_u = sol_2d[1,:,:] / sol_2d_rho
    sol_2d_v = sol_2d[2,:,:] / sol_2d_rho
    sol_2d_w = sol_2d[3,:,:] / sol_2d_rho
    sol_2d_speed = np.sqrt(np.square(sol_2d_u)+np.square(sol_2d_v)+np.square(sol_2d_w))
    sol_2d_p = (gamma-1)*(sol_2d[4,:,:]-0.5*sol_2d_rho*(np.square(sol_2d_speed)))

    print('Solution:')
    print('  Density range: ', np.min(sol_2d_rho), np.max(sol_2d_rho))
    print('  Pressure range: ', np.min(sol_2d_p), np.max(sol_2d_p))
    print('  u range: ', np.min(sol_2d_u), np.max(sol_2d_u))
    print('  v range: ', np.min(sol_2d_v), np.max(sol_2d_v))
    print('  w range: ', np.min(sol_2d_w), np.max(sol_2d_w))

    iblank = np.ones(sol_2d_rho.shape)
    for ii in range(iblank.shape[0]):
      for jj in range(iblank.shape[1]):
        s = np.sqrt(x[ii]*x[ii]+y[jj]*y[jj])
        if s < 1.0:
          iblank[ii,jj] = np.nan

    sol_2d_rho = iblank * sol_2d_rho
    sol_2d_p = iblank * sol_2d_p
    sol_2d_u = iblank * sol_2d_u
    sol_2d_v = iblank * sol_2d_v
    sol_2d_w = iblank * sol_2d_w
    sol_2d_speed = iblank * sol_2d_speed
    sol_2d_u_wake = np.ma.masked_greater(sol_2d_u, 0)

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes.set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_rho = axes.pcolor(x2d, y2d, sol_2d_rho, cmap=colormap)
    axes.set_title('Density, t={:.3}'.format(t_final))
    fig.colorbar(plot_rho, ax=axes)
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_density_'+f'{s:02d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_density.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes.set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_p = axes.pcolor(x2d, y2d, sol_2d_p, cmap=colormap)
    axes.set_title('Pressure, t={:.3}'.format(t_final))
    fig.colorbar(plot_p, ax=axes)
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_pressure_'+f'{s:02d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_pressure.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes.set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_vel_magn = axes.pcolor(x2d, y2d, sol_2d_speed, cmap=colormap)
    axes.set_title('Velocity magnitude, t={:.3}'.format(t_final))
    fig.colorbar(plot_vel_magn, ax=axes)
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_velocity_magnitude_'+f'{s:02d}'+'.png'
    else:
      plt_fname = plt_dir_name+'/fig_velocity_magnitude.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes.set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_u = axes.pcolor(x2d, y2d, sol_2d_u, cmap=colormap)
    axes.set_title('x-Velocity, t={:.3}'.format(t_final))
    fig.colorbar(plot_u, ax=axes)
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_u_'+f'{s:02d}'+'.png'
    else:
      plt_fname = plt_dir_name+'/fig_u.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes.set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_v = axes.pcolor(x2d, y2d, sol_2d_v, cmap=colormap)
    axes.set_title('y-Velocity, t={:.3}'.format(t_final))
    fig.colorbar(plot_v, ax=axes)
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_v_'+f'{s:02d}'+'.png'
    else:
      plt_fname = plt_dir_name+'/fig_v.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes.set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_p = axes.pcolor(x2d, y2d, sol_2d_p, cmap='RdBu')
    plot_u_wake = axes.pcolor(x2d, y2d, sol_2d_u_wake, cmap=colormap)
    axes.set_title('Pressure, x-Velocity in wake (u < 0), t={:.3}'.format(t_final))
    cbar_p = fig.colorbar(plot_p, ax=axes)
    cbar_u_wake = fig.colorbar(plot_u_wake, ax=axes)
    cbar_p.set_label("Pressure")
    cbar_u_wake.set_label("u")
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_u_wake_'+f'{s:02d}'+'.png'
    else:
      plt_fname = plt_dir_name+'/fig_u_wake.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

