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

fig_nv=3
fig_nh=1
figsize=(1.6*fig_nh*10,fig_nv*10)
#vel_figsize=(20,10)
plt_dir_name='plots'

# streamline plot properties
#linewidth = 1
#seed_points = np.array([ [2, 2, 2, 2, 2],
#                          [-0.2, -0.1, 0.0, 0.1, 0.2] ])

# axis to slice (0 = x, 1 = y, 2 = z)
slice_axis = 2
slice_loc = 0.5
# plot area
plt_xmin = -4
plt_xmax = 10
plt_ymin = -6
plt_ymax = 6
#vel_plt_xmin = -3
#vel_plt_xmax = 7
#vel_plt_ymin = -2
#vel_plt_ymax = 2

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
  grid, soln_fom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                              nsims,
                                                              n_snapshots,
                                                              ndims,
                                                              nvars,
                                                              size )
  soln_fom_snapshots = np.float32(soln_fom_snapshots)
  grid, soln_rom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                              nsims,
                                                              n_snapshots,
                                                              ndims,
                                                              nvars,
                                                              size,
                                                              op_root='op_rom')
  soln_rom_snapshots = np.float32(soln_rom_snapshots)
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
      soln_fom_snapshots_sim = soln_fom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
      sol_fom_3d = np.transpose(soln_fom_snapshots_sim.reshape(n_snapshots,size[2],size[1],size[0],nvars))[:,:,:,:,i]
      print('3D solution shape: ', sol_fom_3d.shape)
      if slice_axis == 0:
        print('Slicing along x at %f' % slice_loc)
        sol_fom_2d = np.squeeze(sol_fom_3d[:,int(slice_loc*sol_fom_3d.shape[1]),:,:])
        y2d, x2d = np.meshgrid(z, y)
      elif slice_axis == 1:
        print('Slicing along y at %f' % slice_loc)
        sol_fom_2d = np.squeeze(sol_fom_3d[:,:,int(slice_loc*sol_fom_3d.shape[2]),:])
        y2d, x2d = np.meshgrid(z, x)
      else:
        print('Slicing along z at %f' % slice_loc)
        sol_fom_2d = np.squeeze(sol_fom_3d[:,:,:,int(slice_loc*sol_fom_3d.shape[3])])
        y2d, x2d = np.meshgrid(y, x)
      print('2D slice shape: ', sol_fom_2d.shape)

      sol_fom_2d_rho = sol_fom_2d[0,:,:]
      sol_fom_2d_u = sol_fom_2d[1,:,:] / sol_fom_2d_rho
      sol_fom_2d_v = sol_fom_2d[2,:,:] / sol_fom_2d_rho
      sol_fom_2d_w = sol_fom_2d[3,:,:] / sol_fom_2d_rho
      sol_fom_2d_speed = np.sqrt(np.square(sol_fom_2d_u)+np.square(sol_fom_2d_v)+np.square(sol_fom_2d_w))
      sol_fom_2d_p = (gamma-1)*(sol_fom_2d[4,:,:]-0.5*sol_fom_2d_rho*(np.square(sol_fom_2d_speed)))

      print('FOM Solution:')
      print('  Density range: ', np.min(sol_fom_2d_rho), np.max(sol_fom_2d_rho))
      print('  Pressure range: ', np.min(sol_fom_2d_p), np.max(sol_fom_2d_p))
      print('  u range: ', np.min(sol_fom_2d_u), np.max(sol_fom_2d_u))
      print('  v range: ', np.min(sol_fom_2d_v), np.max(sol_fom_2d_v))
      print('  w range: ', np.min(sol_fom_2d_w), np.max(sol_fom_2d_w))

      soln_rom_snapshots_sim = soln_rom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
      sol_rom_3d = np.transpose(soln_rom_snapshots_sim.reshape(n_snapshots,size[2],size[1],size[0],nvars))[:,:,:,:,i]
      print('3D solution shape: ', sol_rom_3d.shape)
      if slice_axis == 0:
        print('Slicing along x at %f' % slice_loc)
        sol_rom_2d = np.squeeze(sol_rom_3d[:,int(slice_loc*sol_rom_3d.shape[1]),:,:])
        y2d, x2d = np.meshgrid(z, y)
      elif slice_axis == 1:
        print('Slicing along y at %f' % slice_loc)
        sol_rom_2d = np.squeeze(sol_rom_3d[:,:,int(slice_loc*sol_rom_3d.shape[2]),:])
        y2d, x2d = np.meshgrid(z, x)
      else:
        print('Slicing along z at %f' % slice_loc)
        sol_rom_2d = np.squeeze(sol_rom_3d[:,:,:,int(slice_loc*sol_rom_3d.shape[3])])
        y2d, x2d = np.meshgrid(y, x)
      print('2D slice shape: ', sol_rom_2d.shape)

      sol_rom_2d_rho = sol_rom_2d[0,:,:]
      sol_rom_2d_u = sol_rom_2d[1,:,:] / sol_rom_2d_rho
      sol_rom_2d_v = sol_rom_2d[2,:,:] / sol_rom_2d_rho
      sol_rom_2d_w = sol_rom_2d[3,:,:] / sol_rom_2d_rho
      sol_rom_2d_speed = np.sqrt(np.square(sol_rom_2d_u)+np.square(sol_rom_2d_v)+np.square(sol_rom_2d_w))
      sol_rom_2d_p = (gamma-1)*(sol_rom_2d[4,:,:]-0.5*sol_rom_2d_rho*(np.square(sol_rom_2d_speed)))

      print('ROM Solution:')
      print('  Density range: ', np.min(sol_rom_2d_rho), np.max(sol_rom_2d_rho))
      print('  Pressure range: ', np.min(sol_rom_2d_p), np.max(sol_rom_2d_p))
      print('  u range: ', np.min(sol_rom_2d_u), np.max(sol_rom_2d_u))
      print('  v range: ', np.min(sol_rom_2d_v), np.max(sol_rom_2d_v))
      print('  w range: ', np.min(sol_rom_2d_w), np.max(sol_rom_2d_w))

      iblank = np.ones(sol_rom_2d_rho.shape)
      for ii in range(iblank.shape[0]):
        for jj in range(iblank.shape[1]):
          s = np.sqrt(x[ii]*x[ii]+y[jj]*y[jj])
          if s < 1.0:
            iblank[ii,jj] = np.nan

      sol_fom_2d_rho = iblank * sol_fom_2d_rho
      sol_fom_2d_p = iblank * sol_fom_2d_p
      sol_fom_2d_u = iblank * sol_fom_2d_u
      sol_fom_2d_v = iblank * sol_fom_2d_v
      sol_fom_2d_w = iblank * sol_fom_2d_w
      sol_fom_2d_speed = iblank * sol_fom_2d_speed
      sol_fom_2d_u_wake = np.ma.masked_greater(sol_fom_2d_u, 0)

      sol_rom_2d_rho = iblank * sol_rom_2d_rho
      sol_rom_2d_p = iblank * sol_rom_2d_p
      sol_rom_2d_u = iblank * sol_rom_2d_u
      sol_rom_2d_v = iblank * sol_rom_2d_v
      sol_rom_2d_w = iblank * sol_rom_2d_w
      sol_rom_2d_speed = iblank * sol_rom_2d_speed
      sol_rom_2d_u_wake = np.ma.masked_greater(sol_rom_2d_u, 0)

      sol_diff_2d_rho = sol_fom_2d_rho - sol_rom_2d_rho
      sol_diff_2d_p = sol_fom_2d_p - sol_rom_2d_p
      sol_diff_2d_u = sol_fom_2d_u - sol_rom_2d_u
      sol_diff_2d_v = sol_fom_2d_v - sol_rom_2d_v
      sol_diff_2d_w = sol_fom_2d_w - sol_rom_2d_w
      sol_diff_2d_speed = sol_fom_2d_speed - sol_rom_2d_speed
      sol_diff_2d_u_wake = sol_fom_2d_u_wake - sol_rom_2d_u_wake

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes[n].set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_rho_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_rho, cmap=colormap)
      axes[0].set_title('FOM Density, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_rho_fom, ax=axes[0])
      plot_rho_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_rho, cmap=colormap)
      axes[1].set_title('ROM Density, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_rho_rom, ax=axes[1])
      plot_rho_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_rho, cmap=colormap)
      axes[2].set_title('Diff Density, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_rho_diff, ax=axes[2])
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_density_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
          plt_fname = plt_dir_name+'/fig_density_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes[n].set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_p_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_p, cmap=colormap)
      axes[0].set_title('FOM Pressure, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_p_fom, ax=axes[0])
      plot_p_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_p, cmap=colormap)
      axes[1].set_title('ROM Pressure, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_p_rom, ax=axes[1])
      plot_p_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_p, cmap=colormap)
      axes[2].set_title('Diff Pressure, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_p_diff, ax=axes[2])
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_pressure_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
          plt_fname = plt_dir_name+'/fig_pressure_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes[n].set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_vel_magn_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_speed, cmap=colormap)
      axes[0].set_title('FOM Velocity magnitude, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_vel_magn_fom, ax=axes[0])
      plot_vel_magn_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_speed, cmap=colormap)
      axes[1].set_title('ROM Velocity magnitude, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_vel_magn_rom, ax=axes[1])
      plot_vel_magn_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_speed, cmap=colormap)
      axes[2].set_title('Diff Velocity magnitude, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_vel_magn_diff, ax=axes[2])
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_velocity_magnitude_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_velocity_magnitude_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes[n].set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_u_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_u, cmap=colormap)
      axes[0].set_title('FOM x-Velocity, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_u_fom, ax=axes[0])
      plot_u_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_u, cmap=colormap)
      axes[1].set_title('ROM x-Velocity, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_u_rom, ax=axes[1])
      plot_u_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_u, cmap=colormap)
      axes[2].set_title('Diff x-Velocity, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_u_diff, ax=axes[2])
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_u_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_u_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes[n].set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_v_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_v, cmap=colormap)
      axes[0].set_title('FOM y-Velocity, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_v_fom, ax=axes[0])
      plot_v_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_v, cmap=colormap)
      axes[1].set_title('ROM y-Velocity, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_v_rom, ax=axes[1])
      plot_v_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_v, cmap=colormap)
      axes[2].set_title('Diff x-Velocity, t={:.3}'.format(i*dt_snapshots))
      fig.colorbar(plot_v_diff, ax=axes[2])
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_v_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_v_'+f'{i:05d}'+'.png'
      print('Saving %s' % plt_fname)
      plt.savefig(plt_fname)
      plt.close()

      fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
      for n in range(fig_nv*fig_nh):
        axes[n].set(  xlim=(plt_xmin, plt_xmax),
                      ylim=(plt_ymin, plt_ymax) )
      plot_p_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_p, cmap='RdBu')
      plot_u_wake_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_u_wake, cmap=colormap)
      axes[0].set_title('FOM Pressure, x-Velocity in wake (u < 0), t={:.3}'.format(i*dt_snapshots))
      cbar_p_fom = fig.colorbar(plot_p_fom, ax=axes[0])
      cbar_u_wake_fom = fig.colorbar(plot_u_wake_fom, ax=axes[0])
      cbar_p_fom.set_label("Pressure")
      cbar_u_wake_fom.set_label("u")
      plot_p_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_p, cmap='RdBu')
      plot_u_wake_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_u_wake, cmap=colormap)
      axes[1].set_title('ROM Pressure, x-Velocity in wake (u < 0), t={:.3}'.format(i*dt_snapshots))
      cbar_p_rom = fig.colorbar(plot_p_rom, ax=axes[1])
      cbar_u_wake_rom = fig.colorbar(plot_u_wake_rom, ax=axes[1])
      cbar_p_rom.set_label("Pressure")
      cbar_u_wake_rom.set_label("u")
      plot_p_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_p, cmap='RdBu')
      plot_u_wake_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_u_wake, cmap=colormap)
      axes[2].set_title('Diff Pressure, x-Velocity in wake (u < 0), t={:.3}'.format(i*dt_snapshots))
      cbar_p_diff = fig.colorbar(plot_p_diff, ax=axes[2])
      cbar_u_wake_diff = fig.colorbar(plot_u_wake_diff, ax=axes[2])
      cbar_p_diff.set_label("Pressure")
      cbar_u_wake_diff.set_label("u")
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
  grid, soln_fom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                              nsims,
                                                              n_snapshots,
                                                              ndims,
                                                              nvars,
                                                              size )
  soln_fom_snapshots = np.float32(soln_fom_snapshots)
  grid, soln_rom_snapshots = hyparutils.getSolutionSnapshots( sim_path,
                                                              nsims,
                                                              n_snapshots,
                                                              ndims,
                                                              nvars,
                                                              size,
                                                              op_root='op_rom')
  soln_rom_snapshots = np.float32(soln_rom_snapshots)
  x = grid[:size[0]]
  y = grid[size[0]:(size[0]+size[1])]
  z = grid[(size[0]+size[1]):]
  print('2D Domain:');
  print(' x: ', np.min(x), np.max(x))
  print(' x.shape: ', x.shape)
  print(' y: ', np.min(y), np.max(y))
  print(' y.shape: ', y.shape)

  for s in range(nsims):
    soln_fom_snapshots_sim = soln_fom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
    sol_fom_3d = np.transpose(soln_fom_snapshots_sim.reshape(n_snapshots,size[2],size[1],size[0],nvars))[:,:,:,:]
    print('3D solution shape: ', sol_fom_3d.shape)
    if slice_axis == 0:
      print('Slicing along x at %f' % slice_loc)
      sol_fom_2d = np.squeeze(sol_fom_3d[:,int(slice_loc*sol_fom_3d.shape[1]),:,:])
      y2d, x2d = np.meshgrid(z, y)
    elif slice_axis == 1:
      print('Slicing along y at %f' % slice_loc)
      sol_fom_2d = np.squeeze(sol_fom_3d[:,:,int(slice_loc*sol_fom_3d.shape[2]),:])
      y2d, x2d = np.meshgrid(z, x)
    else:
      print('Slicing along z at %f' % slice_loc)
      sol_fom_2d = np.squeeze(sol_fom_3d[:,:,:,int(slice_loc*sol_fom_3d.shape[3])])
      y2d, x2d = np.meshgrid(y, x)
    print('2D slice shape: ', sol_fom_2d.shape)

    sol_fom_2d_rho = sol_fom_2d[0,:,:]
    sol_fom_2d_u = sol_fom_2d[1,:,:] / sol_fom_2d_rho
    sol_fom_2d_v = sol_fom_2d[2,:,:] / sol_fom_2d_rho
    sol_fom_2d_w = sol_fom_2d[3,:,:] / sol_fom_2d_rho
    sol_fom_2d_speed = np.sqrt(np.square(sol_fom_2d_u)+np.square(sol_fom_2d_v)+np.square(sol_fom_2d_w))
    sol_fom_2d_p = (gamma-1)*(sol_fom_2d[4,:,:]-0.5*sol_fom_2d_rho*(np.square(sol_fom_2d_speed)))

    print('FOM Solution:')
    print('  Density range: ', np.min(sol_fom_2d_rho), np.max(sol_fom_2d_rho))
    print('  Pressure range: ', np.min(sol_fom_2d_p), np.max(sol_fom_2d_p))
    print('  u range: ', np.min(sol_fom_2d_u), np.max(sol_fom_2d_u))
    print('  v range: ', np.min(sol_fom_2d_v), np.max(sol_fom_2d_v))
    print('  w range: ', np.min(sol_fom_2d_w), np.max(sol_fom_2d_w))

    soln_rom_snapshots_sim = soln_rom_snapshots[s*n_snapshots:(s+1)*n_snapshots]
    sol_rom_3d = np.transpose(soln_rom_snapshots_sim.reshape(n_snapshots,size[2],size[1],size[0],nvars))[:,:,:,:]
    print('3D solution shape: ', sol_rom_3d.shape)
    if slice_axis == 0:
      print('Slicing along x at %f' % slice_loc)
      sol_rom_2d = np.squeeze(sol_rom_3d[:,int(slice_loc*sol_rom_3d.shape[1]),:,:])
      y2d, x2d = np.meshgrid(z, y)
    elif slice_axis == 1:
      print('Slicing along y at %f' % slice_loc)
      sol_rom_2d = np.squeeze(sol_rom_3d[:,:,int(slice_loc*sol_rom_3d.shape[2]),:])
      y2d, x2d = np.meshgrid(z, x)
    else:
      print('Slicing along z at %f' % slice_loc)
      sol_rom_2d = np.squeeze(sol_rom_3d[:,:,:,int(slice_loc*sol_rom_3d.shape[3])])
      y2d, x2d = np.meshgrid(y, x)
    print('2D slice shape: ', sol_rom_2d.shape)

    sol_rom_2d_rho = sol_rom_2d[0,:,:]
    sol_rom_2d_u = sol_rom_2d[1,:,:] / sol_rom_2d_rho
    sol_rom_2d_v = sol_rom_2d[2,:,:] / sol_rom_2d_rho
    sol_rom_2d_w = sol_rom_2d[3,:,:] / sol_rom_2d_rho
    sol_rom_2d_speed = np.sqrt(np.square(sol_rom_2d_u)+np.square(sol_rom_2d_v)+np.square(sol_rom_2d_w))
    sol_rom_2d_p = (gamma-1)*(sol_rom_2d[4,:,:]-0.5*sol_rom_2d_rho*(np.square(sol_rom_2d_speed)))

    print('ROM Solution:')
    print('  Density range: ', np.min(sol_rom_2d_rho), np.max(sol_rom_2d_rho))
    print('  Pressure range: ', np.min(sol_rom_2d_p), np.max(sol_rom_2d_p))
    print('  u range: ', np.min(sol_rom_2d_u), np.max(sol_rom_2d_u))
    print('  v range: ', np.min(sol_rom_2d_v), np.max(sol_rom_2d_v))
    print('  w range: ', np.min(sol_rom_2d_w), np.max(sol_rom_2d_w))

    iblank = np.ones(sol_rom_2d_rho.shape)
    for ii in range(iblank.shape[0]):
      for jj in range(iblank.shape[1]):
        s = np.sqrt(x[ii]*x[ii]+y[jj]*y[jj])
        if s < 1.0:
          iblank[ii,jj] = np.nan

    sol_fom_2d_rho = iblank * sol_fom_2d_rho
    sol_fom_2d_p = iblank * sol_fom_2d_p
    sol_fom_2d_u = iblank * sol_fom_2d_u
    sol_fom_2d_v = iblank * sol_fom_2d_v
    sol_fom_2d_w = iblank * sol_fom_2d_w
    sol_fom_2d_speed = iblank * sol_fom_2d_speed
    sol_fom_2d_u_wake = np.ma.masked_greater(sol_fom_2d_u, 0)

    sol_rom_2d_rho = iblank * sol_rom_2d_rho
    sol_rom_2d_p = iblank * sol_rom_2d_p
    sol_rom_2d_u = iblank * sol_rom_2d_u
    sol_rom_2d_v = iblank * sol_rom_2d_v
    sol_rom_2d_w = iblank * sol_rom_2d_w
    sol_rom_2d_speed = iblank * sol_rom_2d_speed
    sol_rom_2d_u_wake = np.ma.masked_greater(sol_rom_2d_u, 0)

    sol_diff_2d_rho = sol_fom_2d_rho - sol_rom_2d_rho
    sol_diff_2d_p = sol_fom_2d_p - sol_rom_2d_p
    sol_diff_2d_u = sol_fom_2d_u - sol_rom_2d_u
    sol_diff_2d_v = sol_fom_2d_v - sol_rom_2d_v
    sol_diff_2d_w = sol_fom_2d_w - sol_rom_2d_w
    sol_diff_2d_speed = sol_fom_2d_speed - sol_rom_2d_speed
    sol_diff_2d_u_wake = sol_fom_2d_u_wake - sol_rom_2d_u_wake

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes[n].set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_rho_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_rho, cmap=colormap)
    axes[0].set_title('FOM Density, t={:.3}'.format(t_final))
    fig.colorbar(plot_rho_fom, ax=axes[0])
    plot_rho_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_rho, cmap=colormap)
    axes[1].set_title('ROM Density, t={:.3}'.format(t_final))
    fig.colorbar(plot_rho_rom, ax=axes[1])
    plot_rho_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_rho, cmap=colormap)
    axes[2].set_title('Diff Density, t={:.3}'.format(t_final))
    fig.colorbar(plot_rho_diff, ax=axes[2])
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_density_'+f'{s:02d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_density.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes[n].set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_p_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_p, cmap=colormap)
    axes[0].set_title('FOM Pressure, t={:.3}'.format(t_final))
    fig.colorbar(plot_p_fom, ax=axes[0])
    plot_p_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_p, cmap=colormap)
    axes[1].set_title('ROM Pressure, t={:.3}'.format(t_final))
    fig.colorbar(plot_p_rom, ax=axes[1])
    plot_p_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_p, cmap=colormap)
    axes[2].set_title('Diff Pressure, t={:.3}'.format(t_final))
    fig.colorbar(plot_p_diff, ax=axes[2])
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_pressure_'+f'{s:02d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_pressure.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes[n].set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_vel_magn_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_speed, cmap=colormap)
    axes[0].set_title('FOM Velocity magnitude, t={:.3}'.format(t_final))
    fig.colorbar(plot_vel_magn_fom, ax=axes[0])
    plot_vel_magn_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_speed, cmap=colormap)
    axes[1].set_title('ROM Velocity magnitude, t={:.3}'.format(t_final))
    fig.colorbar(plot_vel_magn_rom, ax=axes[1])
    plot_vel_magn_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_speed, cmap=colormap)
    axes[2].set_title('Diff Velocity magnitude, t={:.3}'.format(t_final))
    fig.colorbar(plot_vel_magn_diff, ax=axes[2])
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_velocity_magnitude_'+f'{s:02d}'+'.png'
    else:
      plt_fname = plt_dir_name+'/fig_velocity_magnitude.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes[n].set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_u_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_u, cmap=colormap)
    axes[0].set_title('FOM x-Velocity, t={:.3}'.format(t_final))
    fig.colorbar(plot_u_fom, ax=axes[0])
    plot_u_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_u, cmap=colormap)
    axes[1].set_title('ROM x-Velocity, t={:.3}'.format(t_final))
    fig.colorbar(plot_u_rom, ax=axes[1])
    plot_u_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_u, cmap=colormap)
    axes[2].set_title('Diff x-Velocity, t={:.3}'.format(t_final))
    fig.colorbar(plot_u_diff, ax=axes[2])
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_u_'+f'{s:02d}'+'.png'
    else:
      plt_fname = plt_dir_name+'/fig_u.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes[n].set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_v_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_v, cmap=colormap)
    axes[0].set_title('FOM y-Velocity, t={:.3}'.format(t_final))
    fig.colorbar(plot_v_fom, ax=axes[0])
    plot_v_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_v, cmap=colormap)
    axes[1].set_title('ROM y-Velocity, t={:.3}'.format(t_final))
    fig.colorbar(plot_v_rom, ax=axes[1])
    plot_v_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_v, cmap=colormap)
    axes[2].set_title('Diff x-Velocity, t={:.3}'.format(t_final))
    fig.colorbar(plot_v_diff, ax=axes[2])
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_v_'+f'{s:02d}'+'.png'
    else:
      plt_fname = plt_dir_name+'/fig_v.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

    fig, axes = plt.subplots(fig_nv,fig_nh,figsize=figsize)
    for n in range(fig_nv*fig_nh):
      axes[n].set(  xlim=(plt_xmin, plt_xmax),
                    ylim=(plt_ymin, plt_ymax) )
    plot_p_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_p, cmap='RdBu')
    plot_u_wake_fom = axes[0].pcolor(x2d, y2d, sol_fom_2d_u_wake, cmap=colormap)
    axes[0].set_title('FOM Pressure, x-Velocity in wake (u < 0), t={:.3}'.format(t_final))
    cbar_p_fom = fig.colorbar(plot_p_fom, ax=axes[0])
    cbar_u_wake_fom = fig.colorbar(plot_u_wake_fom, ax=axes[0])
    cbar_p_fom.set_label("Pressure")
    cbar_u_wake_fom.set_label("u")
    plot_p_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_p, cmap='RdBu')
    plot_u_wake_rom = axes[1].pcolor(x2d, y2d, sol_rom_2d_u_wake, cmap=colormap)
    axes[1].set_title('ROM Pressure, x-Velocity in wake (u < 0), t={:.3}'.format(t_final))
    cbar_p_rom = fig.colorbar(plot_p_rom, ax=axes[1])
    cbar_u_wake_rom = fig.colorbar(plot_u_wake_rom, ax=axes[1])
    cbar_p_rom.set_label("Pressure")
    cbar_u_wake_rom.set_label("u")
    plot_p_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_p, cmap='RdBu')
    plot_u_wake_diff = axes[2].pcolor(x2d, y2d, sol_diff_2d_u_wake, cmap=colormap)
    axes[2].set_title('Diff Pressure, x-Velocity in wake (u < 0), t={:.3}'.format(t_final))
    cbar_p_diff = fig.colorbar(plot_p_diff, ax=axes[2])
    cbar_u_wake_diff = fig.colorbar(plot_u_wake_diff, ax=axes[2])
    cbar_p_diff.set_label("Pressure")
    cbar_u_wake_diff.set_label("u")
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_u_wake_'+f'{s:02d}'+'.png'
    else:
      plt_fname = plt_dir_name+'/fig_u_wake.png'
    print('Saving %s' % plt_fname)
    plt.savefig(plt_fname)
    plt.close()

