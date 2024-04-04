'''
Python script to create plots of atmoshpheric flow variables
of a HyPar simulation of the rising thermal bubble case.
This particular file is for the
3D Navier-Stokes/Euler physics, where the solution vector
components are (rho, rho*u, rho*v, e).

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
import numpy.linalg as LA
import matplotlib
import matplotlib.pyplot as plt

hypar_dir = os.environ.get('HYPAR_DIR')
hypar_dir_python = hypar_dir + '/Examples/Python'
sys.path.append(hypar_dir_python)

import modHyParUtils as hyparutils

font = {'size':22}
matplotlib.rc('font', **font)
colormap='jet'

plt_dir_name='plots'
fig_nv=3
fig_nh=1
figsize=(1.4*fig_nh*10,fig_nv*10)

'''
Set up the simulation parameters
'''
sim_path = '.'
sim_inp_data = hyparutils.readHyParInpFile(sim_path+'/simulation.inp')
solver_inp_data = hyparutils.readHyParInpFile(sim_path+'/solver.inp')
physics_inp_data = hyparutils.readHyParInpFile(sim_path+'/physics.inp')

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

R = 287.058
T_ref = 300.0
try:
  grav = np.float32(physics_inp_data['gravity'])
except:
  grav = np.float32(np.array([0.0, 9.8, 0.0]))
try:
  gamma = np.float32(physics_inp_data['gamma'])[0]
except:
  gamma = 1.4
try:
  p_ref = np.float32(physics_inp_data['p_ref'])[0]
except:
  p_ref = 100000.0
rho_ref = p_ref / (R*T_ref)
Cp = gamma / (gamma-1.0) * R

print('Atmospheric Flow Paramters (SI Units):')
print('  R        = ', R)
print('  T_ref    = ', T_ref)
print('  p_ref    = ', p_ref)
print('  rho_ref  = ', rho_ref)
print('  gamma    = ', gamma)
print('  gravity  = ', grav)
print('  Cp       = ', Cp)

def convertVar( u_conserved: np.ndarray,
                y: np.ndarray ) -> np.ndarray:
  theta0 = T_ref
  p_exner = 1.0 - (grav[1]*y)/(Cp*T_ref)
  rho0 = (p_ref/(R*theta0)) * np.power(p_exner, (1.0/(gamma-1.0)))
  p0 = p_ref * np.power(p_exner, (gamma/(gamma-1.0)))

  rho = u_conserved[0,:,:]
  uvel = u_conserved[1,:,:] / rho
  vvel = u_conserved[2,:,:] / rho
  ene = u_conserved[3,:,:]
  ke = 0.5 * (np.square(uvel) + np.square(vvel))
  pressure = (gamma-1.0) * (ene - rho*ke)
  theta = ((gamma-1.0)/R) * (ene - rho*ke) / (rho*p_exner)

  retval = np.stack([ rho,
                      uvel,
                      vvel,
                      pressure,
                      theta,
                      rho0,
                      p0,
                      p_exner,
                      theta0*np.ones(rho.shape)],
                    axis=0)
  return retval

def plotData( data_fom_2d: np.ndarray,
              data_rom_2d: np.ndarray,
              x2d: np.ndarray,
              y2d: np.ndarray,
              fig_filename,
              varname,
              t_solution ):

  data_del_2d = data_fom_2d - data_rom_2d

  fom_sol_norm = LA.norm(data_fom_2d.ravel())
  diff_norm = LA.norm(data_del_2d.ravel())
  if fom_sol_norm > 1e-14:
    diff_norm = diff_norm / fom_sol_norm

  fig = plt.figure(figsize=figsize)
  ax = fig.add_subplot(fig_nv, fig_nh, 1)
  plot_contour2= ax.pcolor( x2d, y2d, data_fom_2d, cmap=colormap )
  cb = fig.colorbar(plot_contour2, ax=ax)
  cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
  ax.set_xlabel("x (m)")
  ax.set_ylabel("y (m)")
  ax.set_title("FOM {:}, t={:.1f}".format(varname, t_solution))
  ax.set_xlim(np.min(x2d), np.max(x2d))
  ax.set_ylim(np.min(y2d), np.max(y2d))

  ax = fig.add_subplot(fig_nv, fig_nh, 2)
  plot_contour2= ax.pcolor( x2d, y2d, data_rom_2d, cmap=colormap )
  cb = fig.colorbar(plot_contour2, ax=ax)
  cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
  ax.set_xlabel("x (m)")
  ax.set_ylabel("y (m)")
  ax.set_title("ROM {:}, t={:.1f}".format(varname, t_solution))
  ax.set_xlim(np.min(x2d), np.max(x2d))
  ax.set_ylim(np.min(y2d), np.max(y2d))

  ax = fig.add_subplot(fig_nv, fig_nh, 3)
  plot_contour2= ax.pcolor( x2d, y2d, data_del_2d, cmap=colormap )
  cb = fig.colorbar(plot_contour2, ax=ax)
  cb.ax.set_yticklabels(["{:.1e}".format(i) for i in cb.get_ticks()])
  ax.set_xlabel("x (m)")
  ax.set_ylabel("z (m)")
  ax.set_title("Diff {:}, t={:.1f}, norm={:.3e}".format(varname,
                                                      t_solution,
                                                      diff_norm))
  ax.set_xlim(np.min(x2d), np.max(x2d))
  ax.set_ylim(np.min(y2d), np.max(y2d))

  print('Saving %s' % fig_filename)
  plt.savefig(fig_filename)
  plt.close()

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

  for i in range(n_snapshots):
    for sim in range(nsims):

      '''
      Load simulation data (solution snapshots)
      '''
      op_root='op'
      if nsims >= 100:
        fname = sim_path + '/'+op_root+'_'+f'{sim:02d}'+'_'+f'{i:05d}'+'.bin'
      elif nsims >= 10:
        fname = sim_path + '/'+op_root+'_'+f'{sim:02d}'+'_'+f'{i:05d}'+'.bin'
      elif nsims > 1:
        fname = sim_path + '/'+op_root+'_'+f'{sim:01d}'+'_'+f'{i:05d}'+'.bin'
      else:
        fname = sim_path + '/'+op_root+'_'+f'{i:05d}'+'.bin'
      grid,solution_fom = hyparutils.readOpFile(fname, ndims, nvars, size)
      solution_fom = np.float32(solution_fom)
      op_root='op_rom'
      if nsims >= 100:
        fname = sim_path + '/'+op_root+'_'+f'{sim:02d}'+'_'+f'{i:05d}'+'.bin'
      elif nsims >= 10:
        fname = sim_path + '/'+op_root+'_'+f'{sim:02d}'+'_'+f'{i:05d}'+'.bin'
      elif nsims > 1:
        fname = sim_path + '/'+op_root+'_'+f'{sim:01d}'+'_'+f'{i:05d}'+'.bin'
      else:
        fname = sim_path + '/'+op_root+'_'+f'{i:05d}'+'.bin'
      grid,solution_rom = hyparutils.readOpFile(fname, ndims, nvars, size)
      solution_rom = np.float32(solution_rom)

      x = grid[:size[0]]
      y = grid[size[0]:]
      y2d, x2d = np.meshgrid(y, x)

      u_fom_cons_2d = np.transpose(solution_fom.reshape(size[1],size[0],nvars))
      u_fom_prim_2d = convertVar(u_fom_cons_2d, y2d)

      u_rom_cons_2d = np.transpose(solution_rom.reshape(size[1],size[0],nvars))
      u_rom_prim_2d = convertVar(u_rom_cons_2d, y2d)

      theta_fom = u_fom_prim_2d[4,:,:]
      theta_rom = u_rom_prim_2d[4,:,:]
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_theta_'+f'{s:03d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_theta_'+f'{i:05d}'+'.png'
      plotData( theta_fom,
                theta_rom,
                x2d, y2d,
                plt_fname,
                'theta',
                i*dt_snapshots)

      uvel_fom = u_fom_prim_2d[1,:,:]
      uvel_rom = u_rom_prim_2d[1,:,:]
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_uvel_'+f'{s:03d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_uvel_'+f'{i:05d}'+'.png'
      plotData( uvel_fom,
                uvel_rom,
                x2d, y2d,
                plt_fname,
                'uvel',
                i*dt_snapshots)

      vvel_fom = u_fom_prim_2d[2,:,:]
      vvel_rom = u_rom_prim_2d[2,:,:]
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_vvel_'+f'{s:03d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_vvel_'+f'{i:05d}'+'.png'
      plotData( vvel_fom,
                vvel_rom,
                x2d, y2d,
                plt_fname,
                'vvel',
                i*dt_snapshots)

else:

  niter = int(solver_inp_data['n_iter'][0])
  dt = float(solver_inp_data['dt'][0])
  t_final = dt*niter

  print('Simulation parameters:')
  print('  ndims = ', ndims)
  print('  nvars = ', nvars)
  print('  grid size = ', size)
  print('  final time = ', t_final)

  for sim in range(nsims):

    '''
    Load simulation data (solution snapshots)
    '''
    op_root='op'
    if nsims >= 100:
      fname = sim_path + '/'+op_root+'_'+f'{sim:03d}'+'.bin'
    elif nsims >= 10:
      fname = sim_path + '/'+op_root+'_'+f'{sim:02d}'+'.bin'
    elif nsims > 1:
      fname = sim_path + '/'+op_root+'_'+f'{sim:01d}'+'.bin'
    else:
      fname = sim_path + '/'+op_root+'.bin'
    grid,solution_fom = hyparutils.readOpFile(fname, ndims, nvars, size)
    solution_fom = np.float32(solution_fom)
    op_root='op_rom'
    if nsims >= 100:
      fname = sim_path + '/'+op_root+'_'+f'{sim:03d}'+'.bin'
    elif nsims >= 10:
      fname = sim_path + '/'+op_root+'_'+f'{sim:02d}'+'.bin'
    elif nsims > 1:
      fname = sim_path + '/'+op_root+'_'+f'{sim:01d}'+'.bin'
    else:
      fname = sim_path + '/'+op_root+'.bin'
    grid,solution_rom = hyparutils.readOpFile(fname, ndims, nvars, size)
    solution_rom = np.float32(solution_rom)

    x = grid[:size[0]]
    y = grid[size[0]:]
    y2d, x2d = np.meshgrid(y,x)

    u_fom_cons_2d = np.transpose(solution_fom.reshape(size[1],size[0],nvars))
    u_fom_prim_2d = convertVar(u_fom_cons_2d, y2d)

    u_rom_cons_2d = np.transpose(solution_rom.reshape(size[1],size[0],nvars))
    u_rom_prim_2d = convertVar(u_rom_cons_2d, y2d)

    theta_fom = u_fom_prim_2d[4,:,:]
    theta_rom = u_rom_prim_2d[4,:,:]
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_theta_'+f'{s:03d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_theta.png'
    plotData( theta_fom,
              theta_rom,
              x2d, y2d,
              plt_fname, 'theta',
              t_final )

    uvel_fom = u_fom_prim_2d[1,:,:]
    uvel_rom = u_rom_prim_2d[1,:,:]
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_uvel_'+f'{s:03d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_uvel.png'
    plotData( uvel_fom,
              uvel_rom,
              x2d, y2d,
              plt_fname, 'uvel',
              t_final )

    vvel_fom = u_fom_prim_2d[2,:,:]
    vvel_rom = u_rom_prim_2d[2,:,:]
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_vvel_'+f'{s:03d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_vvel.png'
    plotData( vvel_fom,
              vvel_rom,
              x2d, y2d,
              plt_fname, 'vvel',
              t_final )

