'''
Python script to create plots of atmoshpheric flow variables
of a HyPar simulation of the rising thermal bubble case.
This particular file is for the
2D Navier-Stokes/Euler physics, where the solution vector
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
fig_nv=1
fig_nh=1
figsize=(2.0*fig_nh*10,fig_nv*10)

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
try:
  BV = np.float32(physics_inp_data['HB'])[1]
except:
  BV = 0.01
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
print('  B-V freq = ', BV)

def convertVar( u_conserved: np.ndarray,
                y: np.ndarray ) -> np.ndarray:
  p_exner = 1.0+((grav[1]*grav[1])/(Cp*T_ref*BV*BV))*(np.exp(-BV*BV*y/grav[1])-1.0)
  theta0 = T_ref * np.exp(BV*BV*y/grav[1])
  p0 = p_ref * np.power(p_exner, (gamma/(gamma-1.0)))
  rho0 = rho_ref * np.power(p_exner, (1.0/(gamma-1.0))) * np.exp(-BV*BV*y/grav[1])

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
                      theta0],
                    axis=0)
  return retval

def plotData( data_2d: np.ndarray,
              x2d: np.ndarray,
              y2d: np.ndarray,
              fig_filename,
              varname,
              t_solution ):

  fig = plt.figure(figsize=figsize)
  ax = fig.add_subplot(fig_nv, fig_nh, 1)
  plot_contour2= ax.pcolor( x2d, y2d, data_2d, cmap=colormap )
  cb = fig.colorbar(plot_contour2, ax=ax)
  cb.ax.set_yticklabels(["{:.4f}".format(i) for i in cb.get_ticks()])
  ax.set_xlabel("x (m)")
  ax.set_ylabel("y (m)")
  ax.set_title("{:}, t={:.1f}".format(varname, t_solution))
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
      grid,solution = hyparutils.readOpFile(fname, ndims, nvars, size)
      solution = np.float32(solution)

      x = grid[:size[0]]
      y = grid[size[0]:]
      y2d, x2d = np.meshgrid(y, x)

      u_cons_2d = np.transpose(solution.reshape(size[1],size[0],nvars))
      u_prim_2d = convertVar(u_cons_2d, y2d)

      dtheta = u_prim_2d[4,:,:]-u_prim_2d[8,:,:]
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_dtheta_'+f'{s:03d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_dtheta_'+f'{i:05d}'+'.png'
      plotData( dtheta,
                x2d, y2d,
                plt_fname,
                'dtheta',
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
    grid,solution = hyparutils.readOpFile(fname, ndims, nvars, size)
    solution = np.float32(solution)

    x = grid[:size[0]]
    y = grid[size[0]:]
    y2d, x2d = np.meshgrid(y,x)

    u_cons_2d = np.transpose(solution.reshape(size[1],size[0],nvars))
    u_prim_2d = convertVar(u_cons_2d, y2d)

    dtheta = u_prim_2d[4,:,:]-u_prim_2d[8,:,:]
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_dtheta_'+f'{s:03d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_dtheta.png'
    plotData( dtheta,
              x2d, y2d,
              plt_fname, 'dtheta',
              t_final )
