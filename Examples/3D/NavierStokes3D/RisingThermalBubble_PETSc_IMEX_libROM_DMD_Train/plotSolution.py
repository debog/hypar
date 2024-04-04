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

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

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

density = 4
linewidth = 1
slice_loc = 0.5

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

  rho = u_conserved[0,:,:,:]
  uvel = u_conserved[1,:,:,:] / rho
  vvel = u_conserved[2,:,:,:] / rho
  wvel = u_conserved[3,:,:,:] / rho
  ene = u_conserved[4,:,:,:]
  ke = 0.5 * (np.square(uvel) + np.square(vvel) + np.square(wvel))
  pressure = (gamma-1.0) * (ene - rho*ke)
  theta = ((gamma-1.0)/R) * (ene - rho*ke) / (rho*p_exner)

  retval = np.stack([ rho,
                      uvel,
                      vvel,
                      wvel,
                      pressure,
                      theta,
                      rho0,
                      p0,
                      p_exner,
                      theta0*np.ones(rho.shape)],
                    axis=0)
  return retval

def plotData( data_fom_3d: np.ndarray,
              data_rom_3d: np.ndarray,
              x3d: np.ndarray,
              y3d: np.ndarray,
              z3d: np.ndarray,
              fig_filename,
              varname,
              t_solution ):

  # For ease of visualization, switch y and z. In the simulation,
  # y is the vertical axis, but that makes things difficult for a
  # 3D plot!
  data_fom_3d = np.transpose(data_fom_3d,[0,2,1])
  data_rom_3d = np.transpose(data_rom_3d,[0,2,1])
  data_del_3d = data_fom_3d - data_rom_3d

  fom_sol_norm = LA.norm(data_fom_3d.ravel())
  diff_norm = LA.norm(data_del_3d.ravel())
  if fom_sol_norm > 1e-14:
    diff_norm = diff_norm / fom_sol_norm

  fig = plt.figure(figsize=figsize)

  kw1 =  {'vmin': data_fom_3d.min(),
          'vmax': data_fom_3d.max(),
          'levels': np.linspace((data_fom_3d.min()+0.1*(data_fom_3d.max()-data_fom_3d.min())),
                                (data_fom_3d.max()-0.1*(data_fom_3d.max()-data_fom_3d.min())),
                                25),
          'antialiased' : True,
          'alpha' : 0.8,
          'cmap' : colormap,
        }
  kw2 =  {'vmin': data_fom_3d.min(),
          'vmax': data_fom_3d.max(),
          'levels': np.linspace(data_fom_3d.min(),
                                data_fom_3d.max(),
                                100),
          'antialiased' : True,
          'alpha' : 0.6,
          'cmap' : colormap,
        }
  ax = fig.add_subplot(fig_nv, fig_nh, 1, projection='3d')
  plot_contour1= ax.contour(  x3d[:, 0, :],
                              data_fom_3d[:,int(slice_loc*data_fom_3d.shape[1]),:],
                              z3d[:, 0, :],
                              zdir='y', offset=y.max()/2, **kw1 )
  plot_contour2= ax.contourf( data_fom_3d[int(slice_loc*data_fom_3d.shape[0]), :, :],
                              y3d[0, :, :],
                              z3d[0, :, :],
                              zdir='x', offset=x.max()/2, **kw2 )
  cb = fig.colorbar(plot_contour2, ax=ax)
  cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])

  ax.set_xlabel("x (m)")
  ax.set_ylabel("z (m)")
  ax.set_zlabel("y (m)")
  ax.set_title("FOM {:}, t={:.1f}".format(varname, t_solution))
  ax.set_xlim(np.min(x3d), np.max(x3d))
  ax.set_ylim(np.min(y3d), np.max(y3d))
  ax.set_zlim(np.min(z3d), np.max(z3d))
  ax.view_init(10, 210)
  ax.dist = 9

  kw1 =  {'vmin': data_rom_3d.min(),
          'vmax': data_rom_3d.max(),
          'levels': np.linspace((data_rom_3d.min()+0.1*(data_rom_3d.max()-data_rom_3d.min())),
                                (data_rom_3d.max()-0.1*(data_rom_3d.max()-data_rom_3d.min())),
                                25),
          'antialiased' : True,
          'alpha' : 0.8,
          'cmap' : colormap,
        }
  kw2 =  {'vmin': data_rom_3d.min(),
          'vmax': data_rom_3d.max(),
          'levels': np.linspace(data_rom_3d.min(),
                                data_rom_3d.max(),
                                100),
          'antialiased' : True,
          'alpha' : 0.6,
          'cmap' : colormap,
        }
  ax = fig.add_subplot(fig_nv, fig_nh, 2, projection='3d')
  plot_contour1= ax.contour(  x3d[:, 0, :],
                              data_rom_3d[:,int(slice_loc*data_rom_3d.shape[1]),:],
                              z3d[:, 0, :],
                              zdir='y', offset=y.max()/2, **kw1 )
  plot_contour2= ax.contourf( data_rom_3d[int(slice_loc*data_rom_3d.shape[0]), :, :],
                              y3d[0, :, :],
                              z3d[0, :, :],
                              zdir='x', offset=x.max()/2, **kw2 )
  cb = fig.colorbar(plot_contour2, ax=ax)
  cb.ax.set_yticklabels(["{:.2f}".format(i) for i in cb.get_ticks()])
  ax.set_xlabel("x (m)")
  ax.set_ylabel("z (m)")
  ax.set_zlabel("y (m)")
  ax.set_title("ROM {:}, t={:.1f}".format(varname, t_solution))
  ax.set_xlim(np.min(x3d), np.max(x3d))
  ax.set_ylim(np.min(y3d), np.max(y3d))
  ax.set_zlim(np.min(z3d), np.max(z3d))
  ax.view_init(10, 210)
  ax.dist = 9

  kw1 =  {'vmin': data_del_3d.min(),
          'vmax': data_del_3d.max(),
          'levels': np.linspace((data_del_3d.min()+0.1*(data_del_3d.max()-data_del_3d.min())),
                                (data_del_3d.max()-0.1*(data_del_3d.max()-data_del_3d.min())),
                                25),
          'antialiased' : True,
          'alpha' : 0.8,
          'cmap' : colormap,
        }
  kw2 =  {'vmin': data_del_3d.min(),
          'vmax': data_del_3d.max(),
          'levels': np.linspace(data_del_3d.min(),
                                data_del_3d.max(),
                                100),
          'antialiased' : True,
          'alpha' : 0.6,
          'cmap' : colormap,
        }
  ax = fig.add_subplot(fig_nv, fig_nh, 3, projection='3d')
  plot_contour1= ax.contour(  x3d[:, 0, :],
                              data_del_3d[:,int(slice_loc*data_del_3d.shape[1]),:],
                              z3d[:, 0, :],
                              zdir='y', offset=y.max()/2, **kw1 )
  plot_contour2= ax.contourf( data_del_3d[int(slice_loc*data_del_3d.shape[0]), :, :],
                              y3d[0, :, :],
                              z3d[0, :, :],
                              zdir='x', offset=x.max()/2, **kw2 )
  cb = fig.colorbar(plot_contour2, ax=ax)
  cb.ax.set_yticklabels(["{:.1e}".format(i) for i in cb.get_ticks()])
  ax.set_xlabel("x (m)")
  ax.set_ylabel("z (m)")
  ax.set_zlabel("y (m)")
  ax.set_title("Diff {:}, t={:.1f}, norm={:.3e}".format(varname,
                                                      t_solution,
                                                      diff_norm))
  ax.set_xlim(np.min(x3d), np.max(x3d))
  ax.set_ylim(np.min(y3d), np.max(y3d))
  ax.set_zlim(np.min(z3d), np.max(z3d))
  ax.view_init(10, 210)
  ax.dist = 9

  plt.tight_layout()
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
        fname = sim_path + '/'+op_root+'_'+f'{sim:03d}'+'_'+f'{i:05d}'+'.bin'
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
        fname = sim_path + '/'+op_root+'_'+f'{sim:03d}'+'_'+f'{i:05d}'+'.bin'
      elif nsims >= 10:
        fname = sim_path + '/'+op_root+'_'+f'{sim:02d}'+'_'+f'{i:05d}'+'.bin'
      elif nsims > 1:
        fname = sim_path + '/'+op_root+'_'+f'{sim:01d}'+'_'+f'{i:05d}'+'.bin'
      else:
        fname = sim_path + '/'+op_root+'_'+f'{i:05d}'+'.bin'
      grid,solution_rom = hyparutils.readOpFile(fname, ndims, nvars, size)
      solution_rom = np.float32(solution_rom)

      x = grid[:size[0]]
      y = grid[size[0]:(size[0]+size[1])]
      z = grid[(size[0]+size[1]):]
      y3d, x3d, z3d = np.meshgrid(x,y,z)

      u_fom_cons_3d = np.transpose(solution_fom.reshape(size[2],size[1],size[0],nvars))
      u_fom_prim_3d = convertVar(u_fom_cons_3d, y3d)

      u_rom_cons_3d = np.transpose(solution_rom.reshape(size[2],size[1],size[0],nvars))
      u_rom_prim_3d = convertVar(u_rom_cons_3d, y3d)

      theta_fom = u_fom_prim_3d[5,:,:,:]
      theta_rom = u_rom_prim_3d[5,:,:,:]
      if nsims > 1:
        plt_fname = plt_dir_name+'/fig_theta_'+f'{s:02d}'+'_'+f'{i:05d}'+'.png'
      else:
        plt_fname = plt_dir_name+'/fig_theta_'+f'{i:05d}'+'.png'
      plotData( theta_fom,
                theta_rom,
                x3d, y3d, z3d,
                plt_fname,
                'theta',
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
    y = grid[size[0]:(size[0]+size[1])]
    z = grid[(size[0]+size[1]):]
    y3d, x3d, z3d = np.meshgrid(x,y,z)

    u_fom_cons_3d = np.transpose(solution_fom.reshape(size[2],size[1],size[0],nvars))
    u_fom_prim_3d = convertVar(u_fom_cons_3d, y3d)

    u_rom_cons_3d = np.transpose(solution_rom.reshape(size[2],size[1],size[0],nvars))
    u_rom_prim_3d = convertVar(u_rom_cons_3d, y3d)

    theta_fom = u_fom_prim_3d[5,:,:,:]
    theta_rom = u_rom_prim_3d[5,:,:,:]
    if nsims > 1:
      plt_fname = plt_dir_name+'/fig_theta_'+f'{s:02d}'+'.png'
    else:
        plt_fname = plt_dir_name+'/fig_theta.png'
    plotData( theta_fom,
              theta_rom,
              x3d, y3d, z3d,
              plt_fname, 'theta',
              t_final )

