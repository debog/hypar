'''
Import whatever needs to be imported
'''

import os
import time
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy.linalg as LA
from itertools import product

hypar_dir = os.environ.get('HYPAR_DIR')
hypar_dir_python = hypar_dir + '/Examples/Python'
sys.path.append(hypar_dir_python)

lasdi_dir = os.environ.get('LASDI_DIR')
sys.path.append(lasdi_dir)

import modHyParUtils as hyparutils
import modLaSDIUtils as lasdiutils
import modAutoEncoder as autoencoder

from LaSDI import LaSDI

'''
Set up the simulation parameters
'''
train_sim_path = 'training_data'
sim_inp_data = hyparutils.readHyParInpFile(train_sim_path+'/simulation.inp')
solver_inp_data = hyparutils.readHyParInpFile(train_sim_path+'/solver.inp')

if not sim_inp_data:
    nsims = 1
else:
    nsims = int(sim_inp_data['nsims'][0])

print('Number of training simulations: ', nsims)

if solver_inp_data['op_overwrite'] != 'no':
    raise Exception("op_overwrite must be 'no' in solver.inp.")

if solver_inp_data['op_file_format'] != 'binary':
    raise Exception("op_file_format must be 'binary' in solver.inp.")


ndims = int(solver_inp_data['ndims'][0])
nvars = int(solver_inp_data['nvars'][0])
size = np.int32(solver_inp_data['size'][0])
ndof = np.product(size)*nvars
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
print('  ndof = ', ndof)
print('  (training data) simulation dt = ', dt)
print('  (training data) niter = ', niter)
print('  (training data) final time = ', t_final)
print('  (training data) snapshot dt = ', dt_snapshots)
print('  (training data) number of snapshots = ', n_snapshots)

param_data = hyparutils.readHyParInpFile(train_sim_path+'/parameters.dat')
rhoL = np.float32(param_data['rhoL'])
rhoR = np.float32(param_data['rhoR'])
pL = np.float32(param_data['pL'])
pR = np.float32(param_data['pR'])
P = list(product(rhoL,rhoR,pL,pR))
P = np.array(P)
print('Parameter space:\n', P)

'''
Load simulation data (solution snapshots)
'''

x,solution_snapshots = hyparutils.getSolutionSnapshots(train_sim_path, nsims, n_snapshots, ndims, nvars, size)
solution_snapshots = np.float32(solution_snapshots)
print('Training data solution snapshot shape: ', solution_snapshots.shape)

# set device
device = autoencoder.getDevice()
print('Using device: ', device)

# set encoder and decoder types, activation function, etc.
encoder_class = autoencoder.Encoder
decoder_class = autoencoder.Decoder
f_activation = autoencoder.SiLU

# set the number of nodes in each layer
m = ndof
f = 12
b = 36
db = 12
M2 = b + (m-1)*db
M1 = 2*m
mask = lasdiutils.create_mask_1d(m,b,db)

'''
Generate LaSDI-NM Model
'''
AE_fname = 'AE_git.tar'
encoder, decoder = autoencoder.readAEFromFile(  encoder_class,
                                                decoder_class,
                                                f_activation,
                                                mask,
                                                m, f, M1, M2,
                                                device,
                                                AE_fname )

latent_space_snapshots = autoencoder.encodedSnapshots(  encoder,
                                                        solution_snapshots,
                                                        niter+1,
                                                        device )
LaSDI_model = LaSDI(encoder, decoder, NN = True, device = device)
LaSDI_model.train_dynamics(latent_space_snapshots, P, dt_snapshots)

'''
Load test simulation data (solution snapshots)
'''
test_sim_path = 'test'
solver_inp_data = hyparutils.readHyParInpFile(test_sim_path+'/solver.inp')

if solver_inp_data['op_overwrite'] != 'no':
    raise Exception("op_overwrite must be 'no' in solver.inp.")

if solver_inp_data['op_file_format'] != 'binary':
    raise Exception("op_file_format must be 'binary' in solver.inp.")

niter = int(solver_inp_data['n_iter'])
dt = float(solver_inp_data['dt'])
t_final = dt*niter
op_write_iter = int(solver_inp_data['file_op_iter'])
dt_snapshots = op_write_iter*dt
n_snapshots = int(niter/op_write_iter) + 1

print('Test simulation parameters:')
print('  simulation dt = ', dt)
print('  niter = ', niter)
print('  final time = ', t_final)
print('  snapshot dt = ', dt_snapshots)
print('  number of snapshots = ', n_snapshots)

param_data = hyparutils.readHyParInpFile(test_sim_path+'/parameters.dat')
rhoL = np.float32(param_data['rhoL'])
rhoR = np.float32(param_data['rhoR'])
pL = np.float32(param_data['pL'])
pR = np.float32(param_data['pR'])
P = list(product(rhoL,rhoR,pL,pR))
P = np.squeeze(np.array(P))
print('Test sim parameters:\n', P)

x,solution_snapshots_FOM = hyparutils.getSolutionSnapshots(test_sim_path, 1, n_snapshots, ndims, nvars, size)
solution_snapshots_FOM = np.float32(solution_snapshots_FOM)

FOM_err, FOM_time, FOM_totaltime = hyparutils.readHyParErrFile(test_sim_path, ndims)

t_arr = np.linspace(0, t_final, int(niter/op_write_iter)+1)

start = time.time()
FOM_recon = LaSDI_model.generate_ROM( solution_snapshots_FOM[0],
                                      P,
                                      t_arr )
LaSDI_time = time.time()-start

print('FOM wctime:', FOM_time)
print('LaSDI wctime: ', LaSDI_time)

fig = plt.figure()
ax = plt.axes()
fig.suptitle('Reconstruction of FOM from LaSDI-NM', y = 1.05)
ax.set_title('Relative Error: {:.3}  Speedup: {:.3} times'.format(LA.norm(FOM_recon[-1,::3]-solution_snapshots_FOM[-1,::3])/LA.norm(solution_snapshots_FOM[-1,::3]), FOM_time/LaSDI_time))
ax.plot(x,solution_snapshots_FOM[-1,::3], label = 'FOM')
ax.plot(x, FOM_recon[-1,::3],'--', label = 'LaSDI')
ax.legend()
ax.set_xlim(np.min(x),np.max(x))
plt.savefig('FOM_vs_LaSDI_solution_density.png')

fig = plt.figure()
ax = plt.axes()
fig.suptitle('Reconstruction of FOM from LaSDI-NM', y = 1.05)
ax.set_title('Relative Error: {:.3}  Speedup: {:.3} times'.format(LA.norm(FOM_recon[-1,1::3]-solution_snapshots_FOM[-1,1::3])/LA.norm(solution_snapshots_FOM[-1,1::3]), FOM_time/LaSDI_time))
ax.plot(x,solution_snapshots_FOM[-1,1::3], label = 'FOM')
ax.plot(x, FOM_recon[-1,1::3],'--', label = 'LaSDI')
ax.legend()
ax.set_xlim(np.min(x),np.max(x))
plt.savefig('FOM_vs_LaSDI_solution_momentum.png')

fig = plt.figure()
ax = plt.axes()
fig.suptitle('Reconstruction of FOM from LaSDI-NM', y = 1.05)
ax.set_title('Relative Error: {:.3}  Speedup: {:.3} times'.format(LA.norm(FOM_recon[-1,2::3]-solution_snapshots_FOM[-1,2::3])/LA.norm(solution_snapshots_FOM[-1,2::3]), FOM_time/LaSDI_time))
ax.plot(x,solution_snapshots_FOM[-1,2::3], label = 'FOM')
ax.plot(x, FOM_recon[-1,2::3],'--', label = 'LaSDI')
ax.legend()
ax.set_xlim(np.min(x),np.max(x))
plt.savefig('FOM_vs_LaSDI_solution_energy.png')

FOM_re = np.empty([3,n_snapshots])
for i in range(n_snapshots):
    FOM_re[0,i] = LA.norm(FOM_recon[i,0::3]-solution_snapshots_FOM[i,0::3])/LA.norm(solution_snapshots_FOM[i,0::3])
    if LA.norm(solution_snapshots_FOM[i,1::3]) > 0:
      FOM_re[1,i] = LA.norm(FOM_recon[i,1::3]-solution_snapshots_FOM[i,1::3])/LA.norm(solution_snapshots_FOM[i,1::3])
    else:
      FOM_re[1,i] = LA.norm(FOM_recon[i,1::3]-solution_snapshots_FOM[i,1::3])
    FOM_re[2,i] = LA.norm(FOM_recon[i,2::3]-solution_snapshots_FOM[i,2::3])/LA.norm(solution_snapshots_FOM[i,2::3])

fig = plt.figure()
fig.suptitle('Relative Error for FOM reconstruction via LaSDI-NM', y = 1.05)
ax = plt.axes()
ax.set_title('Max Relative Error: {:.3}'.format(np.amax(FOM_re)))
ax.plot(t_arr, FOM_re[0,:], label='density')
ax.plot(t_arr, FOM_re[1,:], label='momentum')
ax.plot(t_arr, FOM_re[2,:], label='energy')
ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Relative Error')
#plt.show()
plt.savefig('LaSDI_error_vs_time.png')
