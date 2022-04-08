'''
Import whatever needs to be imported
'''

#Enable interactive plot
#%matplotlib notebook

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

print('Number of simulations: ', nsims)

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
print('  simulation dt = ', dt)
print('  niter = ', niter)
print('  final time = ', t_final)
print('  snapshot dt = ', dt_snapshots)
print('  number of snapshots = ', n_snapshots)

'''
Load simulation data (solution snapshots)
'''

x,solution_snapshots = hyparutils.getSolutionSnapshots(train_sim_path, nsims, n_snapshots, ndims, nvars, size)
solution_snapshots = np.float32(solution_snapshots)
print('Solution snapshot shape: ', solution_snapshots.shape)

def animatedensity(i):
    for s in range(nsims):
        solution_snapshots_sim = solution_snapshots[s*n_snapshots:(s+1)*n_snapshots]
        lines[s].set_ydata(solution_snapshots_sim[i, 0::3])

def animatemomentum(i):
    for s in range(nsims):
        solution_snapshots_sim = solution_snapshots[s*n_snapshots:(s+1)*n_snapshots]
        lines[s].set_ydata(solution_snapshots_sim[i, 1::3])

def animateenergy(i):
    for s in range(nsims):
        solution_snapshots_sim = solution_snapshots[s*n_snapshots:(s+1)*n_snapshots]
        lines[s].set_ydata(solution_snapshots_sim[i, 2::3])

fig = plt.figure(figsize=(10,6))
ax = plt.axes()
ax.set(xlim=(np.min(x), np.max(x)), ylim=(np.min(solution_snapshots[:,0::3]), np.max(solution_snapshots[:,0::3])))
lines = [None]*nsims
for s in range(nsims):
    solution_snapshots_sim = solution_snapshots[s*n_snapshots:(s+1)*n_snapshots]
    lines[s] = ax.plot(x, solution_snapshots_sim[0, 0::3], lw=1)[0]

anim = FuncAnimation(fig, animatedensity, interval=10, frames=n_snapshots, blit=False)
anim.save('training_simulation_density.mp4')

fig = plt.figure(figsize=(10,6))
ax = plt.axes()
ax.set(xlim=(np.min(x), np.max(x)), ylim=(np.min(solution_snapshots[:,1::3]), np.max(solution_snapshots[:,1::3])))
lines = [None]*nsims
for s in range(nsims):
    solution_snapshots_sim = solution_snapshots[s*n_snapshots:(s+1)*n_snapshots]
    lines[s] = ax.plot(x, solution_snapshots_sim[0, 1::3], lw=1)[0]

anim = FuncAnimation(fig, animatemomentum, interval=10, frames=n_snapshots, blit=False)
anim.save('training_simulation_momentum.mp4')

fig = plt.figure(figsize=(10,6))
ax = plt.axes()
ax.set(xlim=(np.min(x), np.max(x)), ylim=(np.min(solution_snapshots[:,2::3]), np.max(solution_snapshots[:,2::3])))
lines = [None]*nsims
for s in range(nsims):
    solution_snapshots_sim = solution_snapshots[s*n_snapshots:(s+1)*n_snapshots]
    lines[s] = ax.plot(x, solution_snapshots_sim[0, 2::3], lw=1)[0]

anim = FuncAnimation(fig, animateenergy, interval=10, frames=n_snapshots, blit=False)
anim.save('training_simulation_energy.mp4')

# define testset and trainset indices
test_frac = 0.4
print('Test fraction: ', test_frac)
test_ind = np.random.permutation(np.arange(solution_snapshots.shape[0]))[:int(test_frac*solution_snapshots.shape[0])]
train_ind = np.setdiff1d(np.arange(solution_snapshots.shape[0]),test_ind)

# set trainset and testset
trainset = solution_snapshots[train_ind]
testset = solution_snapshots[test_ind]
print('Training data shape: ', trainset.shape)
print('Test data shape: ', testset.shape)

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

# set batch_size, number of epochs, patience for early stop
batch_size = 20
num_epochs = 600
num_epochs_print = 10

# filename
AE_fname = 'AE_git.tar'

encoder, decoder = autoencoder.createAE(encoder_class,
                                        decoder_class,
                                        f_activation,
                                        mask,
                                        m, f, M1, M2,
                                        device )

autoencoder.trainAE(encoder,
                    decoder,
                    trainset,
                    testset,
                    batch_size,
                    num_epochs,
                    num_epochs_print,
                    AE_fname )
