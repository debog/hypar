import numpy as np

'''
Read in a HyPar text input file (eg, solver.inp, physics.inp, etc)

Input arguments:
  fname: name of the file (eg: "solver.inp")

The function will return a dictionary with the keywords found
in the file.
'''
def readHyParInpFile(fname: str):
  input_data = {}

  try:
    with open(fname,'r') as solver_file:
      begin = False
      for line in solver_file:
        if line.startswith("end"):
          break
        elif begin:
          fields = line.split()
          input_data[fields[0]] = np.array(fields[1:])
        elif line.startswith("begin"):
          begin = True
  except:
    print('File ', fname, ' is absent.')

  return input_data

def readHyParErrFile( path: str,
                      ndims: int):

  fname = path + '/errors.dat'
  with open(fname,'r') as file:
    for line in file:
      words = line.split()
      total_wctime = float(words[-1])
      solver_wctime = float(words[-2])
      errors = np.float32(np.array(words[-5:-3]))
  return errors, solver_wctime, total_wctime


'''
Read in a binary solution output file

Input arguments:
  fname: name of file (eg. "op.bin", "op_00010.bin")
  a_ndims: number of spatial dimensions
  a_nvars: number of vector components at each grid point
  a_size: integer array with grid size in each dimension

Returns a tuple with the grid x and the solution u.
'''
def readOpFile( fname: str,
                a_ndims: int,
                a_nvars: int,
                a_size: np.array):
  with open(fname, "rb") as f_op:
    ndims = np.fromfile(f_op, dtype=np.int32, count=1)[0]
    if ndims != a_ndims:
      raise Exception("ndims did not match in op file!")
    nvars = np.fromfile(f_op, dtype=np.int32, count=1)[0]
    if nvars != a_nvars:
      raise Exception("nvars did not match in op file!")
    size = np.fromfile(f_op, dtype=np.int32, count=ndims)
    if np.any(np.not_equal(size,a_size)):
      print('size in op file: ', size)
      print('specified size: ', a_size)
      raise Exception("size did not match in op file!")

    x = np.fromfile(f_op, dtype=np.float64, count=np.sum(size))
    u = np.fromfile(f_op, dtype=np.float64, count=np.prod(size)*nvars)

    return x, u

'''
For an unsteady simulation, this function reads in the solution
from each op_sim_idx_<index>.bin file and creates a 2D matrix whose rows
are the serialized solution vector at a particular simulation time.

Input arguments:
  path        : directory path where the op_* files are located
  nsims       : total number of simulation domains
  sim_idx     : the index of the simulation for which to
                read the solution files (0 <= sim_idx < nsims)
  n_op_files  : number of snapshot files to read
  ndims       : number of spatial dimensions
  nvars       : number of vector components at each grid point
  size        : integer array with grid size in each dimension
(Optional)
  op_root     : filename root for output files
'''
def getSimulationSnapshots( path: str,
                            nsims: int,
                            sim_idx: int,
                            n_op_files: int,
                            ndims: int,
                            nvars: int,
                            size: np.ndarray,
                            op_root: str ='op'):
  ndof = nvars * np.prod(size)
  snapshots = np.empty((0,ndof),np.float64)
  if n_op_files > 1:
    for i in range(n_op_files):
      if nsims >= 100:
        fname = path + '/'+op_root+'_'+f'{sim_idx:03d}'+'_'+f'{i:05d}'+'.bin'
      elif nsims >= 10:
        fname = path + '/'+op_root+'_'+f'{sim_idx:02d}'+'_'+f'{i:05d}'+'.bin'
      elif nsims > 1:
        fname = path + '/'+op_root+'_'+f'{sim_idx:01d}'+'_'+f'{i:05d}'+'.bin'
      else:
        fname = path + '/'+op_root+'_'+f'{i:05d}'+'.bin'

      try:
        x, u = readOpFile(fname, ndims, nvars, size)
        snapshots = np.concatenate((snapshots,np.expand_dims(u,axis=0)),axis=0)
      except:
        pass

  else:
    if nsims >= 100:
      fname = path + '/'+op_root+'_'+f'{sim_idx:03d}'+'.bin'
    elif nsims >= 10:
      fname = path + '/'+op_root+'_'+f'{sim_idx:02d}'+'.bin'
    elif nsims > 1:
      fname = path + '/'+op_root+'_'+f'{sim_idx:01d}'+'.bin'
    else:
      fname = path + '/'+op_root+'.bin'

    try:
      x, u = readOpFile(fname, ndims, nvars, size)
      snapshots = np.concatenate((snapshots,np.expand_dims(u,axis=0)),axis=0)
    except:
      pass

  return x,snapshots

'''
For an unsteady simulation, this function reads in the solution
from each op_<index>.bin file and creates a 2D matrix whose rows
are the serialized solution vector at a particular simulation time.

Input arguments:
  path        : directory path where the op_* files are located
  nsims       : total number of simulation domains
  n_op_files  : number of snapshot files to read
  ndims       : number of spatial dimensions
  nvars       : number of vector components at each grid point
  size        : integer array with grid size in each dimension
(Optional)
  op_root     : filename root for output files

NOTE: for nsims > 1, the grid size for each simulation domain
must be exactly the same!
'''
def getSolutionSnapshots( path: str,
                          nsims: int,
                          n_op_files: int,
                          ndims: int,
                          nvars: int,
                          size: np.ndarray,
                          op_root: str ='op' ):
  if nsims > 1:
    ndof = nvars * np.prod(size[0,:])
  else:
    ndof = nvars * np.prod(size)
  snapshots = np.empty((0,ndof),np.float64)
  for sim in range(nsims):
    if nsims > 1:
      size_sim = size[sim,:]
    else:
      size_sim = size
    if n_op_files > 1:
      for i in range(n_op_files):
        if nsims >= 100:
          fname = path + '/'+op_root+'_'+f'{sim:03d}'+'_'+f'{i:05d}'+'.bin'
        elif nsims >= 10:
          fname = path + '/'+op_root+'_'+f'{sim:02d}'+'_'+f'{i:05d}'+'.bin'
        elif nsims > 1:
          fname = path + '/'+op_root+'_'+f'{sim:01d}'+'_'+f'{i:05d}'+'.bin'
        else:
          fname = path + '/'+op_root+'_'+f'{i:05d}'+'.bin'

        try:
          x, u = readOpFile(fname, ndims, nvars, size_sim)
          snapshots = np.concatenate((snapshots,np.expand_dims(u,axis=0)),axis=0)
        except:
          pass

    else:
      if nsims >= 100:
        fname = path + '/'+op_root+'_'+f'{sim:03d}'+'.bin'
      elif nsims >= 10:
        fname = path + '/'+op_root+'_'+f'{sim:02d}'+'.bin'
      elif nsims > 1:
        fname = path + '/'+op_root+'_'+f'{sim:01d}'+'.bin'
      else:
        fname = path + '/'+op_root+'.bin'

      try:
        x, u = readOpFile(fname, ndims, nvars, size_sim)
        snapshots = np.concatenate((snapshots,np.expand_dims(u,axis=0)),axis=0)
      except:
        pass

  return x,snapshots

