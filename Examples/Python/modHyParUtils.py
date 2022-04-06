import numpy as np

'''
Read in a HyPar text input file (eg, solver.inp, physics.inp, etc)

Input arguments:
  fname: name of the file (eg: "solver.inp")

The function will return a dictionary with the keywords found
in the file.
'''
def readHyParInpFile(fname):
  input_data = {}
  
  with open(fname,'r') as solver_file:
    begin = False
    for line in solver_file:
      if line.startswith("end"):
        break
      elif begin:
        fields = line.split()
        input_data[fields[0]] = fields[1]
      elif line.startswith("begin"):
        begin = True
        
  return input_data

'''
Read in a binary solution output file

Input arguments:
  fname: name of file (eg. "op.bin", "op_00010.bin")
  a_ndims: number of spatial dimensions
  a_nvars: number of vector components at each grid point
  a_size: integer array with grid size in each dimension

Returns a tuple with the grid x and the solution u.
'''
def readOpFile(fname, a_ndims, a_nvars, a_size):
  with open(fname, "rb") as f_op:
    ndims = np.fromfile(f_op, dtype=np.int32, count=1)[0]
    if ndims != a_ndims:
      raise Exception("ndims did not match in op file!")
    nvars = np.fromfile(f_op, dtype=np.int32, count=1)[0]
    if nvars != a_nvars:
      raise Exception("nvars did not match in op file!")
    size = np.fromfile(f_op, dtype=np.int32, count=ndims)
    if np.array_equal(size,a_size):
      raise Exception("size did not match in op file!")

    x = np.fromfile(f_op, dtype=np.float64, count=np.sum(size))
    u = np.fromfile(f_op, dtype=np.float64, count=np.prod(size)*nvars)

    return x, u
'''
For an unsteady simulation, this function reads in the solution
from each op_<index>.bin file and creates a 2D matrix whose rows
are the serialized solution vector at a particular simulation time.
'''
def getSolutionSnapshots(n_op_files, ndims, nvars, size):
  ndof = nvars * np.prod(size)
  snapshots = np.empty((0,ndof),np.float64)
  for i in range(n_op_files):
    fname = 'op_'+f'{i:05d}'+'.bin'
    x, u = readOpFile(fname, ndims, nvars, size)
    snapshots = np.concatenate((snapshots,np.expand_dims(u,axis=0)),axis=0)

  return snapshots

