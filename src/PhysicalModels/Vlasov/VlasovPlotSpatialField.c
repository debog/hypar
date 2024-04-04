/*! @file VlasovPlotSpatialField.c
    @author Debojyoti Ghosh
    @brief Plot out the field to file
*/

#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <io.h>
#include <mpivars.h>
#include <hypar.h>
#include <physicalmodels/vlasov.h>

#ifdef with_python
#include <Python.h>
#ifdef with_python_numpy
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#endif
#endif

/*! Plot out a spatial field or variable to file */
int VlasovPlotSpatialField( void*   s,         /*!< Solver object of type #HyPar */
                            void*   m,         /*!< MPI object of type #MPIVariables */
                            double* a_field,   /*!< Vector field to write */
                            double  a_time,    /*!< Current simulation time */
                            char*   fname_root /*!< Filename root (extension is added automatically).
                                                     For unsteady output, a numerical index is added
                                                     that is the same as for the solution output files. */)
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  Vlasov        *param  = (Vlasov*)       solver->physics;

  if (solver->nsims > 1) {
    char index[_MAX_STRING_SIZE_];
    GetStringFromInteger(solver->my_idx, index, (int)log10(solver->nsims)+1);
    strcat(fname_root, "_");
    strcat(fname_root, index);
    strcat(fname_root, "_");
  }

  char filename[_MAX_STRING_SIZE_] = "";
  strcat(filename,fname_root);
  if (!strcmp(solver->op_overwrite,"no")) {
    strcat(filename,"_");
    strcat(filename,solver->filename_index);
  }
  strcat(filename,".png");

  int d,
      ghosts = solver->ghosts,
      ndims_x = param->ndims_x;

  int dim_global_x[ndims_x];
  _ArrayCopy1D_(solver->dim_global, dim_global_x, ndims_x);
  int dim_local_x[ndims_x];
  _ArrayCopy1D_(solver->dim_local, dim_local_x, ndims_x);

  /* gather the spatial coordinates into a global array */
  double *xg;
  {
    int size_g = param->npts_global_x;
    xg = (double*) calloc (size_g, sizeof(double));
    _ArraySetValue_(xg, size_g, 0.0);

    int offset_global, offset_local;
    offset_global = offset_local = 0;
    for (d=0; d<ndims_x; d++) {
      IERR MPIGatherArray1D(  mpi,
                              (mpi->rank?NULL:&xg[offset_global]),
                              &solver->x[offset_local+ghosts],
                              mpi->is[d],
                              mpi->ie[d],
                              solver->dim_local[d],
                              0);
      offset_global += dim_global_x[d];
      offset_local  += dim_local_x [d] + 2*ghosts;
    }
  }

  /* gather the field into a global array */
  double *field_g;
  {
    int size_g = param->npts_global_x * param->ndims_x;
    field_g = (double*) calloc (size_g, sizeof(double));
    _ArraySetValue_(field_g, size_g, 0.0);

    if (param->ndims_x > 1) {
      if (!mpi->rank) {
        fprintf(stderr,"Warning in VlasovPlotSpatialField():\n");
        fprintf(stderr,"  field plotting not yet supported for >1 spatial dimensions.\n");
      }
    } else {
      IERR MPIGatherArray1D(  mpi,
                              (mpi->rank ? NULL : field_g),
                              (a_field+ghosts),
                              mpi->is[0],
                              mpi->ie[0],
                              solver->dim_local[0],
                              0);
    }
  }

  if (!mpi->rank) {
#ifdef with_python
#ifdef with_python_numpy
    import_array();
    PyObject* py_plt_func = (PyObject*) solver->py_plt_func;
    PyObject* py_plt_func_args = (PyObject*) solver->py_plt_func_args;
    py_plt_func_args = PyTuple_New(7);
    {
      PyObject* py_obj = Py_BuildValue("i", ndims_x);
      PyTuple_SetItem(py_plt_func_args, 0, py_obj);
    }
    {
      PyObject* py_obj = Py_BuildValue("i", 1);
      PyTuple_SetItem(py_plt_func_args, 1, py_obj);
    }
    {
      npy_intp shape[1] = {ndims_x};
      PyObject* size_arr = PyArray_SimpleNewFromData(1,shape,NPY_INT,dim_global_x);
      PyTuple_SetItem(py_plt_func_args, 2, size_arr);
    }
    {
      PyObject* py_obj = Py_BuildValue("d", a_time);
      PyTuple_SetItem(py_plt_func_args, 3, py_obj);
    }
    {
      npy_intp shape[1] = {param->npts_global_x};
      PyObject* x_arr = PyArray_SimpleNewFromData(1,shape,NPY_DOUBLE,xg);
      PyTuple_SetItem(py_plt_func_args, 4, x_arr);
    }
    {
      npy_intp shape[1] = {param->npts_global_x * param->ndims_x};
      PyObject* u_arr = PyArray_SimpleNewFromData(1,shape,NPY_DOUBLE,field_g);
      PyTuple_SetItem(py_plt_func_args, 5, u_arr);
    }
    {
      PyObject* py_obj = Py_BuildValue("s", filename);
      PyTuple_SetItem(py_plt_func_args, 6, py_obj);
    }
    if (!py_plt_func) {
      fprintf(stderr,"Error in PlotArraySerial(): solver->py_plt_func is NULL!\n");
    } else {
      PyObject_CallObject(py_plt_func, py_plt_func_args);
    }
#else
    fprintf(stderr,"Error in PlotArraySerial(): HyPar needs to be compiled with Numpy support (-Dwith_python_numpy) to use this feature..\n");
#endif
#else
    fprintf(stderr,"Error in PlotArraySerial(): HyPar needs to be compiled with Python support (-Dwith_python) to use this feature..\n");
#endif
  }

  /* free up arrays */
  free(xg);
  free(field_g);

  return 0;
}
