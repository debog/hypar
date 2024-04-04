/*! @file PlotArray.c
    @author Debojyoti Ghosh
    @brief Plot a vector field, stored as an array, and save to file

    Contains functions to plot out a vector field, stored
    as an array and save the figure to a file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <plotfuncs.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef with_python
#include <Python.h>
#ifdef with_python_numpy
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#endif
#endif

/* Function declarations */
static int PlotArraySerial   (int,int,int*,int*,int,double*,double*,double,void*,void*,char*);

/*! Plot a vector field, stored as an array, and save figure to file: this is a
    wrapper function that calls PlotArraySerial(). */
int PlotArray(int     a_ndims,      /*!< Number of spatial dimensions */
              int     a_nvars,      /*!< Number of variables per grid point */
              int*    a_dim_global, /*!< Integer array of size a_ndims with global grid size in each dimension */
              int*    a_dim_local,  /*!< Integer array of size a_ndims with local  grid size in each dimension */
              int     a_ghosts,     /*!< Number of ghost points */
              double* a_x,          /*!< Array of spatial coordinates (i.e. the grid) */
              double* a_u,          /*!< Vector field to write */
              double  a_time,       /*!< Current simulation time */
              void*   a_s,          /*!< Solver object of type #HyPar */
              void*   a_m,          /*!< MPI object of type #MPIVariables */
              char*   a_fname_root  /*!< Filename root (extension is added automatically). For unsteady output,
                                         a numerical index is added that is the same as for the solution output files. */ )
{
  PlotArraySerial(  a_ndims,
                    a_nvars,
                    a_dim_global,
                    a_dim_local,
                    a_ghosts,
                    a_x,
                    a_u,
                    a_time,
                    a_s,
                    a_m,
                    a_fname_root  );

  return 0;
}

/*! Function to plot a vector field, stored as an array, and
    save the figure to a file. It will allocate the global domain on
    rank 0, so do not use for big problems for which the entire global
    domain will not fit on one node. This approach is also not very
    scalable.
    + Needs HyPar to be compiled with Python.
*/
int PlotArraySerial(int     a_ndims,      /*!< Number of spatial dimensions */
                    int     a_nvars,      /*!< Number of variables per grid point */
                    int*    a_dim_global, /*!< Integer array of size a_ndims with global grid size in each dimension */
                    int*    a_dim_local,  /*!< Integer array of size a_ndims with local  grid size in each dimension */
                    int     a_ghosts,     /*!< Number of ghost points */
                    double* a_x,          /*!< Array of spatial coordinates (i.e. the grid) */
                    double* a_u,          /*!< Vector field to write */
                    double  a_time,       /*!< Current simulation time */
                    void*   a_s,          /*!< Solver object of type #HyPar */
                    void*   a_m,          /*!< MPI object of type #MPIVariables */
                    char*   a_fname_root  /*!< Filename root (extension is added automatically). For unsteady output,
                                               a numerical index is added that is the same as for the solution output files. */ )
{
  HyPar         *solver = (HyPar*)       a_s;
  MPIVariables  *mpi    = (MPIVariables*)a_m;

  int size_global_x = 0;
  int size_global_u = 0;

  /* root process: allocate global output arrays */
  double *ug, *xg;
  if (!mpi->rank) {

    size_global_u = 1;
    for (int d=0; d<a_ndims; d++) size_global_u *= a_dim_global[d];
    ug = (double*) calloc (size_global_u*a_nvars,sizeof(double));
    _ArraySetValue_(ug,size_global_u*a_nvars,0.0);

    size_global_x = 0;
    for (int d=0; d<a_ndims; d++) size_global_x += a_dim_global[d];
    xg = (double*) calloc (size_global_x,sizeof(double));
    _ArraySetValue_(xg,size_global_x,0.0);

  } else {

    /* null pointers on non-root processes */
    ug = xg = NULL;

  }

  /* Assemble the local output arrays into the global output arrays */
  MPIGatherArraynD( a_ndims,
                    mpi,
                    ug,
                    a_u,
                    a_dim_global,
                    a_dim_local,
                    a_ghosts,
                    a_nvars );
  int offset_global, offset_local;
  offset_global = offset_local = 0;
  for (int d=0; d<a_ndims; d++) {
    MPIGatherArray1D( mpi,
                      (mpi->rank?NULL:&xg[offset_global]),
                      &a_x[offset_local+a_ghosts],
                      mpi->is[d],
                      mpi->ie[d],
                      a_dim_local[d],
                      0 );
    offset_global += a_dim_global[d];
    offset_local  += a_dim_local [d] + 2*a_ghosts;
  }

  if (!mpi->rank) {

    char filename[_MAX_STRING_SIZE_] = "";
    strcat(filename,a_fname_root);
    if (!strcmp(solver->op_overwrite,"no")) {
      strcat(filename,"_");
      strcat(filename,solver->filename_index);
    }
    strcat(filename,solver->plotfilename_extn);

#ifdef with_python
#ifdef with_python_numpy
    import_array();
    PyObject* py_plt_func = (PyObject*) solver->py_plt_func;
    PyObject* py_plt_func_args = (PyObject*) solver->py_plt_func_args;
    py_plt_func_args = PyTuple_New(7);
    {
      PyObject* py_obj = Py_BuildValue("i", a_ndims);
      PyTuple_SetItem(py_plt_func_args, 0, py_obj);
    }
    {
      PyObject* py_obj = Py_BuildValue("i", a_nvars);
      PyTuple_SetItem(py_plt_func_args, 1, py_obj);
    }
    {
      npy_intp shape[1] = {a_ndims};
      PyObject* size_arr = PyArray_SimpleNewFromData(1,shape,NPY_INT,a_dim_global);
      PyTuple_SetItem(py_plt_func_args, 2, size_arr);
    }
    {
      PyObject* py_obj = Py_BuildValue("d", a_time);
      PyTuple_SetItem(py_plt_func_args, 3, py_obj);
    }
    {
      npy_intp shape[1] = {size_global_x};
      PyObject* x_arr = PyArray_SimpleNewFromData(1,shape,NPY_DOUBLE,xg);
      PyTuple_SetItem(py_plt_func_args, 4, x_arr);
    }
    {
      npy_intp shape[1] = {size_global_u*a_nvars};
      PyObject* u_arr = PyArray_SimpleNewFromData(1,shape,NPY_DOUBLE,ug);
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

    /* Clean up output arrays */
    free(xg);
    free(ug);
  }

  return 0;
}
