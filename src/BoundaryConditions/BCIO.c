/*! @file BCIO.c
    @author Debojyoti Ghosh
    @brief I/O functions for boundary conditions that need to read in data.
*/

#include <stdlib.h>
#include <stdio.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <mpivars.h>

/*! Read in the turbulent inflow data: The turbulent inflow data needs to be provided
    as a binary file. For parallel runs, only rank 0 reads the file, and then
    distributes the data to the other processors.
    \n\n
    This function needs to be better documented.
*/
int BCReadTurbulentInflowData(void *b,void *m,int ndims,int nvars,int *DomainSize)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;

  char    *filename     = boundary->UnsteadyDirichletFilename;
  int     *inflow_size  = NULL;
  double  *inflow_data  = NULL;
  double  *buffer       = NULL;

  int dim = boundary->dim;
  int face= boundary->face;
  int d;

  if (!mpi->rank) {

    printf("Reading turbulent inflow boundary data from %s.\n",filename);

    FILE *in;
    int  ferr;

    /* calculate the number of processors that sit on unsteady boundary */
    int nproc = 1;
    for (d=0; d<ndims; d++) nproc *= mpi->iproc[d]; nproc /= mpi->iproc[dim];

    in = fopen(filename,"rb");
    if (!in) {
      fprintf(stderr,"Error in BCReadTurbulentInflowData(): cannot open unsteady boundary data file %s.\n",filename);
      return(1);
    }
    int count = 0;
    while ((!feof(in)) && (count < nproc)) {
      int rank[ndims], size[ndims];
      ferr = fread(rank,sizeof(int),ndims,in);
      if (ferr != ndims) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (1) in file reading, count %d.\n",count);
        return(1);
      }
      if (rank[dim] != (face > 0 ? 0 : mpi->iproc[dim]-1) ) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (2) in file reading, count %d.\n",count);
        return(1);
      }
      ferr = fread(size,sizeof(int),ndims,in);
      if (ferr != ndims) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (3) in file reading, count %d.\n",count);
        return(1);
      }
      int flag = 1;
      for (d=0; d<ndims; d++) if ((d != dim) && (size[d] != DomainSize[d])) flag = 0;
      if (!flag) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (4) (dimension mismatch) in file reading, count %d.\n",count);
        return(1);
      }

      int data_size = nvars;
      for (d=0; d<ndims; d++) data_size *= size[d];
      buffer = (double*) calloc (data_size,sizeof(double));
      ferr = fread(buffer,sizeof(double),data_size,in);
      if (ferr != data_size) {
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): Error (6) in file reading, count %d.\n",count);
        return(1);
      }

      int rank1D = MPIRank1D(ndims,mpi->iproc,rank);

      if (!rank1D) {

        int index[ndims];
        inflow_size = (int*) calloc (ndims, sizeof(int));
        _ArrayCopy1D_(size,inflow_size,ndims);
        inflow_data = (double*) calloc (data_size, sizeof(double));
        ArrayCopynD(ndims,buffer,inflow_data,size,0,0,index,nvars);

      } else {
#ifndef serial
        MPI_Request req[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
        MPI_Isend(size,ndims,MPI_INT,rank1D,2152,mpi->world,&req[0]);
        MPI_Isend(buffer,data_size,MPI_DOUBLE,rank1D,2153,mpi->world,&req[1]);
        MPI_Status status_arr[3];
        MPI_Waitall(2,&req[0],status_arr);
#else
        fprintf(stderr,"Error in BCReadTurbulentInflowData(): This is a serial run. Invalid (non-zero) rank read.\n");
#endif
      }

      free(buffer);
      count++;
    }

    if (count < nproc) {
      fprintf(stderr,"Error in BCReadTurbulentInflowData(): missing data in unsteady boundary data file %s.\n",filename);
      fprintf(stderr,"Error in BCReadTurbulentInflowData(): should contain data for %d processors, ", nproc);
      fprintf(stderr,"Error in BCReadTurbulentInflowData(): but contains data for %d processors!\n", count);
      return(1);
    }

    fclose(in);

  } else {
#ifndef serial
    if (mpi->ip[dim] == (face > 0 ? 0 : mpi->iproc[dim]-1) ) {
      MPI_Request req = MPI_REQUEST_NULL;
      inflow_size = (int*) calloc (ndims,sizeof(int));
      MPI_Irecv(inflow_size,ndims,MPI_INT,0,2152,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
      int data_size = nvars;
      for (d=0; d<ndims; d++) data_size *= inflow_size[d];
      inflow_data = (double*) calloc (data_size,sizeof(double));
      MPI_Irecv(inflow_data,data_size,MPI_DOUBLE,0,2153,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
    }
#else
    fprintf(stderr,"Error in BCReadTurbulentInflowData(): Serial code should not be here!.\n");
#endif
  }

  boundary->UnsteadyDirichletSize = inflow_size;
  boundary->UnsteadyDirichletData = inflow_data;

  return(0);
}

/*! Read in the temperature data: The temperature data needs to be provided
    as a binary file. For parallel runs, only rank 0 reads the file, and then
    distributes the data to the other processors.
    \n\n
    This function needs to be better documented.
*/
int BCReadTemperatureData(void *b,void *m,int ndims,int nvars,int *DomainSize)
{
  DomainBoundary *boundary = (DomainBoundary*) b;
  MPIVariables   *mpi      = (MPIVariables*)   m;

  char    *filename = boundary->UnsteadyTemperatureFilename;
  int     *temperature_field_size  = NULL;
  double  *time_level_data         = NULL;
  double  *temperature_field_data  = NULL;
  double  *time_buffer = NULL;
  double  *data_buffer = NULL;

  int dim = boundary->dim;
  int face= boundary->face;
  int d;

  if (!mpi->rank) {

    printf("Reading boundary temperature data from %s.\n",filename);

    FILE *in;
    int  ferr;

    /* calculate the number of processors that sit on this temperature boundary */
    int nproc = 1;
    for (d=0; d<ndims; d++) nproc *= mpi->iproc[d]; nproc /= mpi->iproc[dim];

    in = fopen(filename,"rb");
    if (!in) {
      fprintf(stderr,"Error in BCReadTemperatureData(): cannot open boundary temperature data file %s.\n",filename);
      return(1);
    }

    int count = 0;
    while ((!feof(in)) && (count < nproc)) {

      int rank[ndims], size[ndims];
      ferr = fread(rank,sizeof(int),ndims,in);
      if (ferr != ndims) {
        fprintf(stderr,"Error in BCReadTemperatureData(): Error (1) in file reading, count %d.\n",count);
        return(1);
      }
      if (rank[dim] != (face > 0 ? 0 : mpi->iproc[dim]-1) ) {
        fprintf(stderr,"Error in BCReadTemperatureData(): Error (2) in file reading, count %d.\n",count);
        return(1);
      }
      ferr = fread(size,sizeof(int),ndims,in);
      if (ferr != ndims) {
        fprintf(stderr,"Error in BCReadTemperatureData(): Error (3) in file reading, count %d.\n",count);
        return(1);
      }

      int n_data = size[dim]; /* number of time levels for which temperature data is provided */
      time_buffer = (double*) calloc (n_data,sizeof(double));
      ferr = fread(time_buffer,sizeof(double),n_data,in);
      if (ferr != n_data) {
        fprintf(stderr,"Error in BCReadTemperatureData(): Error (5) in file reading, count %d.\n",count);
        return(1);
      }

      int data_size = 1;
      for (d=0; d<ndims; d++) data_size *= size[d];
      data_buffer = (double*) calloc (data_size,sizeof(double));
      ferr = fread(data_buffer,sizeof(double),data_size,in);
      if (ferr != data_size) {
        fprintf(stderr,"Error in BCReadTemperatureData(): Error (6) in file reading, count %d.\n",count);
        return(1);
      }

      int rank1D = MPIRank1D(ndims,mpi->iproc,rank);

      if (!rank1D) {

        int index[ndims];

        temperature_field_size = (int*) calloc (ndims, sizeof(int));
        _ArrayCopy1D_(size,temperature_field_size,ndims);

        time_level_data = (double*) calloc (size[dim], sizeof(double));
        _ArrayCopy1D_(time_buffer,time_level_data,size[dim]);

        temperature_field_data = (double*) calloc (data_size, sizeof(double));
        ArrayCopynD(ndims,data_buffer,temperature_field_data,size,0,0,index,1);

      } else {

#ifndef serial
        MPI_Request req[3] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL,MPI_REQUEST_NULL};
        MPI_Isend(size,ndims,MPI_INT,rank1D,2152,mpi->world,&req[0]);
        MPI_Isend(time_buffer,size[dim],MPI_DOUBLE,rank1D,2154,mpi->world,&req[2]);
        MPI_Isend(data_buffer,data_size,MPI_DOUBLE,rank1D,2153,mpi->world,&req[1]);
        MPI_Status status_arr[3];
        MPI_Waitall(3,&req[0],status_arr);
#else
        fprintf(stderr,"Error in BCReadTemperatureData(): This is a serial run. Invalid (non-zero) rank read.\n");
#endif

      }

      free(time_buffer);
      free(data_buffer);
      count++;
    }

    if (count < nproc) {
      fprintf(stderr,"Error in BCReadTemperatureData(): missing data in unsteady boundary data file %s.\n",filename);
      fprintf(stderr,"Error in BCReadTemperatureData(): should contain data for %d processors, ", nproc);
      fprintf(stderr,"Error in BCReadTemperatureData(): but contains data for %d processors!\n", count);
      return(1);
    }

    fclose(in);

  } else {

#ifndef serial
    if (mpi->ip[dim] == (face > 0 ? 0 : mpi->iproc[dim]-1) ) {

      MPI_Request req = MPI_REQUEST_NULL;

      temperature_field_size = (int*) calloc (ndims,sizeof(int));
      MPI_Irecv(temperature_field_size,ndims,MPI_INT,0,2152,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);

      int flag = 1;
      for (d=0; d<ndims; d++) if ((d != dim) && (temperature_field_size[d] != DomainSize[d])) flag = 0;
      if (!flag) {
        fprintf(stderr,"Error in BCReadTemperatureData(): Error (4) (dimension mismatch) in file reading, rank %d.\n",mpi->rank);
        return(1);
      }

      time_level_data = (double*) calloc (temperature_field_size[dim],sizeof(double));
      MPI_Irecv(time_level_data, temperature_field_size[dim], MPI_DOUBLE,0,2154,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);

      int data_size = 1;
      for (d=0; d<ndims; d++) data_size *= temperature_field_size[d];
      temperature_field_data = (double*) calloc (data_size,sizeof(double));
      MPI_Irecv(temperature_field_data,data_size,MPI_DOUBLE,0,2153,mpi->world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);

    }
#else
    fprintf(stderr,"Error in BCReadTemperatureData(): Serial code should not be here!.\n");
#endif

  }

  boundary->UnsteadyTemperatureSize = temperature_field_size;
  boundary->UnsteadyTimeLevels      = time_level_data;
  boundary->UnsteadyTemperatureData = temperature_field_data;

  return(0);
}
