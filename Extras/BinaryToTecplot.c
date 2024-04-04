/*
 * This code converts binary output from HyPar to
 * 2D/3D Tecplot format files.
 *
 * If the user provides an input file "bin2tec.inp"
 * that contains the filename (without the extension
 * ".bin") of the file they want converted, this
 * code will convert that file.
 *
 * If this input file doesn't exist, then it will
 * try to find and convert HyPar solution files.
 * For this to happen, it needs to be run at the
 * same location where the input files for the
 * simulation exist. It needs to read in stuff from
 * solver.inp and check the existence of
 * sparse_grids.inp.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define _MAX_STRING_SIZE_ 500

#define _ArraySetValue_(x,size,value) \
  { \
    int arraycounter; \
    for (arraycounter = 0; arraycounter < (size); arraycounter++)  x[arraycounter] = (value);\
  }

#define _ArrayIndex1D_(N,imax,i,ghost,index)  \
  { \
    index = i[N-1]+(ghost); \
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost))); \
    } \
  }

#define _ArrayIncrementIndex_(N,imax,i,done) \
  { \
    int arraycounter = 0; \
    while (arraycounter < (N)) { \
      if (i[arraycounter] == imax[arraycounter]-1) { \
        i[arraycounter] = 0; \
        arraycounter++; \
      } else { \
        i[arraycounter]++; \
        break; \
      } \
    } \
    if (arraycounter == (N)) done = 1; \
    else          done = 0; \
  }

void IncrementIndex(char *f)
{
  if (f[4] == '9') {
    f[4] = '0';
    if (f[3] == '9') {
      f[3] = '0';
      if (f[2] == '9') {
        f[2] = '0';
        if (f[1] == '9') {
          f[1] = '0';
          if (f[0] == '9') {
            f[0] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[0]++;
          }
        } else {
          f[1]++;
        }
      } else {
        f[2]++;
      }
    } else {
      f[3]++;
    }
  } else {
    f[4]++;
  }
}

int ReadBinary( const char* const a_fname,
                int* const        a_ndims,
                int* const        a_nvars,
                int** const       a_dims,
                double** const    a_x,
                double** const    a_u )
{
  FILE* in = fopen(a_fname,"rb");
  if (!in) return 1;

  printf("Reading file %s.\n",a_fname);
  size_t ferr;

  /* read the file headers */
  int ndims, nvars;
  ferr = fread(&ndims,sizeof(int),1,in);
  if (ferr != 1) {
    fprintf(stderr,
            "Error while reading %s: unable to read ndims.\n",
            a_fname);
    return 1;
  }
  ferr = fread(&nvars,sizeof(int),1,in);
  if (ferr != 1) {
    fprintf(stderr,
            "Error while reading %s: unable to read nvars.\n",
            a_fname);
    return 1;
  }
  *a_ndims = ndims;
  *a_nvars = nvars;
  printf("  ndims: %d\n",ndims);
  printf("  nvars: %d\n",nvars);

  /* some checks */
  if ((ndims != 2) && (ndims != 3)) {
    printf("Error: ndims in %s not equal to 2 or 3!\n", a_fname);
    return 1;
  }

  /* read dimensions */
  int* dims = calloc(ndims, sizeof(int));
  ferr = fread(dims,sizeof(int),ndims,in);
  if (ferr != ndims) {
    fprintf(stderr,
            "Error while reading %s: unable to read dims.\n",
            a_fname);
    return 1;
  }
  *a_dims = dims;
  if      (ndims == 2) printf("  dimensions: %d x %d\n",dims[0],dims[1]);
  else if (ndims == 3) printf("  dimensions: %d x %d x %d\n",dims[0],dims[1],dims[2]);

  double *U,*x;
  /* allocate grid and solution arrays */
  int d;
  int sizex = 0;      for (d=0; d<ndims; d++) sizex += dims[d];
  int sizeu = nvars;  for (d=0; d<ndims; d++) sizeu *= dims[d];
  x = (double*) calloc (sizex,sizeof(double));
  U = (double*) calloc (sizeu,sizeof(double));

  /* read grid and solution */
  ferr = fread(x,sizeof(double),sizex,in);
  if (ferr != sizex) {
    fprintf(stderr,
            "Error while reading %s: unable to read x.\n",
            a_fname);
    return 1;
  }
  ferr = fread(U,sizeof(double),sizeu,in);
  if (ferr != sizeu) {
    fprintf(stderr,
            "Error while reading %s: unable to read u.\n",
            a_fname);
    return 1;
  }
  /* done reading */
  fclose(in);

  *a_x = x;
  *a_u = U;

  return 0;
}

void WriteTecplot2D(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  if (ndims !=2) {
    fprintf(stderr,"Error in WriteTecplot2D(): This functions is hardcoded for 2-dimensional ");
    fprintf(stderr,"problems only. Instead, ndims=%d.\n",ndims);
    return;
  }
  int i;
  int imax = dim[0];
  int jmax = dim[1];

  printf("  Writing tecplot solution file %s.\n",f);
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return;
  }

  /* writing tecplot data file headers */
  fprintf(out,"VARIABLES=\"I\",\"J\",\"X\",\"Y\",");
  char varname[3] = "00";
  for (i = 0; i < nvars; i++) {
    fprintf(out,"\"%s\",",varname);
    if (varname[1] == '9') { varname[0]++; varname[1] = '0'; }
    else                     varname[1]++;
  }
  fprintf(out,"\n");
  fprintf(out,"ZONE I=%d,J=%d,F=POINT\n",imax,jmax);

  /* writing the data */
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int i, p;
    _ArrayIndex1D_(ndims,dim,index,0,p);
    for (i=0; i<ndims; i++) fprintf(out,"%4d ",index[i]);
    for (i=0; i<ndims; i++) {
      int j,offset = 0; for (j=0; j<i; j++) offset += dim[j];
      fprintf(out,"%+E ",x[offset+index[i]]);
    }
    for (i=0; i<nvars; i++) fprintf(out,"%+E ",u[nvars*p+i]);
    fprintf(out,"\n");
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  fclose(out);
  return;
}

int WriteTecplot3D(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  if (ndims !=3) {
    fprintf(stderr,"Error in WriteTecplot3D(): This functions is hardcoded for 3-dimensional ");
    fprintf(stderr,"problems only. Instead, ndims=%d.\n",ndims);
    return(1);
  }
  int i;
  int imax = dim[0];
  int jmax = dim[1];
  int kmax = dim[2];

  printf("  Writing tecplot solution file %s.\n",f);
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  /* writing tecplot data file headers */
  fprintf(out,"VARIABLES=\"I\",\"J\",\"K\",\"X\",\"Y\",\"Z\",");
  char varname[3] = "00";
  for (i = 0; i < nvars; i++) {
    fprintf(out,"\"%s\",",varname);
    if (varname[1] == '9') { varname[0]++; varname[1] = '0'; }
    else                     varname[1]++;
  }
  fprintf(out,"\n");
  fprintf(out,"ZONE I=%d,J=%d,K=%d,F=POINT\n",imax,jmax,kmax);

  /* writing the data */
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int i, p;
    _ArrayIndex1D_(ndims,dim,index,0,p);
    for (i=0; i<ndims; i++) fprintf(out,"%4d ",index[i]);
    for (i=0; i<ndims; i++) {
      int j,offset = 0; for (j=0; j<i; j++) offset += dim[j];
      fprintf(out,"%+E ",x[offset+index[i]]);
    }
    for (i=0; i<nvars; i++) fprintf(out,"%+E ",u[nvars*p+i]);
    fprintf(out,"\n");
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  fclose(out);
  return(0);
}

int convertFile(const char* const fname_root)
{
  char file_in[_MAX_STRING_SIZE_],
       file_out[_MAX_STRING_SIZE_];

  strcpy(file_in, fname_root);
  strcat(file_in, ".");
  strcat(file_in, "bin");

  strcpy(file_out, fname_root);
  strcat(file_out, ".");
  strcat(file_out, "dat");

  int     ndims, nvars;
  int*    dims = NULL;
  double* x = NULL;
  double* U = NULL;

  int ierr;
  ierr = ReadBinary(file_in, &ndims, &nvars, &dims, &x, &U);

  if (!ierr) {

    /* write Tecplot file */
    int ind[ndims];
    if      (ndims == 2) WriteTecplot2D(2,nvars,dims,x,U,file_out,&ind[0]);
    else if (ndims == 3) WriteTecplot3D(3,nvars,dims,x,U,file_out,&ind[0]);

    /* clean up */
    free(dims);
    free(U);
    free(x);

  }

  return ierr;
}

int main()
{
  char fname_root[_MAX_STRING_SIZE_];
  char overwrite[_MAX_STRING_SIZE_];

  char input_fname_user[_MAX_STRING_SIZE_] = "bin2tec.inp";
  char input_fname_solver[_MAX_STRING_SIZE_] = "solver.inp";
  char input_fname_sg[_MAX_STRING_SIZE_] = "sparse_grids.inp";

  FILE* inputs;
  inputs = fopen(input_fname_user,"r");
  if (inputs) {
    printf("Reading filename root from %s.\n", input_fname_user);
    fscanf(inputs, "%s", fname_root);
    printf("  filename root is %s.\n", fname_root);
    fclose(inputs);
    strcpy(overwrite,"yes");
  } else {
    printf("Reading simulation inputs from %s.\n", input_fname_solver);
    inputs = fopen(input_fname_solver,"r");
    char op_file_format[_MAX_STRING_SIZE_];
    if (!inputs) {
      fprintf(stderr,"Error: File %s not found.\n", input_fname_solver);
      return 1;
    } else {
      char word[100];
      fscanf(inputs,"%s",word);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          fscanf(inputs,"%s",word);
           if (!strcmp(word, "op_file_format"   ))  fscanf(inputs,"%s" ,op_file_format);
           else if (!strcmp(word, "op_overwrite"     ))  fscanf(inputs,"%s" ,overwrite      );
        }
      }
      fclose(inputs);
    }
    if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
      printf("Error: solution output needs to be in binary files.\n");
      return 1;
    }
    inputs = fopen(input_fname_sg, "r");
    if (inputs) {
      fclose(inputs);
      printf("  Found %s - this is a sparse grids simulation.\n", input_fname_sg);
      strcpy(fname_root, "op_fg");
    } else {
      strcpy(fname_root, "op");
    }
    printf("  filename root is %s.\n", fname_root);
  }

  if (!strcmp(overwrite,"no")) {

    char fname_idx[6] = "00000";

    while(1) {

      char filename[_MAX_STRING_SIZE_];
      strcpy(filename,fname_root);
      strcat(filename,"_");
      strcat(filename,fname_idx);

      int ierr = convertFile(filename);
      if (ierr) {
        printf("No more files found (%s). Exiting.\n", filename);
        break;
      }

      IncrementIndex(fname_idx);
    }

  } else if (!strcmp(overwrite,"yes")) {

    char filename[_MAX_STRING_SIZE_];
    strcpy(filename,fname_root);

    int ierr = convertFile(filename);
    if (ierr) {
      printf("File not found or not in correct format (%s). Exiting.\n", filename);
    }

  }

  return 0;
}
