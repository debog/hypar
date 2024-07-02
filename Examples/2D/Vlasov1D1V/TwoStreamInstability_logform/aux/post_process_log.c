/*
  Two-Stream Instability (Log form):
  The code takes a binary solution file (that contains the
  log of the distribution function as its input
  and writes out the distribution function to a Tecplot or
  text file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void IncrementFilename(char *f)
{
  if (f[7] == '9') {
    f[7] = '0';
    if (f[6] == '9') {
      f[6] = '0';
      if (f[5] == '9') {
        f[5] = '0';
        if (f[4] == '9') {
          f[4] = '0';
          if (f[3] == '9') {
            f[3] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[3]++;
          }
        } else {
          f[4]++;
        }
      } else {
        f[5]++;
      }
    } else {
      f[6]++;
    }
  } else {
    f[7]++;
  }
}

void WriteTecplot2D(int nvars,int imax, int jmax,double *x,double *u,char *f)
{
  printf("\tWriting tecplot solution file %s.\n",f);
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return;
  }

  double *X = x;
  double *Y = x+imax;

  /* writing tecplot data file headers */
  fprintf(out,"VARIABLES=\"I\",\"J\",\"X\",\"Y\",\"F\",\n");
  fprintf(out,"ZONE I=%d,J=%d,F=POINT\n",imax,jmax);

  /* writing the data */
  int i,j;
  for (j=0; j<jmax; j++) {
    for (i=0; i<imax; i++) {
      int v, p = i + imax*j;
      fprintf(out,"%4d %4d ",i,j);
      fprintf(out,"%1.16E %1.16E ",X[i],Y[j]);
      for (v=0; v<nvars; v++) fprintf(out,"%1.16E ",u[nvars*p+v]);
      fprintf(out,"\n");
    }
  }
  fclose(out);
  return;
}

void WriteText2D(int nvars,int imax, int jmax,double *x,double *u,char *f)
{
  printf("\tWriting text solution file %s.\n",f);
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return;
  }

  double *X = x;
  double *Y = x+imax;

  /* writing the data */
  int i,j;
  for (j=0; j<jmax; j++) {
    for (i=0; i<imax; i++) {
      int v, p = i + imax*j;
      fprintf(out,"%4d %4d ",i,j);
      fprintf(out,"%1.16E %1.16E ",X[i],Y[j]);
      for (v=0; v<nvars; v++) fprintf(out,"%1.16E ",u[nvars*p+v]);
      fprintf(out,"\n");
    }
  }
  fclose(out);
  return;
}

int PostProcess(char *fname, char *oname, int flag)
{
  FILE *in; in = fopen(fname,"rb");

  if (!in) return(-1);

  printf("Reading file %s.\n",fname);
  int ndims, nvars;
  double *U,*x;

  /* read the file headers */
  fread(&ndims,sizeof(int),1,in);
  fread(&nvars,sizeof(int),1,in);

  /* some checks */
  if (ndims != 2) {
    printf("Error: ndims in %s not equal to 2!\n",fname);
    return(1);
  }
  if (nvars != 1) {
    printf("Error: nvars in %s not equal to 4!\n",fname);
    return(1);
  }

  /* read dimensions */
  int dims[ndims];
  fread(dims,sizeof(int),ndims,in);
  printf("Dimensions: %d x %d\n",dims[0],dims[1]);
  printf("Nvars     : %d\n",nvars);

  /* allocate grid and solution arrays */
  x = (double*) calloc (dims[0]+dims[1]       ,sizeof(double));
  U = (double*) calloc (dims[0]*dims[1]*nvars ,sizeof(double));

  /* read grid and solution */
  fread(x,sizeof(double),dims[0]+dims[1]      ,in);
  fread(U,sizeof(double),dims[0]*dims[1]*nvars,in);
  /* done reading */
  fclose(in);

  int imax = dims[0];
  int jmax = dims[1];

  /* allocate distribution function array */
  double *Q = (double*) calloc (imax*jmax,sizeof(double));

  /* calculate primitive variables */
  int i, j;
  double *X           = x;
  double *Y           = x+imax;

  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      int p = i + imax*j;
      Q[p] = exp(U[p]);
    }
  }

  /* write Tecplot/Text file */
  if (flag) WriteTecplot2D(1,imax,jmax,x,Q,oname);
  else      WriteText2D   (1,imax,jmax,x,Q,oname);

  /* clean up */
  free(U);
  free(Q);
  free(x);
}

int main()
{
  FILE *out1, *out2, *in, *inputs;
  char filename[50], op_file_format[50], tecfile[50], overwrite[50];

  int flag;
  printf("Write tecplot file (1) or plain text file (0): ");
  scanf("%d",&flag);

  if ((flag != 1) && (flag != 0)) {
    printf("Error: Invalid input. Should be 1 or 0.\n");
    return(0);
  }

  printf("Reading solver.inp.\n");
  inputs = fopen("solver.inp","r");
  if (!inputs) {
    fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
    return(1);
  } else {
    char word[100];
    fscanf(inputs,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(inputs,"%s",word);
         if      (!strcmp(word, "op_file_format"   ))  fscanf(inputs,"%s" ,op_file_format);
         else if (!strcmp(word, "op_overwrite"     ))  fscanf(inputs,"%s" ,overwrite      );
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  int use_log_form = 0; /* default in HyPar is 0 */
  printf("Reading physics.inp.\n");
  inputs = fopen("physics.inp","r");
  if (!inputs) {
    fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
    return(1);
  } else {
    char word[100];
    fscanf(inputs,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(inputs,"%s",word);
        if      (!strcmp(word, "use_log_form"))    fscanf(inputs,"%d",&use_log_form);
      }
    }
    fclose(inputs);
  }
  if (!use_log_form) {
    printf("Error: \"use_log_form\" is 0 or unspecified in physics.inp.\n");
    return(0);
  }

  if (!strcmp(overwrite,"no")) {
    strcpy(filename,"op_00000.bin");
    while(1) {
      /* set filename */
      strcpy(tecfile,filename);
      tecfile[9]  = 'd';
      tecfile[10] = 'a';
      tecfile[11] = 't';
      int err = PostProcess(filename, tecfile, flag);
      if (err == -1) {
        printf("No more files found. Exiting.\n");
        break;
      }
      IncrementFilename(filename);
    }
  } else if (!strcmp(overwrite,"yes")) {
    strcpy(filename,"op.bin");
    /* set filename */
    strcpy(tecfile,filename);
    tecfile[3] = 'd';
    tecfile[4] = 'a';
    tecfile[5] = 't';
    int err = PostProcess(filename, tecfile, flag);
    if (err == -1) {
      printf("Error: op.bin not found.\n");
      return(0);
    }
  }

  return(0);
}
