/*
 * This code implements a diff between two HyPar binary
 * solutions.
 *
 * If the two solutions are identical within the specified
 * tolerances, there will be no output.
 *
 * Else, the differences will be printed to the screen.
 *
 * Please compile it using the C99 standard.
 *   gcc -std=c99 hyparDiff.c -lm -o HyParDiff
 *
 * Usage:
 *   ./HyParDiff [OPTIONS] FILENAME1 FILENAME2
 *
 * For help/usage:
 *   ./HyParDiff --help
 *   ./HyParDiff --usage
*/

#include <stdio.h>
#include <argp.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define _MAX_STRING_SIZE_ 500

#define _ArraySetValue_(x,size,value)\
  {\
    int arraycounter;\
    for (arraycounter = 0; arraycounter < (size); arraycounter++)  x[arraycounter] = (value);\
  }

#define _ArrayCopy1D_(x,y,size) \
  { \
    int arraycounter; \
    for (arraycounter = 0; arraycounter < size; arraycounter++) y[arraycounter] = x[arraycounter];\
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

#define _ArrayAXPY_(x,a,y,size)\
  { \
    int arraycounter; \
    for (arraycounter=0; arraycounter<size; arraycounter++) y[arraycounter] += a*x[arraycounter];\
  }

static char doc[] =
  "HyParDiff -- diff program for HyPar binary solution files";

const char *argp_program_version = "HyParDiff 0.0";
const char *argp_program_bug_address = "HyPar authors";

static char args_doc[] = "FILENAME1 FILENAME2";

struct arguments
{
  char* args[2];
  char* abs_tol;
  char* rel_tol;
};

static struct argp_option options[] =
{
  {"atol",'a', "ABS_TOL", 0, "Absolute tolerance"},
  {"rtol",'r', "REL_TOL", 0, "Relative tolerance"},
  {0}
};

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
   * know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
  {
    case 'a':
      arguments->abs_tol = arg;
      break;
    case 'r':
      arguments->rel_tol = arg;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 2) argp_usage (state);
      arguments->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 2) argp_usage (state);
      break;
    default:
      return ARGP_ERR_UNKNOWN;
  }

  return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

double ArrayMaxnD(int    nvars,
                  int    ndims,
                  int    *dim,
                  int    ghosts,
                  int    *index,
                  double *x )
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v;
    for (v=0; v<nvars; v++) {
      double term = ( x[p*nvars+v]>0 ? x[p*nvars+v] : -x[p*nvars+v] );
      if (term > sum) sum = term;
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

double ArraySumAbsnD(int     nvars,
                     int     ndims,
                     int     *dim,
                     int     ghosts,
                     int     *index,
                     double  *x )
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) sum += ( x[p*nvars+v]>0 ? x[p*nvars+v] : -x[p*nvars+v] );
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

double ArraySumSquarenD(int    nvars,
                        int    ndims,
                        int    *dim,
                        int    ghosts,
                        int    *index,
                        double *x )
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) sum += (x[p*nvars+v]*x[p*nvars+v]);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

typedef struct solution_object {
  int     ndims;
  int     nvars;
  int*    size;
  double* x;
  double* u;
  int*    is_periodic;
} SolObj;

void deleteSolObj(SolObj* const a_solution)
{
  if (a_solution->size) free(a_solution->size);
  if (a_solution->x) free(a_solution->x);
  if (a_solution->u) free(a_solution->u);
  if (a_solution->is_periodic) free(a_solution->is_periodic);
}

void createSolObj(const int         a_ndims,
                  const int         a_nvars,
                  const int* const  a_size,
                  SolObj* const     a_solution )
{
  a_solution->ndims = a_ndims;
  a_solution->nvars = a_nvars;

  a_solution->size = (int*) calloc(a_ndims, sizeof(int));
  _ArrayCopy1D_(a_size, a_solution->size, a_ndims);

  long sizex = 0;
  long sizeu = a_solution->nvars;
  for (int n=0; n<a_solution->ndims; n++) {
    sizex += a_solution->size[n];
    sizeu *= a_solution->size[n];
  }
  a_solution->x = (double*) calloc (sizex, sizeof(double));
  a_solution->u = (double*) calloc (sizeu, sizeof(double));

  a_solution->is_periodic = (int*) calloc (a_ndims, sizeof(int));
  _ArraySetValue_(a_solution->is_periodic, a_ndims, 1);
  return;
}

int readBinaryFile( const char* const a_filename,
                    SolObj* const     a_solution )
{
  FILE* in;
  in = fopen(a_filename,"rb");

  if (!in) {
    fprintf(stderr,"ERROR: unable to read %s.\n",a_filename);
    return 1;
  }

  fread(&a_solution->ndims, sizeof(int), 1, in);
  fread(&a_solution->nvars, sizeof(int), 1, in);

  a_solution->size = (int*) calloc(a_solution->ndims, sizeof(int));
  fread(a_solution->size, sizeof(int), a_solution->ndims, in);

  long sizex = 0;
  long sizeu = a_solution->nvars;
  for (int n=0; n<a_solution->ndims; n++) {
    sizex += a_solution->size[n];
    sizeu *= a_solution->size[n];
  }

  a_solution->x = (double*) calloc (sizex, sizeof(double));
  a_solution->u = (double*) calloc (sizeu, sizeof(double));
  fread(a_solution->x, sizeof(double), sizex, in);
  fread(a_solution->u, sizeof(double), sizeu, in);

  a_solution->is_periodic = (int*) calloc(a_solution->ndims,sizeof(int));
  _ArraySetValue_(a_solution->is_periodic, a_solution->ndims, 1);

  fclose(in);
  return 0;
}

int checkConsistency( const SolObj* const a_s1,
                      const SolObj* const a_s2 )
{
  if (a_s1->ndims != a_s2->ndims) {
    fprintf(stderr,"Error in checkConsistency(): ndims is not the same.\n");
    return 1;
  }
  if (a_s1->nvars != a_s2->nvars) {
    fprintf(stderr,"Error in checkConsistency(): nvars is not the same.\n");
    return 1;
  }
  return 0;
}

int ComputeDiffNorm(  SolObj* const a_s1,
                      SolObj* const a_s2,
                      double* const a_diff_norm_abs,
                      double* const a_diff_norm_rel )
{
  static const double tolerance = 1e-14;
  int ierr;
  int index[a_s2->ndims];

  ierr = checkConsistency(a_s1, a_s2);
  if (ierr) return ierr;

  long npoints = 1;
  long size = a_s2->nvars;
  for (int n=0; n<a_s2->ndims; n++) {
    npoints *= (a_s2->size[n]);
    size *= a_s2->size[n];
  }

  /* calculate solution norms (for relative error) */
  double solution_norm[3] = {0.0,0.0,0.0};

  /* L1 */
  {
    double sum = ArraySumAbsnD ( a_s2->nvars,
                          a_s2->ndims,
                          a_s2->size,
                          0,
                          index,
                          a_s2->u );
    solution_norm[0] = sum / ((double)npoints);
  }
  /* L2 */
  {
    double sum = ArraySumSquarenD( a_s2->nvars,
                            a_s2->ndims,
                            a_s2->size,
                            0,
                            index,
                            a_s2->u );
    solution_norm[1] = sqrt(sum/((double)npoints));
  }
  /* Linf */
  {
    double sum = ArrayMaxnD( a_s2->nvars,
                      a_s2->ndims,
                      a_s2->size,
                      0,
                      index,
                      a_s2->u );
    solution_norm[2] = sum;
  }

  _ArrayAXPY_(a_s1->u,-1.0,a_s2->u,size);

  /* L1 */
  {
    double sum = ArraySumAbsnD ( a_s2->nvars,
                          a_s2->ndims,
                          a_s2->size,
                          0,
                          index,
                          a_s2->u );
    a_diff_norm_abs[0] = sum / ((double)npoints);
  }
  /* L2 */
  {
    double sum = ArraySumSquarenD( a_s2->nvars,
                            a_s2->ndims,
                            a_s2->size,
                            0,
                            index,
                            a_s2->u );
    a_diff_norm_abs[1] = sqrt(sum/((double)npoints));
  }
  /* Linf */
  {
    double sum = ArrayMaxnD( a_s2->nvars,
                      a_s2->ndims,
                      a_s2->size,
                      0,
                      index,
                      a_s2->u );
    a_diff_norm_abs[2] = sum;
  }

  if (    (solution_norm[0] > tolerance)
      &&  (solution_norm[1] > tolerance)
      &&  (solution_norm[2] > tolerance) ) {
    a_diff_norm_rel[0] = a_diff_norm_abs[0] / solution_norm[0];
    a_diff_norm_rel[1] = a_diff_norm_abs[1] / solution_norm[1];
    a_diff_norm_rel[2] = a_diff_norm_abs[2] / solution_norm[2];
  } else {
    a_diff_norm_rel[0] = a_diff_norm_abs[0];
    a_diff_norm_rel[1] = a_diff_norm_abs[1];
    a_diff_norm_rel[2] = a_diff_norm_abs[2];
  }

  return 0;
}

int sameSize( const SolObj* const a_s1,
              const SolObj* const a_s2 )
{
  if (checkConsistency( a_s1, a_s2 ) == 0 ) {
    int retval = 1, i;
    for (i = 0; i < a_s1->ndims; i++) {
      if (a_s1->size[i] != a_s2->size[i]) retval = 0;
    }
    return retval;
  } else {
    return 0;
  }
}

int main(int argc, char** argv)
{
  int ierr;

  struct arguments arguments;

  arguments.abs_tol = "1e-14";
  arguments.rel_tol = "1e-14";

  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  char sol_fname[_MAX_STRING_SIZE_];
  char ref_fname[_MAX_STRING_SIZE_];

  strcpy(sol_fname, arguments.args[0]);
  strcpy(ref_fname, arguments.args[1]);

  char* eptr;
  double abs_tol = strtod(arguments.abs_tol, &eptr);
  double rel_tol = strtod(arguments.rel_tol, &eptr);

  SolObj solution;
  ierr = readBinaryFile(sol_fname, &solution);
  if (ierr) return ierr;

  SolObj reference;
  ierr = readBinaryFile(ref_fname, &reference);
  if (ierr) return ierr;

  if (sameSize(&solution, &reference)) {

    double diff_norm_abs[3];
    double diff_norm_rel[3];
    ierr = ComputeDiffNorm( &solution,
                            &reference,
                            diff_norm_abs,
                            diff_norm_rel );
    if (ierr) return ierr;

    int flag_ok = 1;
    if (diff_norm_abs[0] > abs_tol) flag_ok = 0;
    if (diff_norm_abs[1] > abs_tol) flag_ok = 0;
    if (diff_norm_abs[2] > abs_tol) flag_ok = 0;
    if (diff_norm_rel[0] > rel_tol) flag_ok = 0;
    if (diff_norm_rel[1] > rel_tol) flag_ok = 0;
    if (diff_norm_rel[2] > rel_tol) flag_ok = 0;

    if (!flag_ok) {
    printf( "%s and %s differ.\n",
            sol_fname, ref_fname );
      printf("Norms of the difference between the two solutions:\n");
      printf("  L1  : %1.6e (abs), %1.6e (rel)\n", diff_norm_abs[0], diff_norm_rel[0]);
      printf("  L2  : %1.6e (abs), %1.6e (rel)\n", diff_norm_abs[1], diff_norm_rel[1]);
      printf("  Linf: %1.6e (abs), %1.6e (rel)\n", diff_norm_abs[2], diff_norm_rel[2]);

    }

  } else {

    printf( "%s and %s differ in size.\n",
            sol_fname, ref_fname );

  }


  deleteSolObj(&solution);
  deleteSolObj(&reference);

  return 0;
}
