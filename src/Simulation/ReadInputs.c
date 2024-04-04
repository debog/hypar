/*! @file ReadInputs.c
    @author Debojyoti Ghosh
    @brief Read the input parameters from \b solver.inp
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <timeintegration.h>
#include <mpivars.h>
#include <simulation_object.h>

/*! Read the simulation inputs from the file \b solver.inp.
    Rank 0 reads in the inputs and broadcasts them to all the
    processors.\n\n
    The format of \b solver.inp is as follows:\n

        begin
            <keyword>   <value>
            <keyword>   <value>
            <keyword>   <value>
            ...
            <keyword>   <value>
        end

    where the list of keywords and their type are:\n
    Keyword name       | Type         | Variable                      | Default value
    ------------------ | ------------ | ----------------------------- | ------------------------------------------
    ndims              | int          | #HyPar::ndims                 | 1
    nvars              | int          | #HyPar::nvars                 | 1
    size               | int[ndims]   | #HyPar::dim_global            | must be specified
    iproc              | int[ndims]   | #MPIVariables::iproc          | must be specified (see notes below)
    ghost              | int          | #HyPar::ghosts                | 1
    n_iter             | int          | #HyPar::n_iter                | 0
    restart_iter       | int          | #HyPar::restart_iter          | 0
    time_scheme        | char[]       | #HyPar::time_scheme           | euler
    time_scheme_type   | char[]       | #HyPar::time_scheme_type      | none
    hyp_space_scheme   | char[]       | #HyPar::spatial_scheme_hyp    | 1
    hyp_flux_split     | char[]       | #HyPar::SplitHyperbolicFlux   | no
    hyp_interp_type    | char[]       | #HyPar::interp_type           | characteristic
    par_space_type     | char[]       | #HyPar::spatial_type_par      | nonconservative-1stage
    par_space_scheme   | char[]       | #HyPar::spatial_scheme_par    | 2
    dt                 | double       | #HyPar::dt                    | 0.0
    conservation_check | char[]       | #HyPar::ConservationCheck     | no
    screen_op_iter     | int          | #HyPar::screen_op_iter        | 1
    file_op_iter       | int          | #HyPar::file_op_iter          | 1000
    op_file_format     | char[]       | #HyPar::op_file_format        | text
    ip_file_type       | char[]       | #HyPar::ip_file_type          | ascii
    input_mode         | char[]       | #HyPar::input_mode            | serial
    output_mode        | char[]       | #HyPar::output_mode           | serial
    op_overwrite       | char[]       | #HyPar::op_overwrite          | no
    plot_solution      | char[]       | #HyPar::plot_solution         | no
    model              | char[]       | #HyPar::model                 | must be specified
    immersed_body      | char[]       | #HyPar::ib_filename           | "none"
    size_exact         | int[ndims]   | #HyPar::dim_global_ex         | #HyPar::dim_global
    use_gpu            | char[]       | #HyPar::use_gpu               | no
    gpu_device_no      | int          | #HyPar::gpu_device_no         | -1

    \b Notes:
    + "ndims" \b must be specified \b before "size".
    + the input "iproc" is ignored when running a sparse grids simulation.
    + if "input_mode" or "output_mode" are set to "parallel" or "mpi-io",
      the number of I/O ranks must be specified right after as an integer.
      For example:

          begin
              ...
              input_mode  parallel 4
              ...
          end

      This means that 4 MPI ranks will participate in file I/O (assuming
      total MPI ranks is more than 4) (see ReadArrayParallel(),
      WriteArrayParallel(), ReadArrayMPI_IO() ).
      - The number of I/O ranks specified for "input_mode" and "output_mode"
        \b must \b be \b same. Otherwise, the value for the one specified last
        will be used.
      - The number of I/O ranks must be such that the total number of MPI ranks
        is an integer multiple. Otherwise, the code will use only 1 I/O rank.
    + If any of the keywords are not present, the default value is used, except
      the ones whose default values say "must be specified". Thus, keywords that
      are not required for a particular simulation may be left out of the
      solver.inp input file. For example,
      - a #Euler1D simulation does not need "par_space_type" or "par_space_scheme"
        because it does not have a parabolic term.
      - unless a conservation check is required, "conservation_check" can be left
        out and the code will not check for conservation.
      - "immersed_body" need not be specified if there are no immersed bodies present.
         \b NOTE: However, if it is specified, and a file of that filename does not
         exist, it will result in an error.
*/
int ReadInputs( void  *s,     /*!< Array of simulation objects of type #SimulationObject
                                   of size nsims */
                int   nsims,  /*!< Number of simulation objects */
                int   rank    /*!< MPI rank of this process */
              )
{
  SimulationObject *sim = (SimulationObject*) s;
  int n, ferr    = 0;

  if (sim == NULL) {
    printf("Error: simulation object array is NULL!\n");
    printf("Please consider killing this run.\n");
    return(1);
  }

  if (!rank) {

    /* set some default values for optional inputs */
    for (n = 0; n < nsims; n++) {
      sim[n].solver.ndims           = 1;
      sim[n].solver.nvars           = 1;
      sim[n].solver.ghosts          = 1;
      sim[n].solver.dim_global      = NULL;
      sim[n].solver.dim_local       = NULL;
      sim[n].solver.dim_global_ex   = NULL;
      sim[n].mpi.iproc              = NULL;
      sim[n].mpi.N_IORanks          = 1;
      sim[n].solver.dt              = 0.0;
      sim[n].solver.n_iter          = 0;
      sim[n].solver.restart_iter    = 0;
      sim[n].solver.screen_op_iter  = 1;
      sim[n].solver.file_op_iter    = 1000;
      sim[n].solver.write_residual  = 0;
      sim[n].solver.flag_ib         = 0;
#if defined(HAVE_CUDA)
      sim[n].solver.use_gpu         = 0;
      sim[n].solver.gpu_device_no   = -1;
#endif
      strcpy(sim[n].solver.time_scheme        ,"euler"         );
      strcpy(sim[n].solver.time_scheme_type   ," "             );
      strcpy(sim[n].solver.spatial_scheme_hyp ,"1"             );
      strcpy(sim[n].solver.spatial_type_par   ,_NC_1STAGE_     );
      strcpy(sim[n].solver.spatial_scheme_par ,"2"             );
      strcpy(sim[n].solver.interp_type        ,"characteristic");
      strcpy(sim[n].solver.ip_file_type       ,"ascii"         );
      strcpy(sim[n].solver.input_mode         ,"serial"        );
      strcpy(sim[n].solver.output_mode        ,"serial"        );
      strcpy(sim[n].solver.op_file_format     ,"text"          );
      strcpy(sim[n].solver.op_overwrite       ,"no"            );
      strcpy(sim[n].solver.plot_solution      ,"no"            );
      strcpy(sim[n].solver.model              ,"none"          );
      strcpy(sim[n].solver.ConservationCheck  ,"no"            );
      strcpy(sim[n].solver.SplitHyperbolicFlux,"no"            );
      strcpy(sim[n].solver.ib_filename        ,"none"          );
    }

    /* open the file */
    FILE *in;
    printf("Reading solver inputs from file \"solver.inp\".\n");
    in = fopen("solver.inp","r");
    if (!in) {
      fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
      fprintf(stderr,"Please consider killing this run.\n");
      return(1);
    }

    /* reading solver inputs */
    char word[_MAX_STRING_SIZE_];
    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);

    if (!strcmp(word, "begin")){

      while (strcmp(word, "end")) {

        ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);

        if (!strcmp(word, "ndims")) {

          ferr = fscanf(in,"%d",&(sim[0].solver.ndims)); if (ferr != 1) return(1);
          sim[0].solver.dim_global    = (int*) calloc (sim[0].solver.ndims,sizeof(int));
          sim[0].mpi.iproc            = (int*) calloc (sim[0].solver.ndims,sizeof(int));
          sim[0].solver.dim_global_ex = (int*) calloc (sim[0].solver.ndims,sizeof(int));

          int n;
          for (n = 1; n < nsims; n++) {
            sim[n].solver.ndims = sim[0].solver.ndims;
            sim[n].solver.dim_global    = (int*) calloc (sim[n].solver.ndims,sizeof(int));
            sim[n].mpi.iproc            = (int*) calloc (sim[n].solver.ndims,sizeof(int));
            sim[n].solver.dim_global_ex = (int*) calloc (sim[n].solver.ndims,sizeof(int));
          }

        } else if (!strcmp(word, "nvars")) {

          ferr = fscanf(in,"%d",&(sim[0].solver.nvars));
          for (int n = 1; n < nsims; n++) sim[n].solver.nvars = sim[0].solver.nvars;

        } else if   (!strcmp(word, "size")) {

          for (int n = 0; n < nsims; n++) {
            if (!sim[n].solver.dim_global) {
              fprintf(stderr,"Error in ReadInputs(): dim_global not allocated for n=%d.\n", n);
              fprintf(stderr,"Please specify ndims before dimensions.\n"         );
              return(1);
            } else {
              for (int i=0; i<sim[n].solver.ndims; i++) {
                ferr = fscanf(in,"%d",&(sim[n].solver.dim_global[i]));
                if (ferr != 1) {
                  fprintf(stderr,"Error in ReadInputs() while reading grid sizes for domain %d.\n", n);
                  return(1);
                }
                sim[n].solver.dim_global_ex[i] = sim[n].solver.dim_global[i];
              }
            }
          }

        } else if   (!strcmp(word, "size_exact")) {

          for (int n = 0; n < nsims; n++) {
            if (!sim[n].solver.dim_global_ex) {
              fprintf(stderr,"Error in ReadInputs(): dim_global_ex not allocated for n=%d.\n", n);
              fprintf(stderr,"Please specify ndims before dimensions.\n"         );
              return(1);
            } else {
              for (int i=0; i<sim[n].solver.ndims; i++) {
                ferr = fscanf(in,"%d",&(sim[n].solver.dim_global_ex[i]));
                if (ferr != 1) {
                  fprintf(stderr,"Error in ReadInputs() while reading exact solution grid sizes for domain %d.\n", n);
                  return(1);
                }
              }
            }
          }

        } else if (!strcmp(word, "iproc")) {

          int n;
          for (n = 0; n < nsims; n++) {
            if (!sim[n].mpi.iproc) {
              fprintf(stderr,"Error in ReadInputs(): iproc not allocated for n=%d.\n", n);
              fprintf(stderr,"Please specify ndims before iproc.\n"         );
              return(1);
            } else {
              int i;
              for (i=0; i<sim[n].solver.ndims; i++) {
                ferr = fscanf(in,"%d",&(sim[n].mpi.iproc[i]));
                if (ferr != 1) {
                  fprintf(stderr,"Error in ReadInputs() while reading iproc for domain %d.\n", n);
                  return(1);
                }
              }
            }
          }

        } else if (!strcmp(word, "ghost")) {

          ferr = fscanf(in,"%d",&(sim[0].solver.ghosts));

          int n;
          for (n = 1; n < nsims; n++) sim[n].solver.ghosts = sim[0].solver.ghosts;

        } else if (!strcmp(word, "n_iter")) {

          ferr = fscanf(in,"%d",&(sim[0].solver.n_iter));

          int n;
          for (n = 1; n < nsims; n++) sim[n].solver.n_iter = sim[0].solver.n_iter;

        } else if (!strcmp(word, "restart_iter")) {

          ferr = fscanf(in,"%d",&(sim[0].solver.restart_iter));

          int n;
          for (n = 1; n < nsims; n++) sim[n].solver.restart_iter = sim[0].solver.restart_iter;

        } else if (!strcmp(word, "time_scheme")) {

          ferr = fscanf(in,"%s",sim[0].solver.time_scheme);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.time_scheme, sim[0].solver.time_scheme);

        }  else if (!strcmp(word, "time_scheme_type" )) {

          ferr = fscanf(in,"%s",sim[0].solver.time_scheme_type);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.time_scheme_type, sim[0].solver.time_scheme_type);

        }  else if (!strcmp(word, "hyp_space_scheme")) {

          ferr = fscanf(in,"%s",sim[0].solver.spatial_scheme_hyp);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.spatial_scheme_hyp, sim[0].solver.spatial_scheme_hyp);

        }  else if (!strcmp(word, "hyp_flux_split")) {

          ferr = fscanf(in,"%s",sim[0].solver.SplitHyperbolicFlux);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.SplitHyperbolicFlux, sim[0].solver.SplitHyperbolicFlux);

        }  else if (!strcmp(word, "hyp_interp_type")) {

          ferr = fscanf(in,"%s",sim[0].solver.interp_type);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.interp_type, sim[0].solver.interp_type);

        }  else if (!strcmp(word, "par_space_type")) {

          ferr = fscanf(in,"%s",sim[0].solver.spatial_type_par);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.spatial_type_par, sim[0].solver.spatial_type_par);

        }  else if (!strcmp(word, "par_space_scheme")) {

          ferr = fscanf(in,"%s",sim[0].solver.spatial_scheme_par);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.spatial_scheme_par, sim[0].solver.spatial_scheme_par);

        }  else if (!strcmp(word, "dt")) {

          ferr = fscanf(in,"%lf",&(sim[0].solver.dt));

          int n;
          for (n = 1; n < nsims; n++) sim[n].solver.dt = sim[0].solver.dt;

        }  else if (!strcmp(word, "conservation_check" )) {

          ferr = fscanf(in,"%s",sim[0].solver.ConservationCheck);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.ConservationCheck, sim[0].solver.ConservationCheck);

        }  else if (!strcmp(word, "screen_op_iter")) {

          ferr = fscanf(in,"%d",&(sim[0].solver.screen_op_iter));

          int n;
          for (n = 1; n < nsims; n++) sim[n].solver.screen_op_iter = sim[0].solver.screen_op_iter;

        }  else if (!strcmp(word, "file_op_iter")) {

          ferr = fscanf(in,"%d",&(sim[0].solver.file_op_iter));

          int n;
          for (n = 1; n < nsims; n++) sim[n].solver.file_op_iter = sim[0].solver.file_op_iter;

        }  else if (!strcmp(word, "op_file_format")) {

          ferr = fscanf(in,"%s",sim[0].solver.op_file_format);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.op_file_format, sim[0].solver.op_file_format);

        }  else if (!strcmp(word, "ip_file_type")) {

          ferr = fscanf(in,"%s",sim[0].solver.ip_file_type);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.ip_file_type, sim[0].solver.ip_file_type);

        }  else if (!strcmp(word, "input_mode")) {

          ferr = fscanf(in,"%s",sim[0].solver.input_mode);
          if (strcmp(sim[0].solver.input_mode,"serial")) ferr = fscanf(in,"%d",&(sim[0].mpi.N_IORanks));

          int n;
          for (n = 1; n < nsims; n++) {
            strcpy(sim[n].solver.input_mode, sim[0].solver.input_mode);
            if (strcmp(sim[n].solver.input_mode,"serial")) sim[n].mpi.N_IORanks = sim[0].mpi.N_IORanks;
          }

         } else if (!strcmp(word, "output_mode"))  {

          ferr = fscanf(in,"%s",sim[0].solver.output_mode);
          if (strcmp(sim[0].solver.output_mode,"serial")) ferr = fscanf(in,"%d",&(sim[0].mpi.N_IORanks));

          int n;
          for (n = 1; n < nsims; n++) {
            strcpy(sim[n].solver.output_mode, sim[0].solver.output_mode);
            if (strcmp(sim[n].solver.output_mode,"serial")) sim[n].mpi.N_IORanks = sim[0].mpi.N_IORanks;
          }

        } else if   (!strcmp(word, "op_overwrite")) {

          ferr = fscanf(in,"%s",sim[0].solver.op_overwrite);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.op_overwrite, sim[0].solver.op_overwrite);

        } else if   (!strcmp(word, "plot_solution")) {

          ferr = fscanf(in,"%s",sim[0].solver.plot_solution);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.plot_solution, sim[0].solver.plot_solution);

        }  else if (!strcmp(word, "model")) {

          ferr = fscanf(in,"%s",sim[0].solver.model);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.model, sim[0].solver.model);

        }  else if (!strcmp(word, "immersed_body")) {

          ferr = fscanf(in,"%s",sim[0].solver.ib_filename);

          int n;
          for (n = 1; n < nsims; n++) strcpy(sim[n].solver.ib_filename, sim[0].solver.ib_filename);

        }
#if defined(HAVE_CUDA)
        else if (!strcmp(word, "use_gpu")) {
          ferr = fscanf(in,"%s",word);
          if (!strcmp(word, "yes") || !strcmp(word, "true")) sim[0].solver.use_gpu = 1;

          int n;
          for (n = 1; n < nsims; n++) sim[n].solver.use_gpu = sim[0].solver.use_gpu;
        } else if (!strcmp(word, "gpu_device_no")) {
          ferr = fscanf(in,"%d", &sim[0].solver.gpu_device_no);

          int n;
          for (n = 1; n < nsims; n++) sim[n].solver.gpu_device_no = sim[0].solver.gpu_device_no;
        }
#endif
        else if (strcmp(word, "end")) {

          char useless[_MAX_STRING_SIZE_];
          ferr = fscanf(in,"%s",useless);
          printf("Warning: keyword %s in file \"solver.inp\" with value %s not recognized or extraneous. Ignoring.\n",
                  word,useless);

        }
        if (ferr != 1) return(1);

      }

    } else {

       fprintf(stderr,"Error: Illegal format in file \"solver.inp\".\n");
      return(1);

    }

    /* close the file */
    fclose(in);

    /* some checks */
    for (n = 0; n < nsims; n++) {

      if (sim[n].solver.screen_op_iter <= 0)  sim[n].solver.screen_op_iter = 1;
      if (sim[n].solver.file_op_iter <= 0)    sim[n].solver.file_op_iter   = sim[n].solver.n_iter;

      if ((sim[n].solver.ndims != 3) && (strcmp(sim[n].solver.ib_filename,"none"))) {
        printf("Warning: immersed boundaries not implemented for ndims = %d. ",sim[n].solver.ndims);
        printf("Ignoring input for \"immersed_body\" (%s).\n",sim[n].solver.ib_filename);
        strcpy(sim[n].solver.ib_filename,"none");
      }
      sim[n].solver.flag_ib = strcmp(sim[n].solver.ib_filename,"none");

      /* restart only supported for binary output files */
      if ((sim[n].solver.restart_iter != 0) && strcmp(sim[n].solver.op_file_format,"binary")) {
        if (!sim[n].mpi.rank) fprintf(stderr,"Error in ReadInputs(): Restart is supported only for binary output files.\n");
        return(1);
      }
    }
  }

#ifndef serial
  for (n = 0; n < nsims; n++) {

    /* Broadcast the input parameters */
    MPIBroadcast_integer(&(sim[n].solver.ndims),1,0,&(sim[n].mpi.world));
    if (sim[n].mpi.rank) {
      sim[n].solver.dim_global    = (int*) calloc (sim[n].solver.ndims,sizeof(int));
      sim[n].mpi.iproc            = (int*) calloc (sim[n].solver.ndims,sizeof(int));
      sim[n].solver.dim_global_ex = (int*) calloc (sim[n].solver.ndims,sizeof(int));
    }
    MPIBroadcast_integer(&(sim[n].solver.nvars)         ,1                  ,0,&(sim[n].mpi.world));
    MPIBroadcast_integer( sim[n].solver.dim_global      ,sim[n].solver.ndims,0,&(sim[n].mpi.world));
    MPIBroadcast_integer( sim[n].solver.dim_global_ex   ,sim[n].solver.ndims,0,&(sim[n].mpi.world));
    MPIBroadcast_integer( sim[n].mpi.iproc              ,sim[n].solver.ndims,0,&(sim[n].mpi.world));
    MPIBroadcast_integer(&(sim[n].mpi.N_IORanks)        ,1                  ,0,&(sim[n].mpi.world));
    MPIBroadcast_integer(&(sim[n].solver.ghosts)        ,1                  ,0,&(sim[n].mpi.world));
    MPIBroadcast_integer(&(sim[n].solver.n_iter)        ,1                  ,0,&(sim[n].mpi.world));
    MPIBroadcast_integer(&(sim[n].solver.restart_iter)  ,1                  ,0,&(sim[n].mpi.world));
    MPIBroadcast_integer(&(sim[n].solver.screen_op_iter),1                  ,0,&(sim[n].mpi.world));
    MPIBroadcast_integer(&(sim[n].solver.file_op_iter)  ,1                  ,0,&(sim[n].mpi.world));
    MPIBroadcast_integer(&(sim[n].solver.flag_ib)       ,1                  ,0,&(sim[n].mpi.world));
#if defined(HAVE_CUDA)
    MPIBroadcast_integer(&(sim[n].solver.use_gpu)       ,1                  ,0,&(sim[n].mpi.world));
    MPIBroadcast_integer(&(sim[n].solver.gpu_device_no) ,1                  ,0,&(sim[n].mpi.world));
#endif
    MPIBroadcast_character(sim[n].solver.time_scheme        ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.time_scheme_type   ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.spatial_scheme_hyp ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.interp_type        ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.spatial_type_par   ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.spatial_scheme_par ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.ConservationCheck  ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.SplitHyperbolicFlux,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.op_file_format     ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.ip_file_type       ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.input_mode         ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.output_mode        ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.op_overwrite       ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.plot_solution      ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.model              ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));
    MPIBroadcast_character(sim[n].solver.ib_filename        ,_MAX_STRING_SIZE_,0,&(sim[n].mpi.world));

    MPIBroadcast_double(&(sim[n].solver.dt),1,0,&(sim[n].mpi.world));
  }
#endif

  return 0;
}
