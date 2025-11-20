/*! @file NavierStokes3DIBForces.c
    @brief Compute aerodynamic forces on immersed body, if present
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <immersedboundaries.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

int NavierStokes3DComputePressure(double*, const double* const, void*);
int NavierStokes3DComputeTemperature(double*, const double* const, void*);

/*!
  Calculate the shear forces on the immersed body surface: for each
  "local facet", i.e., facets (of the immersed body surface) that lie within
  the local computational domain of this MPI rank, compute the shear forces
  at the facet centroid from the flow variables at the grid points surrounding
  the facet.

  The array to hold the computed shear forces *must* be NULL when this
  function is called. If the current subdomain contains a part of the
  immersed body, they will point to arrays with the local data. Otherwise,
  they will remain NULL.

  The shear force array will be of the size (4 X #ImmersedBoundary::nfacets_local),
  where the ordering of the data is: x-component, y-component, z-component, magnitude.
*/
static int ComputeShear(void *a_s,              /*!< Solver object of type #HyPar */
                        void *a_m,              /*!< MPI object of type #MPIVariables */
                        const double* const a_u,/*!< Array containing the conserved flow variables */
                        double** const a_sf     /*!< Array for (x,y,z)-components & magnitude of shear */
                       )
{
  HyPar             *solver  = (HyPar*)          a_s;
  MPIVariables      *mpi     = (MPIVariables*)   a_m;
  NavierStokes3D    *physics = (NavierStokes3D*) solver->m_physics;
  ImmersedBoundary  *IB      = (ImmersedBoundary*) solver->m_ib;

  if ((*a_sf) != NULL) {
    fprintf(stderr, "Error in ComputeShear()\n");
    fprintf(stderr, " shear force array is not NULL!\n");
    return 1;
  }

  if (!solver->m_flag_ib) return(0);

  int           nfacets_local = IB->m_nfacets_local;
  FacetMap      *fmap = IB->m_fmap;
  static double v[_MODEL_NVARS_];

  int nv = 4;

  if (nfacets_local > 0) {

    (*a_sf) = (double*) calloc (nv*nfacets_local, sizeof(double));

    if (physics->m_Re > 0) {

      for (int n = 0; n < nfacets_local; n++) {

        double *alpha;
        int    *nodes, j, k;

        alpha = &(fmap[n].m_interp_coeffs[0]);
        nodes = &(fmap[n].m_interp_nodes[0]);
        _ArraySetValue_(v,_MODEL_NVARS_,0.0);
        for (j=0; j<_IB_NNODES_; j++) {
          for (k=0; k<_MODEL_NVARS_; k++) {
            v[k] += ( alpha[j] * a_u[_MODEL_NVARS_*nodes[j]+k] );
          }
        }
        double rho_c, uvel_c, vvel_c, wvel_c, energy_c, pressure_c;
        _NavierStokes3DGetFlowVar_(v,_NavierStokes3D_stride_,rho_c,uvel_c,vvel_c,wvel_c,energy_c,pressure_c,physics->m_gamma);

        alpha = &(fmap[n].m_interp_coeffs_ns[0]);
        nodes = &(fmap[n].m_interp_nodes_ns[0]);
        _ArraySetValue_(v,_MODEL_NVARS_,0.0);
        for (j=0; j<_IB_NNODES_; j++) {
          for (k=0; k<_MODEL_NVARS_; k++) {
            v[k] += ( alpha[j] * a_u[_MODEL_NVARS_*nodes[j]+k] );
          }
        }
        double rho_ns, uvel_ns, vvel_ns, wvel_ns, energy_ns, pressure_ns;
        _NavierStokes3DGetFlowVar_(v,_NavierStokes3D_stride_,rho_ns,uvel_ns,vvel_ns,wvel_ns,energy_ns,pressure_ns,physics->m_gamma);

        double u_x = (uvel_ns - uvel_c) / fmap[n].m_dx;
        double v_x = (vvel_ns - vvel_c) / fmap[n].m_dx;
        double w_x = (wvel_ns - wvel_c) / fmap[n].m_dx;

        double u_y = (uvel_ns - uvel_c) / fmap[n].m_dy;
        double v_y = (vvel_ns - vvel_c) / fmap[n].m_dy;
        double w_y = (wvel_ns - wvel_c) / fmap[n].m_dy;

        double u_z = (uvel_ns - uvel_c) / fmap[n].m_dz;
        double v_z = (vvel_ns - vvel_c) / fmap[n].m_dz;
        double w_z = (wvel_ns - wvel_c) / fmap[n].m_dz;

        double nx = fmap[n].m_facet->m_nx;
        double ny = fmap[n].m_facet->m_ny;
        double nz = fmap[n].m_facet->m_nz;

        double T      = physics->m_gamma*pressure_c/rho_c;
        double mu     = raiseto(T, 0.76);
        double inv_Re = 1.0/physics->m_Re;

        double tau_x = (mu*inv_Re) * (2*u_x*nx + (u_y+v_x)*ny + (u_z+w_x)*nz);
        double tau_y = (mu*inv_Re) * ((v_x+u_y)*nx + 2*v_y*ny + (v_z+w_y)*nz);
        double tau_z = (mu*inv_Re) * ((w_x+u_z)*nx + (w_y+v_z)*ny + 2*w_z*nz);

        (*a_sf)[n*nv+_XDIR_] = tau_x;
        (*a_sf)[n*nv+_YDIR_] = tau_y;
        (*a_sf)[n*nv+_ZDIR_] = tau_z;

        (*a_sf)[n*nv+_ZDIR_+1] = sqrt(tau_x*tau_x + tau_y*tau_y + tau_z*tau_z);
      }

    } else {

      _ArraySetValue_((*a_sf), nv*nfacets_local, 0.0);

    }

  }

  return 0;
}

/*! Write the surface data on the immersed body to a ASCII Tecplot file. */
static int WriteSurfaceData(  void*               a_m,              /*!< MPI object of type #MPIVariables */
                              void*               a_ib,             /*!< Immersed body object of type #ImmersedBoundary */
                              const double* const a_p_surface,      /*!< array with local surface pressure data */
                              const double* const a_T_surface,      /*!< array with local surface temperature data */
                              const double* const a_ngrad_p_surface,/*!< array with local normal gradient of surface pressure data */
                              const double* const a_ngrad_T_surface,/*!< array with local normal gradient of surface temperature data */
                              const double* const a_shear,          /*!< array with local shear data */
                              char*               a_filename        /*!< Name of file to write */
                            )
{
  MPIVariables *mpi = (MPIVariables*) a_m;
  ImmersedBoundary *IB  = (ImmersedBoundary*) a_ib;
  int ierr;

#ifndef serial
  MPI_Status status;
#endif

  /* collect the surface data into global arrays */
  double* p_surface_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, a_p_surface, &p_surface_g, 1);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }
  double* T_surface_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, a_T_surface, &T_surface_g, 1);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }
  double* ngrad_p_surface_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, a_ngrad_p_surface, &ngrad_p_surface_g, 1);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }
  double* ngrad_T_surface_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, a_ngrad_T_surface, &ngrad_T_surface_g, 1);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }
  double* shear_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, a_shear, &shear_g, 4);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }

  /* Rank 0 writes the file */
  if (!mpi->m_rank) {

    int nfacets_global = IB->m_body->m_nfacets;
    const Facet3D* const facets = IB->m_body->m_surface;

    FILE *out;
    out = fopen(a_filename,"w");
    fprintf(out,"TITLE = \"Surface data created by HyPar.\"\n");
    fprintf(out,"VARIABLES = \"X\", \"Y\", \"Z\", ");
    fprintf(out,"\"Surface_Pressure\", ");
    fprintf(out,"\"Surface_Temperature\", ");
    fprintf(out,"\"Normal_Grad_Surface_Pressure\", ");
    fprintf(out,"\"Normal_Grad_Surface_Temperature\", ");
    fprintf(out,"\"Shear_x\", ");
    fprintf(out,"\"Shear_y\", ");
    fprintf(out,"\"Shear_z\", ");
    fprintf(out,"\"Shear_magn\"");
    fprintf(out,"\n");
    fprintf(out,"ZONE N = %d, E = %d, DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n",3*nfacets_global,nfacets_global);

    for (int n = 0; n < nfacets_global; n++) {
      fprintf(  out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                facets[n].m_x1,
                facets[n].m_y1,
                facets[n].m_z1,
                p_surface_g[n],
                T_surface_g[n],
                ngrad_p_surface_g[n],
                ngrad_T_surface_g[n],
                shear_g[4*n+_XDIR_],
                shear_g[4*n+_YDIR_],
                shear_g[4*n+_ZDIR_],
                shear_g[4*n+_ZDIR_+1] );
      fprintf(  out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                facets[n].m_x2,
                facets[n].m_y2,
                facets[n].m_z2,
                p_surface_g[n],
                T_surface_g[n],
                ngrad_p_surface_g[n],
                ngrad_T_surface_g[n],
                shear_g[4*n+_XDIR_],
                shear_g[4*n+_YDIR_],
                shear_g[4*n+_ZDIR_],
                shear_g[4*n+_ZDIR_+1] );
      fprintf(  out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                facets[n].m_x3,
                facets[n].m_y3,
                facets[n].m_z3,
                p_surface_g[n],
                T_surface_g[n],
                ngrad_p_surface_g[n],
                ngrad_T_surface_g[n],
                shear_g[4*n+_XDIR_],
                shear_g[4*n+_YDIR_],
                shear_g[4*n+_ZDIR_],
                shear_g[4*n+_ZDIR_+1] );
    }
    for (int n = 0; n < nfacets_global; n++) fprintf(out,"%d %d %d\n",3*n+1,3*n+2,3*n+3);
    fclose(out);
  }

  if (p_surface_g) free(p_surface_g);
  if (T_surface_g) free(T_surface_g);
  if (ngrad_p_surface_g) free(ngrad_p_surface_g);
  if (ngrad_T_surface_g) free(ngrad_T_surface_g);
  if (shear_g) free(shear_g);

  return 0;
}


/*! Calculate the aerodynamic forces on the immersed body surface and write them
    to file
*/
int NavierStokes3DIBForces( void*   a_s,  /*!< Solver object of type #HyPar */
                            void*   a_m,  /*!< MPI object of type #MPIVariables */
                            double  a_t /*!< Current simulation time */ )
{
  HyPar             *solver  = (HyPar*)          a_s;
  MPIVariables      *mpi     = (MPIVariables*)   a_m;
  NavierStokes3D    *physics = (NavierStokes3D*) solver->m_physics;
  ImmersedBoundary  *IB      = (ImmersedBoundary*) solver->m_ib;
  int ierr;

  if (!solver->m_flag_ib) return(0);

  int npts = solver->m_npoints_local_wghosts;

  double* pressure = (double*) calloc (npts, sizeof(double));
  ierr = NavierStokes3DComputePressure(pressure, solver->m_u, solver);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  NavierStokes3DComputePressure() returned with error.\n");
    return 1;
  }
  double* temperature = (double*) calloc(npts, sizeof(double));
  ierr = NavierStokes3DComputeTemperature(temperature, solver->m_u, solver);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  NavierStokes3DComputeTemperature() returned with error.\n");
    return 1;
  }

  /* Compute surface pressure */
  double* p_surface = NULL;
  ierr = IBComputeFacetVar(solver, mpi, pressure, 1, &p_surface);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  IBComputeFacetVar() returned with error.\n");
    return 1;
  }
  /* Compute surface temperature */
  double* T_surface = NULL;
  ierr = IBComputeFacetVar(solver, mpi, temperature, 1, &T_surface);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  IBComputeFacetVar() returned with error.\n");
    return 1;
  }
  /* Compute normal surface pressure gradient */
  double *ngrad_p_surface = NULL;
  ierr = IBComputeNormalGradient(solver, mpi, pressure, 1, &ngrad_p_surface);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  IBComputeNormalGradient() returned with error.\n");
    return 1;
  }
  /* Compute normal temperature gradient */
  double *ngrad_T_surface = NULL;
  ierr = IBComputeNormalGradient(solver, mpi, temperature, 1, &ngrad_T_surface);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  IBComputeNormalGradient() returned with error.\n");
    return 1;
  }
  /* Compute shear forces */
  double *shear = NULL;
  ierr = ComputeShear(solver, mpi, solver->m_u, &shear);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  ComputeShear() returned with error.\n");
    return 1;
  }

  char surface_filename[_MAX_STRING_SIZE_] = "surface";
  if (solver->m_nsims == 1) {
    if (!strcmp(solver->m_op_overwrite,"no")) {
      strcat(surface_filename,solver->m_filename_index);
    }
  } else {
    char index[_MAX_STRING_SIZE_];
    GetStringFromInteger(solver->m_my_idx, index, (int)log10(solver->m_nsims)+1);
    strcat(surface_filename, "_");
    strcat(surface_filename, index);
    strcat(surface_filename, "_");
  }
  strcat(surface_filename,".dat");
  if (!mpi->m_rank) {
    printf("Writing immersed body surface data file %a_s.\n",surface_filename);
  }
  ierr = WriteSurfaceData(  mpi,
                            IB,
                            p_surface,
                            T_surface,
                            ngrad_p_surface,
                            ngrad_T_surface,
                            shear,
                            surface_filename );
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  WriteSurfaceData() returned with error\n");
    return 1;
  }

  free(pressure);
  free(temperature);
  if (p_surface) free(p_surface);
  if (T_surface) free(T_surface);
  if (ngrad_p_surface) free(ngrad_p_surface);
  if (ngrad_T_surface) free(ngrad_T_surface);
  if (shear) free(shear);

  return 0;
}
