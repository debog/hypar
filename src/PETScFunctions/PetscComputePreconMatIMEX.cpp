/*! @file PetscComputePreconMatIMEX.cpp
    @brief Contains the function to assemble the preconditioning matrix
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <arrayfunctions.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscComputePreconMatIMEX"

/*! Compute and assemble the preconditioning matrix for the implicit-explicit (IMEX) time integration
    of the governing equations: Right now, it just calls PetscComputePreconMatImpl() */
int PetscComputePreconMatIMEX(Mat Pmat,   /*!< Preconditioning matrix to construct */
                              Vec Y,      /*!< Solution vector */
                              void *ctxt  /*!< Application context */ )
{
  /* Same implementation as PetscComputePreconMatImpl() */
  return(PetscComputePreconMatImpl(Pmat,Y,ctxt));
}

#endif
