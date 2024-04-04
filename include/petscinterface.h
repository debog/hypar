/*! @file petscinterface_cpp.h
    @brief Contains C++ function declarations for the PETSc time integration interface
    @author Debojyoti Ghosh
 */

#ifdef with_petsc

#ifndef _PETSC_INTERFACE_CPP_H_
#define _PETSC_INTERFACE_CPP_H_

#include <petscinterface_struct.h>

/* Copy Functions */
int TransferVecToPETSc(const double* const,Vec,void*,const int,const int);
int TransferVecFromPETSc(double* const,const Vec,void*,const int,const int);
int TransferMatToPETSc(void*,Mat,void*);

int PetscRegisterTIMethods (int);

/* Right and left -hand side functions */
PetscErrorCode PetscRHSFunctionExpl (TS,PetscReal,Vec,Vec,void*);
PetscErrorCode PetscRHSFunctionIMEX (TS,PetscReal,Vec,Vec,void*);
PetscErrorCode PetscIFunctionIMEX   (TS,PetscReal,Vec,Vec,Vec,void*);
PetscErrorCode PetscIFunctionImpl   (TS,PetscReal,Vec,Vec,Vec,void*);

PetscErrorCode PetscIJacobianIMEX(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
PetscErrorCode PetscJacobianFunctionIMEX_JFNK       (Mat,Vec,Vec);
PetscErrorCode PetscJacobianFunctionIMEX_Linear     (Mat,Vec,Vec);

PetscErrorCode PetscIJacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
PetscErrorCode PetscJacobianFunction_JFNK  (Mat,Vec,Vec);
PetscErrorCode PetscJacobianFunction_Linear(Mat,Vec,Vec);

int PetscGlobalDOF(void*);
int PetscCleanup(void*);
int PetscCreatePointList(void*);

/* Other functions */
PetscErrorCode PetscPreStage        (TS,PetscReal);
PetscErrorCode PetscPostStage       (TS,PetscReal,PetscInt,Vec*);
PetscErrorCode PetscPreTimeStep     (TS);
PetscErrorCode PetscPostTimeStep    (TS);

/* preconditioning functions */
int PetscComputePreconMatIMEX(Mat,Vec,void*);
int PetscComputePreconMatImpl(Mat,Vec,void*);
int PetscJacobianMatNonzeroEntriesImpl(Mat,int,void*);

/*! Function to compute any error estimates, if available */
PetscErrorCode PetscTimeError (TS);

#ifdef with_librom
/*! Function to compute initial guess of nonlinear solve using libROM */
PetscErrorCode PetscSetInitialGuessROM(SNES,Vec,void*);
#endif

#endif

#endif
