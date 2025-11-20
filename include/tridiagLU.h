/*! @file tridiagLU.h
    @brief Header file for TridiagLU
    @author Debojyoti Ghosh
*/

#ifndef _TRIDIAGLU_H_
#define _TRIDIAGLU_H_

/*

  Parallel direct solver for tridiagonal systems

  TridiagLU(a,b,c,x,n,ns,r,m) - Parallel tridiagonal solver

    Solves the tridiagonal system in parallel by reordering the
    points such that the first point of each subdomain is placed
    at the end.

    The interior points are eliminated in parallel, resulting in
    a reduced system consisting of the first point of each sub-
    domain.

    This reduced system is solved either by the gather-and-
    solve (tridiagLUGS) or the recursive-doubling (tridiagLURD)
    algorithms.

  TridiagLUGS(a,b,c,x,n,ns,r,m) - Tridiagonal solver based on
                                  "gather and solve"

    Each of the "ns" systems is gathered on one processor,
    solved in serial, and the solution scattered back. The
    parallelism lies in solving the "ns" different systems
    by multiple processors (i.e., each processor solves
    ~ns/nproc number of systems in serial).

  Arguments:-
    a   [0,ns-1]x[0,n-1] double*          subdiagonal entries
    b   [0,ns-1]x[0,n-1] double*          diagonal entries
    c   [0,ns-1]x[0,n-1] double*          superdiagonal entries
    x   [0,ns-1]x[0,n-1] double*          right-hand side (solution)
    n                    int              local size of the system
    ns                   int              number of systems to solve
    r                    TridiagLU*       structure containing paramters
                                          for the tridiagonal solve and
                                          the walltimes:
                                            total_time
                                            stage1_time
                                            stage2_time
                                            stage3_time
                                            stage4_time
                        ** Note that these are process-specific. Calling
                           function needs to do something to add/average
                           them to get some global value.
                        ** Can be NULL if runtimes are not needed.
    m                   MPI_Comm*       MPI Communicator

  For a,b,c, and x, [0,ns-1] is the inner loop, i.e., the i-th row of the
  d-th system is a[i*ns+d], b[i*ns+d], c[i*ns+d] and x[i*ns+d].

  Return value (int) -> 0 (successful solve), -1 (singular system)

  Note:-
    x contains the final solution (right-hand side is replaced)
    a,b,c are not preserved
    On rank=0,        a[0*ns+d] has to be zero for all d.
    On rank=nproc-1,  c[(n-1)*ns+d] has to be zero for all d.

  For a serial tridiagonal solver, compile with the flag "-Dserial"
  or send NULL as the argument for the MPI communicator.

*/

/*! Jacobi method \sa TridiagIterJacobi(), BlockTridiagIterJacobi() */
#define _TRIDIAG_JACOBI_  "jacobi"
/*! "Gather-and-solve" method \sa tridiagLUGS */
#define _TRIDIAG_GS_      "gather-and-solve"

/*! \def TridiagLU
    \brief Structure of variables used by TridiagLU

    This structure contains all the variables used by
    TridiagLU.
*/
typedef struct _tridiagLU_ {

  /* Parameters for TridiagLU() */

  /*! Choice of solver for solving the reduced system. May be #_TRIDIAG_JACOBI_
      or #_TRIDIAG_GS_.
  */
  char m_reducedsolvetype[50];

  int     m_evaluate_norm;  /*!< calculate norm at each iteration? (relevant only for iterative solvers) */
  int     m_maxiter;        /*!< maximum number of iterations (relevant only for iterative solvers)      */
  double  m_atol,           /*!< absolute tolerance (relevant only for iterative solvers)                */
          m_rtol;           /*!< relative tolerace (relevant only for iterative solvers)                 */
  int     m_exititer;       /*!< number of iterations it ran for (relevant only for iterative solvers)   */
  double  m_exitnorm;       /*!< error norm at exit (relevant only for iterative solvers)                */
  int     m_verbose;        /*!< print iterations and norms (relevant only for iterative solvers)        */

  double  m_total_time;     /*!< Total wall time in seconds */
  double  m_stage1_time;    /*!< Wall time (in seconds) for stage 1 of TridiagLU() or BlockTridiagLU() */
  double  m_stage2_time;    /*!< Wall time (in seconds) for stage 2 of TridiagLU() or BlockTridiagLU() */
  double  m_stage3_time;    /*!< Wall time (in seconds) for stage 3 of TridiagLU() or BlockTridiagLU() */
  double  m_stage4_time;    /*!< Wall time (in seconds) for stage 4 of TridiagLU() or BlockTridiagLU() */

#ifdef with_scalapack
  int m_blacs_ctxt;         /*!< Context variable for ScaLAPACK (relevant if compiled with ScaLAPACK
                               support (-Dwith_scalapack) \sa tridiagScaLPK */
#endif

} TridiagLU_Params;

int TridiagLU(double*,double*,double*,double*,int,int,void*,void*);
int TridiagLUGS(double*,double*,double*,double*,int,int,void*,void*);
int TridiagIterJacobi(double*,double*,double*,double*,int,int,void*,void*);
int TridiagLUInit(void*,void*);

/* Block solvers */
int BlockTridiagLU(double*,double*,double*,double*,int,int,int,void*,void*);
int BlockTridiagIterJacobi(double*,double*,double*,double*,int,int,int,void*,void*);

#ifdef with_scalapack
/* ScaLAPACK interface */
int TridiagScaLPK(double*,double*,double*,double*,int,int,void*,void*);
#endif

#endif
