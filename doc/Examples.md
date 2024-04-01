Examples
========

\subpage basic_examples :
Some basic examples that are simulated using HyPar. They 
all use explicit time integration, and \b do \b not require HyPar to be compiled
with any external library. 
Most of them can be run on one or a small number of processors.

\subpage gpu_examples : Examples that use NVidia GPUs.

\subpage ib_examples : Examples that use the immersed boundary method to simulate various geometries.

\subpage librom_examples : Examples that use reduced-order modeling capabilities from libROM.

\subpage md_examples : Examples that use the multidomain functionality.

\subpage petsc_examples : 
Some examples that use implicit or semi-implicit (IMEX) time
integration methods implemented in PETSc. To run them, HyPar needs to be compiled \b with \b PETSc.

\subpage sg_examples : Examples that use the sparse grids spatial discretization.

\page basic_examples Basic Examples

The following are some basic examples that are simulated using HyPar. They 
all use explicit time integration, and \b do \b not require HyPar to be compiled
with any external libraries.

\subpage linear_adv_sine \n
\subpage linear_adv_sine_varyingadv \n
\subpage linear_adv_disc \n
\subpage linear_diff_sine 

\subpage burger_1d_sine

\subpage sod_shock_tube  \n
\subpage lax_shock_tube \n
\subpage shu_osher \n
\subpage sod_shock_tube_wgrav

\subpage sw_dambreak

\subpage linear_adv_gauss \n
\subpage linear_adv_2d_sine_varyingadv \n
\subpage linear_diff_sine2d

\subpage burgers_2d_sine

\subpage vlasov_1d1v_prescribed \n
\subpage vlasov_1d1v_twostreaminstability \n
\subpage vlasov_1d1v_landaudamping

\subpage euler2d_riemann4 \n
\subpage euler2d_riemann6 \n
\subpage euler2d_radexp \n
\subpage euler2d_vortex

\subpage euler2d_igwave \n
\subpage euler2d_rtb

\subpage navstok2d_ldsc \n
\subpage navstok2d_flatplate

\subpage sw_circdambreak \n
\subpage sw_latbelt

\subpage linear_adv_3dgauss

\subpage ns3d_isoturb \n
\subpage ns3d_bubble

\subpage numa3d_bubble

\page linear_adv_sine 1D Linear Advection - Sine Wave

Location: \b hypar/Examples/1D/LinearAdvection/SineWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Advection Equation (linearadr.h)

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 1D/LinearAdvection/SineWave/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/SineWave/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/SineWave/physics.inp

\b lusolver.inp (optional)
\include 1D/LinearAdvection/SineWave/lusolver.inp

\b weno.inp (optional)
\include 1D/LinearAdvection/SineWave/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/SineWave/aux/init.c

Output:
-------
After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

Solutions at t=0,0.5,1: The following figure is obtained 
by plotting \a op_00000.dat (initial), \a op_00005.dat (t=0.5),
and \a op_00010.dat (final). 
@image html Solution_1DLinearAdvSine.png

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 1D/LinearAdvection/SineWave/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 1D/LinearAdvection/SineWave/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 1D/LinearAdvection/SineWave/output.log


\page linear_adv_sine_varyingadv 1D Linear Advection - Sine Wave with Spatially-Varying Advection Speed

Location: \b Examples/1D/LinearAdvection/SineWave_NonConstantAdvection
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Linear Advection Equation (linearadr.h)

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Advection speed: \f$a\left(x\right) = 1 + \frac{1}{2}\sin^2\left(2\pi x\right)\f$

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 1D/LinearAdvection/SineWave_NonConstantAdvection/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/SineWave_NonConstantAdvection/boundary.inp

\b physics.inp (specifies filename for advection field)
\include 1D/LinearAdvection/SineWave_NonConstantAdvection/physics.inp
\b Note: Do not include the ".inp" extension in the filename above.

\b weno.inp (optional)
\include 1D/LinearAdvection/SineWave_NonConstantAdvection/weno.inp

To generate \b initial.inp and \b advection.inp, compile and run the 
following code in the run directory. 
\include 1D/LinearAdvection/SineWave_NonConstantAdvection/aux/init.c

Output:
-------
After running the code, there should be 101 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00100.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=2.5\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

The following animation shows the solution vs. grid point index:
@image html Solution_1DLinearAdvSine_VaryingAdv.gif

In addition to the usual output files, the linear advection physics 
module writes out the following files:
+ \b advection_00000.dat, ..., \b advection_00100.dat: These files
  share the same format as the solution output files \b op_*.dat 
  and contains the advection field \f$a\left(x\right)\f$.
The following figure shows the advection speed vs. grid point index:
@image html Solution_1DLinearAdvSine_VaryingAdv.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 1D/LinearAdvection/SineWave_NonConstantAdvection/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 1D/LinearAdvection/SineWave_NonConstantAdvection/output.log

\page linear_adv_disc 1D Linear Advection - Discontinuous Waves

Location: \b hypar/Examples/1D/LinearAdvection/DiscontinuousWaves
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Advection Equation (linearadr.h)

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$-1 \le x \le 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution:
  \f{equation}{
    u\left(x,0\right) = \left\{\begin{array}{lc} 
                          \exp\left(-\log\left(2\right)\frac{\left(x+7\right)^2}{0.0009}\right) & -0.8\le x \le -0.6 \\
                          1 & -0.4\le x \le -0.2 \\
                          1 - \left|10\left(x-0.1\right)\right| & 0\le x \le 0.2 \\
                          \sqrt{1-100\left(x-0.5\right)^2} & 0.4\le x \le 0.6 \\
                          0 & {\rm otherwise}
                        \end{array}\right.
  \f}

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 1D/LinearAdvection/DiscontinuousWaves/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/DiscontinuousWaves/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/DiscontinuousWaves/physics.inp

\b lusolver.inp (optional)
\include 1D/LinearAdvection/DiscontinuousWaves/lusolver.inp

\b weno.inp (optional)
\include 1D/LinearAdvection/DiscontinuousWaves/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/DiscontinuousWaves/aux/init.c

Output:
-------
After running the code, there should be two solution output
files \b op_00000.dat and \b op_00001.dat; the first one is
the initial solution, and the latter is the final solution.
Both these files are ASCII text (#HyPar::op_file_format is
set to \a text in \b solver.inp).

Final solution at t=2.0: The following figure is obtained 
by plotting \a op_00000.dat (initial) and \a op_00001.dat
(final). In both these files, the first column is grid 
index, the second column is x-coordinate, and the third 
column is the solution.
@image html Solution_1DLinearAdvDisc.png

Expected screen output:
\include 1D/LinearAdvection/DiscontinuousWaves/output.log

\page linear_diff_sine 1D Linear Diffusion - Sine Wave

Location: \b hypar/Examples/1D/LinearDiffusion/SineWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Diffusion Equation (linearadr.h)

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_) 
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Numerical Method:
  + Spatial discretization (parabolic): 2nd order (Interp2PrimSecondOrder()),
                                        conservative (ParabolicFunctionCons1Stage())
  + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 1D/LinearDiffusion/SineWave/solver.inp

\b boundary.inp
\include 1D/LinearDiffusion/SineWave/boundary.inp

\b physics.inp
\include 1D/LinearDiffusion/SineWave/physics.inp

To generate \b initial.inp (initial solution) and 
\b exact.inp (exact solution), compile and run the 
following code in the run directory. 
\include 1D/LinearDiffusion/SineWave/aux/exact.c

Output:
-------
After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=10\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

Solutions at t=0,2,4,6,8,10: The following figure is obtained 
by plotting \a op_00000.dat (t=0, initial), \a op_00002.dat (t=2),
\a op_00004.dat (t=4), \a op_00006.dat (t=6), \a op_00008.dat 
(t=8) and \a op_00010.dat (t=10, final). 
@image html Solution_1DLinearDiffSine.png

Since the exact solution is available at the final time 
, the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 1D/LinearDiffusion/SineWave/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 1D/LinearDiffusion/SineWave/output.log

\page burger_1d_sine 1D Inviscid Burgers Equation - Sine Wave

Location: \b hypar/Examples/1D/Burgers/SineWave
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Burgers Equation (burgers.h)

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \frac{1}{2 \pi t_s}\sin\left(2\pi x\right)\f$,
\f$t_s=2\f$ (time to shock formation)

Exact solution: 
  \f$u\left(x,y,t\right) = \frac{1}{2 \pi t_s} \sin\left(2\pi \left(x-u t\right)\right)\f$ (needs to be computed iteratively)\n

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 1D/Burgers/SineWave/solver.inp

\b boundary.inp
\include 1D/Burgers/SineWave/boundary.inp

\b physics.inp
\include 1D/Burgers/SineWave/physics.inp

\b lusolver.inp (optional)
\include 1D/Burgers/SineWave/lusolver.inp

\b weno.inp (optional)
\include 1D/Burgers/SineWave/weno.inp

To generate \b initial.inp (initial solution) 
and \b exact.inp (exact solution), compile and run the 
following code in the run directory. 
\include 1D/Burgers/SineWave/aux/init.c
\b Note: The exact solution is available only if the
final time is less than \f$t_s\f$ above.

Output:
-------
After running the code, there should be 8 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00007.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=2\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

Solutions at t=0,...,2: The following figure is obtained 
by plotting \a op_00000.dat (initial), ..., \a op_00007.dat (t=2).
@image html Solution_1DBurgersSine.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 1D/Burgers/SineWave/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 1D/Burgers/SineWave/output.log



\page sod_shock_tube 1D Euler Equations - Sod Shock Tube

Description: 
-------------------

Location: \b hypar/Examples/1D/Euler1D/SodShockTube 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + G.A. Sod, "A survey of several finite difference methods 
    for systems of nonlinear hyperbolic conservation laws," 
    J. Comput. Phys., 27, 1 (1978).
  + C. B. Laney, "Computational Gasdynamics", Cambridge 
    University Press, 1998.

Domain: \f$0 \le x \le 1.0\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 1, u = 0, p = 1\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.125, u = 0, p = 0.1\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include 1D/Euler1D/SodShockTube/solver.inp

\b boundary.inp
\include 1D/Euler1D/SodShockTube/boundary.inp

\b physics.inp
\include 1D/Euler1D/SodShockTube/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/SodShockTube/aux/init.c

Output:
-------
After running the code, there should be two solution output
files \b op_00000.dat and \b op_00001.dat; the first one is
the initial solution, and the latter is the final solution.
Both these files are ASCII text (#HyPar::op_file_format is
set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

Final solution at t=0.2: The following figure is obtained 
by plotting \a op_00001.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_1DSodShockTube.png

Expected screen output:
\include 1D/Euler1D/SodShockTube/output.log

\page lax_shock_tube 1D Euler Equations - Lax Shock Tube

Description: 
-------------------

Location: \b hypar/Examples/1D/Euler1D/LaxShockTube 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + P.D. Lax, "Weak solutions of nonlinear hyperbolic
    equations and their numerical computation," Comm.
    Pure App. Math., 7, 159 (1954).
  + C. B. Laney, "Computational Gasdynamics", Cambridge 
    University Press, 1998.

Domain: \f$0 \le x \le 1.0\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 0.445, \rho u = 0.311, e = 8.928\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.5, \rho u = 0, e = 1.4275\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include 1D/Euler1D/LaxShockTube/solver.inp

\b boundary.inp
\include 1D/Euler1D/LaxShockTube/boundary.inp

\b physics.inp
\include 1D/Euler1D/LaxShockTube/physics.inp

\b weno.inp (optional)
\include 1D/Euler1D/LaxShockTube/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/LaxShockTube/aux/init.c

Output:
-------
Note that \b iproc = 2 in \b solver.inp, so run this with
2 MPI ranks (or change \b iproc to 1). After running the code, 
there should be two solution output files \b op_00000.dat and 
\b op_00001.dat; the first one is the initial solution, and the 
latter is the final solution. Both these files are ASCII text 
(#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

Final solution at t=0.08: The following figure is obtained 
by plotting \a op_00001.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_1DLaxShockTube.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 1D/Euler1D/LaxShockTube/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for 
each component.

Expected screen output:
\include 1D/Euler1D/LaxShockTube/output.log


\page shu_osher 1D Euler Equations - Shu-Osher Problem

Description: 
------------

Location: \b hypar/Examples/1D/Euler1D/ShuOsherProblem 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + C.-W. Shu and S. Osher, "Efficient implementation of
    essentially non-oscillatory schemes ,II," J. Comput.
    Phys., 83 (1989), pp. 32–78

Domain: \f$-5 \le x \le 5\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ -5 \le x < -4\f$: \f$\rho = 27/7, u = 4\sqrt{35}/7, p = 31/3\f$
  + \f$ -4 \le x \le 5\f$: \f$\rho = 1 + 0.2\sin\left(5x\right), u = 0, p = 1\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include 1D/Euler1D/ShuOsherProblem/solver.inp

\b boundary.inp
\include 1D/Euler1D/ShuOsherProblem/boundary.inp

\b physics.inp
\include 1D/Euler1D/ShuOsherProblem/physics.inp

\b weno.inp (optional)
\include 1D/Euler1D/ShuOsherProblem/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/ShuOsherProblem/aux/init.c

Output:
-------
After running the code, 
there should be two solution output files \b op_00000.dat and 
\b op_00001.dat; the first one is the initial solution, and the 
latter is the final solution. Both these files are ASCII text 
(#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

Final solution at t=1.8: The following figure is obtained 
by plotting \a op_00001.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_1DShuOsherProblem.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 1D/Euler1D/ShuOsherProblem/conservation.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for 
each component.

Expected screen output:
\include 1D/Euler1D/ShuOsherProblem/output.log


\page sod_shock_tube_wgrav 1D Euler Equations - Sod Shock Tube with Gravitational Force

Description: 
-------------------

Location: \b hypar/Examples/1D/Euler1D/SodShockTubeWithGravity 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Euler equations (euler1d.h)

References: 
  + Xing, Y., Shu, C.-W., "High Order Well-Balanced WENO Scheme
    for the Gas Dynamics Equations Under Gravitational Fields",
    Journal of Scientific Computing, 54, 2013, pp. 645-662.

Domain: \f$0 \le x \le 1.0\f$, \a "slip-wall" (#_SLIP_WALL_) 
        boundary conditions (wall velocity is zero), with 
        \b uniform \b gravitational \b force \f$g=1\f$.

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 1, u = 0, p = 1\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.125, u = 0, p = 0.1\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order CRWENO (Interp1PrimFifthOrderCRWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
--------------------
\b solver.inp:
\include 1D/Euler1D/SodShockTubeWithGravity/solver.inp

\b boundary.inp
\include 1D/Euler1D/SodShockTubeWithGravity/boundary.inp

\b physics.inp
\include 1D/Euler1D/SodShockTubeWithGravity/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/SodShockTubeWithGravity/aux/init.c

Output:
-------
After running the code, there should be two solution output
files \b op_00000.dat and \b op_00001.dat; the first one is
the initial solution, and the latter is the final solution.
Both these files are ASCII text (#HyPar::op_file_format is
set to \a text in \b solver.inp).
In these files, the first column is grid index, the second 
column is x-coordinate, and the remaining columns are the 
solution components.

Final solution at t=0.2: The following figure is obtained 
by plotting \a op_00001.dat. Note that the output is in
terms of the conserved variables, so they have to converted
to the primitive variables (density, velocity, and pressure).
@image html Solution_1DSodShockTubeWithGravity.png

Expected screen output:
\include 1D/Euler1D/SodShockTubeWithGravity/output.log

\page sw_dambreak 1D Shallow Water Equations - Dam Breaking over Rectangular Bump

Location: \b hypar/Examples/1D/ShallowWater1D/DamBreakingRectangularBump
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Shallow Water Equations (shallowwater1d.h)

References:
  + Xing, Y., Shu, C.-W., "High order finite difference WENO
    schemes with the exact conservation property for the shallow
    water equations", Journal of Computational Physics, 208, 2005,
    pp. 206-227 (section 4.4).

Domain: \f$ 0 \le x \le 1500\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions

Initial solution:
\f{equation}{
  h\left(x\right) = \left\{\begin{array}{lc} 20 - b\left(x\right) & x <= 750 \\ 15 - b\left(x\right) & {\rm otherwise}\end{array}\right., u\left(x\right) = 0
\f}
where \f$b\left(x\right) = \left\{\begin{array}{lc} 8 & \left|x-750.0\right| \le 1500/8 \\  0 & {\rm  otherwise} \end{array}\right.\f$ is the bottom topography.

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

Note: in addition to the usual input files that HyPar needs, this
physical model needs the following input file(s):
+ \b topography.inp : file containing the bottom topography (same
  format as \b initial.inp).

\b solver.inp
\include 1D/ShallowWater1D/DamBreakingRectangularBump/solver.inp

\b boundary.inp
\include 1D/ShallowWater1D/DamBreakingRectangularBump/boundary.inp

\b physics.inp
\include 1D/ShallowWater1D/DamBreakingRectangularBump/physics.inp

\b weno.inp (optional)
\include 1D/ShallowWater1D/DamBreakingRectangularBump/weno.inp

To generate \b initial.inp and \b topography.inp, compile and run the 
following code in the run directory.
\include 1D/ShallowWater1D/DamBreakingRectangularBump/aux/init.c

Output:
-------
After running the code, there should be 5 solution output files
\b op_00000.dat, ..., \b op_00004.dat; the first one is the solution 
at \f$t=0\f$ and the final one is the solution at \f$t=60\f$. Since
#HyPar::op_overwrite is set to \a no in \b solver.inp, separate files 
are written for solutions at each output time. All the files are ASCII
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third and fourth columns are the two 
solution components.

In addition to the usual output files, the shallow water physics 
module writes out the following files:
+ \b topography_00000.dat, ..., \b topography_00004.dat: These files
  share the same format as the solution output files \b op_*.dat 
  and contains the topography \f$b\left(x\right)\f$.


Solutions at t=0,15,30,45,60: The following figure is obtained 
by plotting \a op_00000.dat (initial), \a op_00005.dat (t=0.5),
and \a op_00010.dat (final). \b Note: the figure plots 
\f$h\left(x\right)+b\left(x\right)\f$, i.e., it adds the third
column in \b op_*.dat and the third column in \b topography_*.dat.
@image html Solution_1DSWDamBreak.png

Expected screen output:
\include 1D/ShallowWater1D/DamBreakingRectangularBump/output.log


\page linear_adv_gauss 2D Linear Advection - Gaussian Pulse

Location: \b hypar/Examples/2D/LinearAdvection/GaussianPulse
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Linear Advection Equation (linearadr.h)

Domain: \f$-6 \le x < 6, -3 \le y < 3\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \exp\left[-\left(\frac{x^2}{2}+\frac{y^2}{2}\right)\right]\f$\n
Exact solution: \f$u\left(x,y,t\right) = u_0\left(x-a_xt,y-a_yt\right)\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/LinearAdvection/GaussianPulse/solver.inp

\b boundary.inp
\include 2D/LinearAdvection/GaussianPulse/boundary.inp

\b physics.inp (specifies \f$a_x\f$ and \f$a_y\f$)
\include 2D/LinearAdvection/GaussianPulse/physics.inp

\b lusolver.inp (optional)
\include 2D/LinearAdvection/GaussianPulse/lusolver.inp

\b weno.inp (optional)
\include 2D/LinearAdvection/GaussianPulse/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 2D/LinearAdvection/GaussianPulse/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 2

in \b solver.inp (i.e., 4 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 8 MPI ranks (or change \b iproc).

After running the code, there should be 21 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00020.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=12\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation was generated from the solution files:
@image html Solution_2DLinearAdvGauss.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/LinearAdvection/GaussianPulse/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/LinearAdvection/GaussianPulse/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 2D/LinearAdvection/GaussianPulse/output.log


\page linear_adv_2d_sine_varyingadv 2D Linear Advection - Sine Wave with Spatially-Varying Advection Speed

Location: \b hypar/Examples/2D/LinearAdvection/SineWave_NonConstantAdvection/case_01
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Linear Advection Equation (linearadr.h)

Domain: \f$0 \le x,y < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \cos\left(4\pi y\right)\f$\n

Advection speed: \f$a_x\left(x,y\right) = \sin\left(4\pi y\right)\f$, \f$a_y\left(x,y\right) = -\cos\left(4\pi x\right)\f$

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 2D/LinearAdvection/SineWave_NonConstantAdvection/case_01/solver.inp

\b boundary.inp
\include 2D/LinearAdvection/SineWave_NonConstantAdvection/case_01/boundary.inp

\b physics.inp (specifies filename for advection field)
\include 2D/LinearAdvection/SineWave_NonConstantAdvection/case_01/physics.inp
\b Note: Do not include the ".inp" extension in the filename above.

\b lusolver.inp (optional)
\include 2D/LinearAdvection/SineWave_NonConstantAdvection/case_01/lusolver.inp

\b weno.inp (optional)
\include 2D/LinearAdvection/SineWave_NonConstantAdvection/case_01/weno.inp

To generate \b initial.inp and \b advection.inp, compile and run the 
following code in the run directory. 
\include 2D/LinearAdvection/SineWave_NonConstantAdvection/case_01/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be 101 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00100.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=4\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation was generated from the solution files:
@image html Solution_2DLinearAdvSine_VaryingAdv.gif

In addition to the usual output files, the linear advection physics 
module writes out the following files:
+ \b advection_00000.dat, ..., \b advection_00100.dat: These files
  share the same format as the solution output files \b op_*.dat 
  and contains the advection field \f$a_x\left(x,y\right), a_y\left(x,y\right)\f$.
The following figure shows the vector plot of the velocity field:
@image html Solution_2DLinearAdvSine_VaryingAdv.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/LinearAdvection/SineWave_NonConstantAdvection/case_01/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 2D/LinearAdvection/SineWave_NonConstantAdvection/case_01/output.log


\page linear_diff_sine2d 2D Linear Diffusion - Sine Wave

Location: \b hypar/Examples/2D/LinearDiffusion/SineWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Linear Diffusion Equation (linearadr.h)

Domain: \f$0 \le x,y < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \sin\left(2\pi x\right)\sin\left(2\pi y\right)\f$\n
Exact solution: \f$u\left(x,y,t\right) = \exp\left[-\pi^2 \left(4\nu_x + 4\nu_y\right) t\right] u0\left(x,y\right)\f$.

Numerical Method:
 + Spatial discretization (parabolic): 2nd order (Interp2PrimSecondOrder()), 
                                       conservative (ParabolicFunctionCons1Stage())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/LinearDiffusion/SineWave/solver.inp

\b boundary.inp
\include 2D/LinearDiffusion/SineWave/boundary.inp

\b physics.inp (specifies \f$\nu_x\f$ and \f$\nu_y\f$)
\include 2D/LinearDiffusion/SineWave/physics.inp

To generate \b initial.inp (initial solution) and 
\b exact.inp (exact solution), compile and run the 
following code in the run directory. 
\include 2D/LinearDiffusion/SineWave/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=10\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation was generated from the solution files:
@image html Solution_2DLinearDiffSine.gif

Since the exact solution is available at the final time 
, the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/LinearDiffusion/SineWave/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 2D/LinearDiffusion/SineWave/output.log

\page burgers_2d_sine 2D Inviscid Burgers Equation - Sine Wave

Location: \b hypar/Examples/2D/Burgers/SineWave
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Burgers Equation (burgers.h)

Domain: \f$1 \le x,y < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: 
  \f$u\left(x,y,0\right) = \frac{1}{2 \pi t_s} \sin\left(2\pi x\right) \sin\left(2\pi y\right)\f$,\n
  \f$t_s=2\f$ (time to shock formation)

Exact solution: 
  \f$u\left(x,y,t\right) = \frac{1}{2 \pi t_s} \sin\left(2\pi \left(x-u t\right)\right) \sin\left(2\pi \left(y-u t\right)\right)\f$ (needs to be computed iteratively)\n

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 2D/Burgers/SineWave/solver.inp

\b boundary.inp
\include 2D/Burgers/SineWave/boundary.inp

\b physics.inp (specifies \f$a_x\f$ and \f$a_y\f$)
\include 2D/Burgers/SineWave/physics.inp

\b lusolver.inp (optional)
\include 2D/Burgers/SineWave/lusolver.inp

\b weno.inp (optional)
\include 2D/Burgers/SineWave/weno.inp

To generate \b initial.inp (initial solution) 
and \b exact.inp (exact solution), compile and run the 
following code in the run directory. 
\include 2D/Burgers/SineWave/aux/init.c
\b Note: The exact solution is available only if the
final time is less than \f$t_s\f$ above.

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 5 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00004.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=2\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following figures show the solutions at each output time:
  + t = 0 (initial)
@image html Solution_2DBurgersSineWave_0.png
  + t = 0.4)
@image html Solution_2DBurgersSineWave_1.png
  + t = 0.8
@image html Solution_2DBurgersSineWave_2.png
  + t = 1.2
@image html Solution_2DBurgersSineWave_3.png
  + t = 1.6 (final)
@image html Solution_2DBurgersSineWave_4.png

Since the exact solution is available at the final time 
the errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/Burgers/SineWave/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/Burgers/SineWave/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 2D/Burgers/SineWave/output.log

\page vlasov_1d1v_prescribed 2D (1D-1V) Vlasov Equation - Prescribed E-Field

Location: \b hypar/Examples/2D/Vlasov1D1V/PrescribedElectricField
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D (1D-1V) Vlasov Equation (vlasov.h)

Domain: 
  + \f$0 \le x < 2\pi\f$, \a "periodic" (#_PERIODIC_)
  + \f$-6 \le y < 6\f$, \a "dirichlet" (#_DIRICHLET_) (\f$f = 0\f$)


Initial solution: \f$f\left(x,v\right) = \left[1 + \frac{1}{10} \cos\left(x\right)\right] \exp\left(- \frac{1}{2}v^2\right)\f$

Prescribed electric field: \f$\frac{1}{10}\cos\left(x\right)\f$

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 2D/Vlasov1D1V/PrescribedElectricField/solver.inp

\b boundary.inp
\include 2D/Vlasov1D1V/PrescribedElectricField/boundary.inp

\b physics.inp
\include 2D/Vlasov1D1V/PrescribedElectricField/physics.inp

To generate \b initial.inp (initial solution), compile and run the 
following code in the run directory. 
\include 2D/Vlasov1D1V/PrescribedElectricField/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be 257 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00257.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=128\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation shows the evolution of the solution:
@image html Solution_1D1VVlasov_PrescribedE.gif

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/Vlasov1D1V/PrescribedElectricField/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 2D/Vlasov1D1V/PrescribedElectricField/output.log

\page vlasov_1d1v_twostreaminstability 2D (1D-1V) Vlasov Equation - Two-Stream Instability

\b Note: this example requires HyPar to be compiled with the FFTW library (https://www.fftw.org/).
While compiling, make sure that the configure step included the <B>\--enable-fftw</B> option
(and if necessary, <B>\--with-fftw-dir=/path/to/fftw</B>).

Location: \b hypar/Examples/2D/Vlasov1D1V/TwoStreamInstability
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D (1D-1V) Vlasov Equation (vlasov.h)

Domain: 
  + \f$0 \le x < 2\pi\f$, \a "periodic" (#_PERIODIC_)
  + \f$-7 \le y < 7\f$, \a "dirichlet" (#_DIRICHLET_) (\f$f = 0\f$)


Initial solution: 
\f$f\left(x,v\right) = \frac{8}{\sqrt{2 \pi T}}\left(1+\frac{1}{10}\cos\left(2k\pi\frac{x}{L}\right)\right)\left[\exp\left(-\frac{\left(v-2\right)^2}{2T}\right) + \exp\left(-\frac{\left(v+2\right)^2}{2T}\right)\right]\f$, \f$k=1,T=1,L=2\pi\f$.

Self-consistent electric field is computed by solving the Poisson equation
in a periodic domain using Fourier transforms. This examples *requires* HyPar
to be compiled with FFTW (http://www.fftw.org/).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 2D/Vlasov1D1V/TwoStreamInstability/solver.inp

\b boundary.inp
\include 2D/Vlasov1D1V/TwoStreamInstability/boundary.inp

\b physics.inp
\include 2D/Vlasov1D1V/TwoStreamInstability/physics.inp

To generate \b initial.inp (initial solution), compile and run the 
following code in the run directory. 
\include 2D/Vlasov1D1V/TwoStreamInstability/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be 201 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00201.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=40\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation shows the evolution of the solution:
@image html Solution_1D1VVlasov_TwoStreamInstability.gif

HyPar will also write out 1D text files named \b efield_00000.dat, ...,
\b efield_00200.dat that contain the electric field at each time
instant the solution was written. In this example, the electric
field is a 1D scalar quantity, and the columns of these files are:
grid index, x-coordinate, electric field value.
The following plot shows the electric field at \f$t=40\f$, obtained
by plotting \b efield_00200.dat:
@image html Solution_E_1D1VVlasov_TwoStreamInstability.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/Vlasov1D1V/TwoStreamInstability/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 2D/Vlasov1D1V/TwoStreamInstability/output.log

\page vlasov_1d1v_landaudamping 2D (1D-1V) Vlasov Equation - Landau Damping

\b Note: this example requires HyPar to be compiled with the FFTW library (https://www.fftw.org/).
While compiling, make sure that the configure step included the <B>\--enable-fftw</B> option
(and if necessary, <B>\--with-fftw-dir=/path/to/fftw</B>).

Location: \b hypar/Examples/2D/Vlasov1D1V/LandauDamping
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D (1D-1V) Vlasov Equation (vlasov.h)

Reference: Finn et al., A Numerical Study of Landau Damping with PETSc-PIC,
           https://doi.org/10.48550/arXiv.2303.12620

Domain: 
  + \f$0 \le x < 2\pi/k\f$, \a "periodic" (#_PERIODIC_)
  + \f$-10 v_{th} \le v < 10 v_{th}\f$, \a "dirichlet" (#_DIRICHLET_) (\f$f = 0\f$)

Initial solution: 
\f$f\left(x,v\right) = \frac{1}{\sqrt{2\pi v_{th}^2}} \exp \left(-\frac{v^2}{2v_{th}^2}\right) \left[1+\alpha\cos\left( kx \right)\right] \f$.

Self-consistent electric field is computed by solving the Poisson equation
in a periodic domain using Fourier transforms. This examples *requires* HyPar
to be compiled with FFTW (http://www.fftw.org/).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 2D/Vlasov1D1V/LandauDamping/solver.inp

\b boundary.inp
\include 2D/Vlasov1D1V/LandauDamping/boundary.inp

\b physics.inp
\include 2D/Vlasov1D1V/LandauDamping/physics.inp

To generate \b initial.inp (initial solution), compile and run the 
following code in the run directory. 
\include 2D/Vlasov1D1V/LandauDamping/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be 151 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00150.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=15\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions (distribution function) at each output time. 
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary().

Along with the solution files, the code will write out the 1D
potential and electric fields in text files: 
\b potential_00000.bin, \b potential_00001.bin, ... \b potential_00150.bin,
and
\b efield_00000.bin, \b efield_00001.bin, ... \b efield_00150.bin.
The columns of these files are: grid index, x-coordinate, and 
potential or electric field value, respectively.

The code \b 2D/Vlasov1D1V/LandauDamping/aux/MaximumEHistory.c reads the
the efield_<nnnnn>.dat files, calculates the maximum absolute value of the
electric field at that time instant, and writes a text file, \b efield_history.dat, 
with two columns: time and maximum absolute value of E-field. The following plot shows
damping of the electric field with time, obtained by plotting \b efield_history.dat:
@image html Solution_E_1D1VVlasov_LandauDamping.png

Expected screen output:
\include 2D/Vlasov1D1V/LandauDamping/output.log

\page euler2d_riemann4 2D Euler Equations - Riemann Problem Case 4

Location: \b hypar/Examples/2D/NavierStokes2D/Riemann2DCase4
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + P. Lax and X.-D. Liu, "Solution of two-dimensional Riemann
    problems of gas dynamics by positive schemes," SIAM J Sci 
    Comp 19 (1998), 319–340.

Domain: \f$-0.5 \le x,y \le 0.5\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions.

Initial solution: see \b Case \b 4 in the reference.

Numerical method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/Riemann2DCase4/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/Riemann2DCase4/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/Riemann2DCase4/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/Riemann2DCase4/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include 2D/NavierStokes2D/Riemann2DCase4/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=0.25\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following plot shows the density contours at the final time t=0.25, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokRiemann4.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/NavierStokes2D/Riemann2DCase4/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) for each component.

Expected screen output:
\include 2D/NavierStokes2D/Riemann2DCase4/output.log


\page euler2d_riemann6 2D Euler Equations - Riemann Problem Case 6

Location: \b hypar/Examples/2D/NavierStokes2D/Riemann2DCase6
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + P. Lax and X.-D. Liu, "Solution of two-dimensional Riemann
    problems of gas dynamics by positive schemes," SIAM J Sci 
    Comp 19 (1998), 319–340.

Domain: \f$-0.5 \le x,y \le 0.5\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions.

Initial solution: see \b Case \b 6 in the reference.

Numerical method:
 + Spatial discretization (hyperbolic): Characteristic-based 3rd order MUSCL (Interp1PrimThirdOrderMUSCLChar())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/Riemann2DCase6/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/Riemann2DCase6/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/Riemann2DCase6/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include 2D/NavierStokes2D/Riemann2DCase6/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=0.3\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following plot shows the density contours at the final time t=0.3, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokRiemann6.png

Expected screen output:
\include 2D/NavierStokes2D/Riemann2DCase6/output.log


\page euler2d_radexp 2D Euler Equations - Radial Expansion Wave

Location: \b hypar/Examples/2D/NavierStokes2D/RadialExpansionWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: http://www.as.dlr.de/hiocfd/case_c1.5.html

Domain: \f$-0.4 \le x,y \le 0.4\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions (supersonic outflow).

Initial solution: see reference above.

Numerical method:
 + Spatial discretization (hyperbolic): Characteristic-based 3rd order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/RadialExpansionWave/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/RadialExpansionWave/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/RadialExpansionWave/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/RadialExpansionWave/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include 2D/NavierStokes2D/RadialExpansionWave/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following plot shows the density contours at the final time t=1, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokRadialExpansion.png

Expected screen output:
\include 2D/NavierStokes2D/RadialExpansionWave/output.log


\page euler2d_vortex 2D Euler Equations - Isentropic Vortex Convection

Location: \b hypar/Examples/2D/NavierStokes2D/InviscidVortexConvection
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.1,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/InviscidVortexConvection/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InviscidVortexConvection/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InviscidVortexConvection/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/InviscidVortexConvection/weno.inp

\b lusolver.inp (optional)
\include 2D/NavierStokes2D/InviscidVortexConvection/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory.
\include 2D/NavierStokes2D/InviscidVortexConvection/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=20\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation shows the density contours as the vortex
convects over the domain:, 
@image html Solution_2DNavStokVortex.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.

Expected screen output:
\include 2D/NavierStokes2D/InviscidVortexConvection/output.log


\page euler2d_igwave 2D Euler Equations (with gravitational force) - Inertia-Gravity Waves

Location: \b hypar/Examples/2D/NavierStokes2D/InertiaGravityWave
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + W. C. Skamarock and J. B. Klemp, "Efficiency and accuracy of 
    the Klemp-Wilhelmson timesplitting technique", Monthly Weather 
    Review, 122 (1994), pp. 2623–2630.
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.1).

Domain: \f$0 \le x \le 300,000\,m, 0 \le y \le 10,000\,m\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions along \f$x\f$, "slip-wall" (#_SLIP_WALL_) boundary conditions along \f$y\f$.

Initial solution: See references above.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes2D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes2D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes2D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes2D::p0)
  + Hydrostatic balance type 3 (#NavierStokes2D::HB)
  + Brunt-Vaisala frequency 0.01 (#NavierStokes2D::N_bv)

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/InertiaGravityWave/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InertiaGravityWave/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InertiaGravityWave/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/InertiaGravityWave/weno.inp
\b Note: \a no_limiting is set to 1, i.e., since this is a 
smooth problem, WENO limiting is turned off, and the spatial
discretization uses a 5th order polynomial.

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/InertiaGravityWave/aux/init.c

Output:
-------
Note that \b iproc is set to 

      12 1

in \b solver.inp (i.e., 12 processors along \a x, and 1
processor along \a y). Thus, this example should be run
with 12 MPI ranks (or change \b iproc).

After running the code, there should be 26 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00025.bin; 
the first one is the solution at \f$t=0s\f$ and the final one
is the solution at \f$t=3000s\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/InertiaGravityWave/aux/PostProcess.c

The following plot shows the potential temperature perturbation
contours at the final time t=3000. It was plotted using VisIt
(https://wci.llnl.gov/simulation/computer-codes/visit/) with 
tecplot2d format chosen in the above postprocessing code.
@image html Solution_2DNavStokIGWave.png

If the postprocessing code above was used to write out files in 
text format, the following MATLAB script can be used to generate 
plots and visualize the solution:
\include 2D/NavierStokes2D/InertiaGravityWave/PlotSolution.m

Expected screen output:
\include 2D/NavierStokes2D/InertiaGravityWave/output.log



\page euler2d_rtb 2D Euler Equations (with gravitational force) - Rising Thermal Bubble

Location: \b hypar/Examples/2D/NavierStokes2D/RisingThermalBubble
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.2).

Domain: \f$0 \le x,y \le 1000\,m\f$, 
        "slip-wall" (#_SLIP_WALL_) boundary conditions on all sides.

Initial solution: See references above.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes2D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes2D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes2D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes2D::p0)
  + Hydrostatic balance type 2 (#NavierStokes2D::HB)

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/RisingThermalBubble/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/RisingThermalBubble/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/RisingThermalBubble/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/RisingThermalBubble/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/RisingThermalBubble/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 3

in \b solver.inp (i.e., 4 processors along \a x, and 3
processor along \a y). Thus, this example should be run
with 12 MPI ranks (or change \b iproc).

After running the code, there should be 41 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00040.bin; 
the first one is the solution at \f$t=0s\f$ and the final one
is the solution at \f$t=700s\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/RisingThermalBubble/aux/PostProcess.c

The following plot shows the potential temperature perturbation
contours at the final time t=700s. It was plotted using VisIt
(https://wci.llnl.gov/simulation/computer-codes/visit/) with 
tecplot2d format chosen in the above postprocessing code.
@image html Solution_2DNavStokRTB.png

The following animation was created using all the transient files:
@image html Solution_2DNavStokRTB.gif

If the postprocessing code above was used to write out files in 
text format, the following MATLAB script can be used to generate 
plots and visualize the solution:
\include 2D/NavierStokes2D/RisingThermalBubble/PlotSolution.m

Expected screen output:
\include 2D/NavierStokes2D/RisingThermalBubble/output.log


\page navstok2d_ldsc 2D Navier-Stokes Equations -  Lid-Driven Square Cavity

Location: \b hypar/Examples/2D/NavierStokes2D/LidDrivenCavity
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Navier-Stokes Equations (navierstokes2d.h)

Reference: 
+ Erturk, E., Corke, T.C., and Gokcol, C., "Numerical Solutions of
  2-D Steady Incompressible Driven Cavity Flow at High Reynolds Numbers",
  International Journal for Numerical Methods in Fluids, 48, 2005,
  http://dx.doi.org/10.1002/fld.953.
+ Ghia, U., Ghia, K.N., Shin, C.T., "High-Re Solutions for Incompressible
  Flow using the Navier-Stokes Equations and a Multigrid Method", Journal
  of Computational Physics, 48, 1982, http://dx.doi.org/10.1016/0021-9991(82)90058-4.

Note that this is an incompressible problem being solved here using the compressible
Navier-Stokes equations in terms of non-dimensional flow variables. The density and 
pressure are taken such that the speed of sound is 1.0, and the flow velocities 
specified in the initial and boundary conditions correspond to a characteristic 
Mach number of 0.1 (thus, they are 0.1 times the values in the above reference).

Domain: \f$0 \le x, y \le 1\f$

Boundary conditions:
+ No-slip wall BC on \f$x=0,1, 0 \le y \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ Moving no-slip wall BC on \f$y=1, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with specified 
  wall velocity of 0.1 in the x-direction).

Initial solution: \f$\rho=1, p=1/\gamma\f$. The velocities are specified according to 
the references above, but scaled by a factor of 0.1 to ensure that the characteristic
Mach number is 0.1.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 100, 1000, 3200\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the cavity length and width.
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.1\f$ (characteristic Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: 4th order, 4-stage Runge-Kutta (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/LidDrivenCavity/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/LidDrivenCavity/boundary.inp

\b physics.inp (\b Note: this file specifies \f$Re = 3200\f$,
change \a Re here for other Reynolds numbers.)
\include 2D/NavierStokes2D/LidDrivenCavity/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/LidDrivenCavity/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processor along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.dat, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
Since #HyPar::op_file_format is set to \a tecplot2d in \b solver.inp,
this file is in the ASCII Tecplot format and can be viewed in any
software that supports this format (e.g. VisIt).

Following plots show the streamlines (colored by the velocity magnitude)
for the solutions with Reynolds number 100, 1000, and 3200, respectively.
@image html Solution_2DNavStokLDSC_Re0100.png
@image html Solution_2DNavStokLDSC_Re1000.png
@image html Solution_2DNavStokLDSC_Re3200.png

Expected screen output (for Reynolds number 3200):
\include 2D/NavierStokes2D/LidDrivenCavity/output.log


\page navstok2d_flatplate 2D Navier-Stokes Equations -  Laminar Flow over Flat Plate

Location: \b hypar/Examples/2D/NavierStokes2D/FlatPlateLaminar/UniformGrid
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h)

Reference: 
+ Hirsch, "Numerical Computation of Internal & External Flows",
  Volume 1 (Fundamentals of Computational Fluid Dynamics), 2nd 
  Edition, Elsevier, Section 12.3.2 (page 618-625).

Domain: \f$-0.25 \le x \le 1\f$, \f$0 \le y \le 0.25\f$

Boundary conditions:
+ Symmetry BC on \f$y=0, -0.25 \le x < 0\f$ (imposed through "slip-wall" 
  #_SLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0
  wall velocity).
+ Subsonic outflow on \f$y=0.25\f$ (#_SUBSONIC_OUTFLOW_) with pressure
  \f$p=1/\gamma\f$.
+ Subsonic inflow on \f$x=0\f$ (#_SUBSONIC_INFLOW_) with density \f$\rho=1\f$,
  and velocity \f$(u,v) = (0.3,0)\f$.
+ Subsonic outflow on \f$x=1\f$ (#_SUBSONIC_OUTFLOW_) with pressure
  \f$p=1/\gamma\f$.

Initial solution: Uniform flow with \f$\rho=1, u=0.3, v=0, p=1/\gamma\f$.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 100,000\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the plate length .
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.3\f$ (freestream Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

This is a steady state problem - the solution converges to a steady laminar flow over a 
flat plate, and the residuals decrease with time. The skin friction coefficient is given by
\f{equation}{
  c_f = \frac {0.664} {\sqrt{Re_x}}, Re_x = \frac {\rho u x} {\mu}
\f}
(Blasius solution).

Input files required:
---------------------

\b solver.inp
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 4

in \b solver.inp (i.e., 2 processors along \a x, and 4
processor along \a y). Thus, this example should be run
with 8 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file (works only for 2D and 3D).
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The following plot was obtained by converting the binary file to the 
Tecplot format, and using VisIt to plot it (it shows the density 
and velocity):
@image html Solution_2DNavStokFlatPlate.png
The following plot shows a magnified view of the boundary layer:
@image html Solution_2DNavStokFlatPlateMagnified.png

The following file computes the skin friction as a function of the 
Reynolds number and writes to to a text file \b SkinFriction.dat with 4 
columns: Reynolds number, computed skin friction coefficient, exact skin 
friction coefficient (\f$0.664/\sqrt{Re_x}\f$), and normal velocity gradient 
on the plate surface (\f$\left.\partial u/\partial y\right|_{y=0}\f$). 
Compile and run it in the run directory.
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/aux/SkinFriction.c
The following figure showing the exact and computed skin friction coefficients
was obtained by plotting \b SkinFriction.dat:
@image html Solution_2DNavStokFlatPlateSkinFriction.png

Expected screen output:
\include 2D/NavierStokes2D/FlatPlateLaminar/UniformGrid/output.log


\page sw_circdambreak 2D Shallow Water Equations - Circular Dam Break

Location: \b hypar/Examples/2D/ShallowWater2D/CircularDamBreak
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Shallow Water Equations (shallowwater2d.h)

References:
  + Delis, Katsaounis, "Numerical solution of the two-dimensional
    shallow water equations by the application of relaxation methods", 
    Applied Mathematical Modelling, 29 (2005), pp. 754-783 (Section 6.3).

Domain: \f$0 \le x,y \le 50 \f$, \a "extrapolate" (#_EXTRAPOLATE_) boundary
        conditions everywhere

Initial solution:
\f{equation}{
  h\left(x,y\right) = \left\{\begin{array}{cc}10 & r \le 11.0\\1 & {\rm otherwise}\end{array}\right., u\left(x,y\right) = 0.
\f}
with the bottom topography as \f$b\left(x,y\right) = 0\f$, and 
\f$r^2=\left(x-25\right)^2+\left(y-25\right)^2\f$.

Other parameters:
 + \f$g=9.8\f$ (#ShallowWater2D::g)

Numerical Method:
  + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
  + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

Note: in addition to the usual input files that HyPar needs, this
physical model needs the following input file(s):
+ \b topography.inp : file containing the bottom topography (same
  format as \b initial.inp).

However, this case has a constant (zero) bottom topography, and thus
this file can be skipped.

\b solver.inp
\include 2D/ShallowWater2D/CircularDamBreak/solver.inp

\b boundary.inp
\include 2D/ShallowWater2D/CircularDamBreak/boundary.inp

\b physics.inp
\include 2D/ShallowWater2D/CircularDamBreak/physics.inp

\b weno.inp (optional)
\include 2D/ShallowWater2D/CircularDamBreak/weno.inp

To generate \b initial.inp and \b topography.inp, compile and run the 
following code in the run directory (note: topography.inp can be skipped
since this case involves a uniform topography).
\include 2D/ShallowWater2D/CircularDamBreak/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 6 solution output files
\b op_00000.dat, ..., \b op_00005.dat; the first one is the solution 
at \f$t=0\f$ and the final one is the solution at \f$t=0.69\f$. Since
#HyPar::op_overwrite is set to \a no in \b solver.inp, separate files 
are written for solutions at each output time. 

#HyPar::op_file_format is set to \a text in \b solver.inp, and
thus, all the files are in plain text (ASCII) format. In these files, 
the data is written out as: the first two columns are grid indices, the 
next two columns are the x and y coordinates, and the remaining columns 
are the three solution components.
#HyPar::op_file_format can also be set to \a tecplot2d to get the solution
files in the Tecplot format (that can be read using any visualization 
software that supports Tecplot format).

In addition to the usual output files, the shallow water physics 
module writes out the following files:
+ \b topography_00000.dat, ..., \b topography_00005.dat: These files
  share the same format as the solution output files \b op_*.dat 
  and contains the topography \f$b\left(x\right)\f$ (for this example,
  they do not matter since the topography is flat).

The following plot shows the final solution (water height):
@image html Solution_2DShallowWater_CircDamBreak.png
It was obtained by using the following MATLAB code to plot \b op_00005.dat:
\include 2D/ShallowWater2D/CircularDamBreak/PlotSolution.m

Expected screen output:
\include 2D/ShallowWater2D/CircularDamBreak/output.log


\page sw_latbelt 2D Shallow Water Equations - Latitude Belt Flow

Location: \b hypar/Examples/2D/ShallowWater2D/LatitudeBeltFlow
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Shallow Water Equations (shallowwater2d.h)

References:
 + Zhu, Et. al., "Variational Data Assimilation with a Variable Resolution
   Finite-Element Shallow-Water Equations Model", Monthly Weather Review,
   122, 1994, pp. 946--965.\n
   Govering equations and approximation of Coriolis forces: Eqns (2.1)-(2.4)\n
   Problem definition: Eqns. (4.1)-(4.3)

Domain: \f$0 \le x \le 6,000,000\,{\rm m}, 0\le y \le 4,400,000\,{\rm m}\f$, 
        \a "periodic" (#_PERIODIC_) boundaries along \a x, 
        \a "slip wall" (#_SW_SLIP_WALL_) boundaries along \a y.

Initial solution: See equations (4.1) and (4.2) in reference,
                  constant (zero) bottom topography

Other parameters:
 + \f$g=10\f$ (#ShallowWater2D::g)
 + \f$\hat{f}=0.0001\f$ (#ShallowWater2D::fhat)
 + \f$\beta = 1.5\times 10^{-11}\f$ (#ShallowWater2D::beta)
 + \f$D=4400000\f$ (#ShallowWater2D::D)

Numerical Method:
  + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
  + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

Note: in addition to the usual input files that HyPar needs, this
physical model needs the following input file(s):
+ \b topography.inp : file containing the bottom topography (same
  format as \b initial.inp).

However, this case has a constant (zero) bottom topography, and thus
this file can be skipped.

\b solver.inp
\include 2D/ShallowWater2D/LatitudeBeltFlow/solver.inp

\b boundary.inp
\include 2D/ShallowWater2D/LatitudeBeltFlow/boundary.inp

\b physics.inp
\include 2D/ShallowWater2D/LatitudeBeltFlow/physics.inp

\b weno.inp (optional)
\include 2D/ShallowWater2D/LatitudeBeltFlow/weno.inp

To generate \b initial.inp and \b topography.inp, compile and run the 
following code in the run directory (note: topography.inp can be skipped
since this case involves a uniform topography).
\include 2D/ShallowWater2D/LatitudeBeltFlow/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 101 solution output files
\b op_00000.dat, ..., \b op_00100.dat; the first one is the solution 
at \f$t=0\f$ and the final one is the solution at \f$t=360000\f$. Since
#HyPar::op_overwrite is set to \a no in \b solver.inp, separate files 
are written for solutions at each output time. 

#HyPar::op_file_format is set to \a text in \b solver.inp, and
thus, all the files are in plain text (ASCII) format. In these files, 
the data is written out as: the first two columns are grid indices, the 
next two columns are the x and y coordinates, and the remaining columns 
are the three solution components.
#HyPar::op_file_format can also be set to \a tecplot2d to get the solution
files in the Tecplot format (that can be read using any visualization 
software that supports Tecplot format).

The following plot shows an animation of the solution: 
@image html Solution_2DShallowWater_LatitudeBelt.gif
It was obtained by using the following MATLAB code:
\include 2D/ShallowWater2D/LatitudeBeltFlow/PlotSolution.m

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/ShallowWater2D/LatitudeBeltFlow/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error of each component(#HyPar::ConservationError).
Note that that conservation error for the water height is zero to machine 
precision (the non-zero conservation "error" for the momentum components
is due to the non-zero source term - Coriolis force).

Expected screen output:
\include 2D/ShallowWater2D/LatitudeBeltFlow/output.log


\page linear_adv_3dgauss 3D Linear Advection - Gaussian Pulse

Location: \b hypar/Examples/3D/LinearAdvection/GaussianPulse
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 3D Linear Advection Equation (linearadr.h)

Domain: \f$-3 \le x,y,z < 3\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \exp\left[-\left(\frac{x^2}{2}+\frac{y^2}{2}+\frac{z^2}{2}\right)\right]\f$\n
Exact solution: \f$u\left(x,y,t\right) = u_0\left(x-a_xt,y-a_yt,z-a_zt\right)\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 3D/LinearAdvection/GaussianPulse/solver.inp

\b boundary.inp
\include 3D/LinearAdvection/GaussianPulse/boundary.inp

\b physics.inp (specifies \f$a_x\f$, \f$a_y\f$, and \f$a_z\f$)
\include 3D/LinearAdvection/GaussianPulse/physics.inp

\b lusolver.inp (optional)
\include 3D/LinearAdvection/GaussianPulse/lusolver.inp

\b weno.inp (optional)
\include 3D/LinearAdvection/GaussianPulse/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 3D/LinearAdvection/GaussianPulse/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2 2

in \b solver.inp (i.e., 2 processors along \a x, \a y, and \a z).
Thus, this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=6\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
All the files are binary (#HyPar::op_file_format is set to \a binary in \b solver.inp).
The code \b hypar/Extras/BinaryToTecplot.c can be used to convert the binary
solution files to 3D Tecplot files that can be visualized in any software
supporting the Tecplot format. Similarly, the code \b hypar/Extras/BinaryToText.c 
can be used to convert the binary solution files to ASCII text files with the 
following data layout: the first three columns are grid indices, the next three
columns are x, y, and z coordinates, and the last column is the solution \a u.

The following animation was generated from the solution files 
(after converting to Tecplot format and plotting the iso-surface
in VisIt):
@image html Solution_3DLinearAdvGauss.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 3D/LinearAdvection/GaussianPulse/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 3D/LinearAdvection/GaussianPulse/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 3D/LinearAdvection/GaussianPulse/output.log


\page ns3d_isoturb 3D Navier-Stokes Equations - Isotropic Turbulence Decay

Location: \b hypar/Examples/3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: \f$0 \le x,y,z < 2\pi\f$, "periodic" (#_PERIODIC_) boundaries 
        everywhere.

Initial solution: Isotropic turbulent flow - The initial solution is 
specified in the Fourier space (with an energy
distribution similar to that of turbulent flow), and then transformed
to the physical space through an inverse transform.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 333.33\f$ (#NavierStokes3D::Re)
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes3D::Pr)
  + \f$M_\infty = 0.3\f$ (turbulence fluctuation Mach number) (#NavierStokes3D::Minf)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/physics.inp

\b weno.inp (optional)
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/weno.inp

\b lusolver.inp (optional)
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/lusolver.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\b Note: this code requires the \b FFTW library installed (http://www.fftw.org/).
To compile:

    gcc -I/path/to/fftw3.h -L/path/to/libfftw3.a -lfftw3 init.c

(see the FFTW website on ways to install it).
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2 1

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=5\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

To generate compute the energy spectrum from a given solution, compile 
and run the following code in the run directory. This code wants to read 
a file called \b op.bin, so make a symbolic link with that name pointing
to the solution file whose energy spectrum is to be computed. (If #HyPar::op_overwrite
is set to \a yes in \b solver.inp, then the solution file itself is called 
\b op.bin). Also note that the solution must be in binary format
(#HyPar::op_file_format must be \a binary in \b solver.inp). It will write out a 
ASCII text file \b spectrum.dat with two columns: \f$k\f$ and \f$E\left(k\right)\f$.
\b Note: this code requires the \b FFTW library installed (http://www.fftw.org/).
To compile:

    gcc -I/path/to/fftw3.h -L/path/to/libfftw3.a -lfftw3 fourier.c

(see the FFTW website on ways to install it). 
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/fourier.c

The following figure shows the initial and final (t=5) energy spectra:
@image html Solution_3DNavStok_IsoTurb_Spectrum.png

The following file computes the kinetic energy as a function of time
from the solution files. It writes out an ASCII text file \b energy.dat
with two colums: time and kinetic energy.
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/aux/kineticenergy.c

The following figure shows the kinetic energy decay:
@image html Solution_3DNavStok_IsoTurb_Energy.png

The code \b hypar/Extras/BinaryToTecplot.c can be used to convert the binary
solution files to 3D Tecplot files that can be visualized in any software
supporting the Tecplot format. Similarly, the code \b hypar/Extras/BinaryToText.c 
can be used to convert the binary solution files to ASCII text files with the 
following data layout: the first three columns are grid indices, the next three
columns are x, y, and z coordinates, and the remaining columns are the solution
components (\f$\rho, \rho u, \rho v, \rho w, e\f$).

The following figure shows the density iso-surface colored by the internal energy
(plotted in ParaView after converting the binary solution to a Tecplot file):
@image html Solution_3DNavStok_IsoTurb.png

Expected screen output:
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay/output.log



\page ns3d_bubble 3D Navier-Stokes Equations - Rising Thermal Bubble

Location: \b hypar/Examples/3D/NavierStokes3D/RisingThermalBubble_Config1
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: \f$0 \le x,y,z < 1000\,{\rm m}\f$, \a "slip-wall" (#_SLIP_WALL_) boundaries 
        everywhere, with zero wall velocity.

Reference:
  + Kelly, J. F., Giraldo, F. X., "Continuous and discontinuous Galerkin methods for a scalable
  three-dimensional nonhydrostatic atmospheric model: Limited-area mode", J. Comput. Phys., 231,
  2012, pp. 7988-8008 (see section 5.1.2).
  + Giraldo, F. X., Kelly, J. F., Constantinescu, E. M., "Implicit-Explicit Formulations of a
  Three-Dimensional Nonhydrostatic Unified Model of the Atmosphere (NUMA)", SIAM J. Sci. Comput., 
  35 (5), 2013, pp. B1162-B1194 (see section 4.1).

Initial solution: A warm bubble in cool ambient atmosphere. Note that in this example, the
gravitational forces and rising of the bubble is along the \a y-axis.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes3D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes3D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes3D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes3D::p0)
  + Hydrostatic balance type 2 (#NavierStokes3D::HB)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSP RK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b solver.inp
\include 3D/NavierStokes3D/RisingThermalBubble_Config1/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/RisingThermalBubble_Config1/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/RisingThermalBubble_Config1/physics.inp

\b weno.inp (optional)
\include 3D/NavierStokes3D/RisingThermalBubble_Config1/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/RisingThermalBubble_Config1/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4 4

in \b solver.inp (i.e., 4 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 64 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=200\,{\rm s}\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The following code (<B>Examples/3D/NavierStokes3D/RisingThermalBubble_Config1/aux/PostProcess.c</B>)
can be used to convert the binary solution file (with conserved variables 
\f$\rho,\rho u,\rho v,\rho w,e\f$) to Tecplot or plain text files with the primitive
and reference variables \f$\rho,u,v,w,P,\theta,\rho_0,P_0,\pi,\theta_0\f$ where the
subscript \f$0\f$ indicates the hydrostatic mean value.
\include 3D/NavierStokes3D/RisingThermalBubble_Config1/aux/PostProcess.c

The following figure shows the potential temperature iso-surface for the initial
and final solutions (plotted in VisIt):
@image html Solution_3DNavStok_Bubble3D.png

The file \b Extras/ExtractSlice.c can be used to extract a slice perpendicular to any dimension
at a specified location along that dimension. The extract slice is written out in the same
binary format as the original solutions files (with the same names op_xxxxx.bin) in a 
subdirectory called \b slices (\b Note: make the subdirectory called \a slices before running
this code). The following code (<B>Examples/3D/NavierStokes3D/RisingThermalBubble_Config1/aux/PostProcessSlice.c</B>)
can then be used (in the \b slices subdirectory) to convert the binary slice solution file (with conserved variables 
\f$\rho,\rho u,\rho v,\rho w,e\f$) to Tecplot or plain text files with the primitive
and reference variables \f$\rho,u,v,w,P,\theta,\rho_0,P_0,\pi,\theta_0\f$.
\b Note that it needs the relevant \b solver.inp and \b physics.inp that can be created
as follows:
+ Copy the original solver.inp and physics.inp to the slices subdirectory.
+ In solver.inp, set \b ndims as \b 2, and remove the component of \b size and \b iproc 
  corresponding to the dimension being eliminated while extracting the slices (in this case,
  it is \a z or the 3rd component).
+ In physics.inp, remove the component of \b gravity corresponding to the dimension perpendicular
  to the slice.

\include 3D/NavierStokes3D/RisingThermalBubble_Config1/aux/PostProcessSlice.c

The following figure shows the potential temperature \f$\theta\f$ along a slice at \f$z=500\,{\rm m}\f$ 
(plotted in VisIt):
@image html Solution_3DNavStok_Bubble.gif

Expected screen output:
\include 3D/NavierStokes3D/RisingThermalBubble_Config1/output.log



\page numa3d_bubble 3D NUMA (Nonhydrostatic Unified Model of the Atmosphere) Equations - Rising Thermal Bubble

Location: \b hypar/Examples/3D/NUMA/RisingThermalBubble
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 3D NUMA (Nonhydrostatic Unified Model of the Atmosphere) Equations (numa3d.h)

Domain: \f$0 \le x,y,z < 1000\,{\rm m}\f$, \a "numa-nfbc" (#_NO_FLUX_BC_) boundaries 
        everywhere.

Reference:
  + Kelly, J. F., Giraldo, F. X., "Continuous and discontinuous Galerkin methods for a scalable
  three-dimensional nonhydrostatic atmospheric model: Limited-area mode", J. Comput. Phys., 231,
  2012, pp. 7988-8008 (see section 5.1.2).
  + Giraldo, F. X., Kelly, J. F., Constantinescu, E. M., "Implicit-Explicit Formulations of a
  Three-Dimensional Nonhydrostatic Unified Model of the Atmosphere (NUMA)", SIAM J. Sci. Comput., 
  35 (5), 2013, pp. B1162-B1194 (see section 4.1).

Initial solution: A warm bubble in cool ambient atmosphere.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#Numa3D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#Numa3D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ (#Numa3D::g)
  + Angular speed of earth \f$\Omega = 0\f$ (#Numa3D::Omega)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#Numa3D::Pref)
  + Reference temperature (at zero altitude) \f$T_{ref} = 300\f$ (#Numa3D::Tref)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b solver.inp
\include 3D/NUMA/RisingThermalBubble/solver.inp

\b boundary.inp
\include 3D/NUMA/RisingThermalBubble/boundary.inp

\b physics.inp
\include 3D/NUMA/RisingThermalBubble/physics.inp

\b weno.inp (optional)
\include 3D/NUMA/RisingThermalBubble/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NUMA/RisingThermalBubble/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2 2

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 2 processor along \a z). Thus, 
this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be 41 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00040.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=200\,{\rm s}\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The code \b hypar/Extras/BinaryToTecplot.c can be used to convert the binary
solution files to 3D Tecplot files that can be visualized in any software
supporting the Tecplot format. Similarly, the code \b hypar/Extras/BinaryToText.c 
can be used to convert the binary solution files to ASCII text files with the 
following data layout: the first three columns are grid indices, the next three
columns are x, y, and z coordinates, and the remaining columns are the solution
components (\f$\rho{'}, \rho u, \rho v, \rho w, \Theta{'}\f$).

The following figure shows the density pertubation iso-surface for the initial
and final solutions (plotted in VisIt):
@image html Solution_3DNUMA_Bubble3D.png

The code \b hypar/Extras/ExtractSlice.c can be used to extract a slice along at a specified location
along any of the dimensions (\b Note: make a subdirectory \a slices before using this).
It will write out 2D slice solutions in binary format in the \a slices subdirectory
with the same names as the original solution files. The codes \b hypar/Extras/BinaryToTecplot.c
and \b hypar/Extras/BinaryToText.c can then be used to convert these to text or Tecplot
formats for visualization. (\b Note: A copy of the original \a solver.inp is required in the \a slices subdirectory
with \b ndims reduced by 1, and \b size and \b iproc with the appropriate component removed such that
they have 2 components.)

The following figure shows the density pertubation along a slice at \f$y=500\,{\rm m}\f$ 
(plotted in VisIt):
@image html Solution_3DNUMA_Bubble.gif

Expected screen output:
\include 3D/NUMA/RisingThermalBubble/output.log



\page petsc_examples PETSc Examples

The following are some examples that use explicit, implicit or semi-implicit (IMEX) time
integration methods implemented in PETSc (https://petsc.org/release/). To run them, 
HyPar needs to be compiled \b with \b PETSc. Familiarity with using PETSc is assumed.

\b Note: 
+ In general, any example or simulation can use PETSc time-integrators (assuming
  HyPar is compiled with PETSc) by specifying the PETSc inputs through a 
  <B>.petscrc</B> file, similar to the ones in the examples below. 
  The following file is an example of a .petscrc file (with explanatory comments). 
  - <B>hypar/Examples/PETScInputs/.petscrc_Example</B> 
+ The PETSc example directories have a file <B>.petscrc</B> and a sym link \b petscrc
pointing to .petscrc. The file .petscrc is the actual input file; \a petscrc is needed to generate
this documentation because Doxygen does not seem to include files with names starting with a dot!
+ The inputs in the .petscrc files (any lines not starting with a #) can also be specified in the
command line, for example, 
    
    /path/to/hypar/bin/HyPar -use-petscts -ts_type rk -ts_rk_type 4 ...


Explicit time integration:
--------------------------
\subpage linear_adv_sine_petsc \n
\subpage linear_adv_disc_petsc (with local truncation error-based adaptive time-step)

Implicit time integration:
--------------------------
\subpage linear_diff_sine_petsc \n
\subpage linear_diff_sine2d_petsc (with local truncation error-based adaptive time-step) \n
\subpage euler2d_vortex_petsc \n
\subpage ns2d_ldsc_petsc_implicit (with local truncation error-based adaptive time-step) \n
\subpage navstok2d_flatplate_petsc

Implicit-Explicit (IMEX) time integration:
------------------------------------------
\subpage euler2d_low_mach_vortex_petsc \n
\subpage euler2d_igwave_petsc \n
\subpage euler2d_rtb_petsc_imex (with local truncation error-based adaptive time-step) \n
\subpage ns2d_ldsc_petsc_imex \n
\subpage ns3d_bubble_petsc


\page linear_adv_sine_petsc 1D Linear Advection - Sine Wave

Location: \b hypar/Examples/1D/LinearAdvection/SineWave_PETSc
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Advection Equation (linearadr.h)

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: PETSc (SolvePETSc()) 
    - Method class: <B>Runge-Kutta</B> (TSRK - https://petsc.org/release/docs/manualpages/TS/TSRK.html)
    - Specific method: <B>Classical 4th-order RK</B> (TSRK4 - https://petsc.org/release/docs/manualpages/TS/TSRK4.html#TSRK4)

Input files required:
---------------------

<B>.petscrc</B>
\include 1D/LinearAdvection/SineWave_PETSc/petscrc

\b solver.inp
\include 1D/LinearAdvection/SineWave_PETSc/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/SineWave_PETSc/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/SineWave_PETSc/physics.inp

\b lusolver.inp (optional)
\include 1D/LinearAdvection/SineWave_PETSc/lusolver.inp

\b weno.inp (optional)
\include 1D/LinearAdvection/SineWave_PETSc/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/SineWave_PETSc/aux/init.c

Output:
-------
After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

Solutions at t=0,0.5,1: The following figure is obtained 
by plotting \a op_00000.dat (initial), \a op_00005.dat (t=0.5),
and \a op_00010.dat (final). 
@image html Solution_1DLinearAdvSinePETSc.png

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 1D/LinearAdvection/SineWave_PETSc/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 1D/LinearAdvection/SineWave_PETSc/output.log


\page linear_adv_disc_petsc 1D Linear Advection - Discontinuous Waves

Location: \b hypar/Examples/1D/LinearAdvection/DiscontinuousWaves_PETSc
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Advection Equation (linearadr.h)

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$-1 \le x \le 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution:
  \f{equation}{
    u\left(x,0\right) = \left\{\begin{array}{lc} 
                          \exp\left(-\log\left(2\right)\frac{\left(x+7\right)^2}{0.0009}\right) & -0.8\le x \le -0.6 \\
                          1 & -0.4\le x \le -0.2 \\
                          1 - \left|10\left(x-0.1\right)\right| & 0\le x \le 0.2 \\
                          \sqrt{1-100\left(x-0.5\right)^2} & 0.4\le x \le 0.6 \\
                          0 & {\rm otherwise}
                        \end{array}\right.
  \f}

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: PETSc (SolvePETSc()) 
    - Method class: <B>Runge-Kutta</B> (TSRK - https://petsc.org/release/docs/manualpages/TS/TSRK.html)
    - Specific method: <B>Dormand-Prince</B> (TSRK5DP - https://petsc.org/release/docs/manualpages/TS/TSRK5DP.html#TSRK5DP)

Input files required:
---------------------

<B>.petscrc</B>
\include 1D/LinearAdvection/DiscontinuousWaves_PETSc/petscrc

\b solver.inp
\include 1D/LinearAdvection/DiscontinuousWaves_PETSc/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/DiscontinuousWaves_PETSc/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/DiscontinuousWaves_PETSc/physics.inp

\b lusolver.inp (optional)
\include 1D/LinearAdvection/DiscontinuousWaves_PETSc/lusolver.inp

\b weno.inp (optional)
\include 1D/LinearAdvection/DiscontinuousWaves_PETSc/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/DiscontinuousWaves_PETSc/aux/init.c

Output:
-------
After running the code, there should be two solution output
files \b op_00000.dat and \b op_00001.dat; the first one is
the initial solution, and the latter is the final solution.
Both these files are ASCII text (#HyPar::op_file_format is
set to \a text in \b solver.inp).

Final solution at t=2.0: The following figure is obtained 
by plotting \a op_00000.dat (initial) and \a op_00001.dat
(final). In both these files, the first column is grid 
index, the second column is x-coordinate, and the third 
column is the solution.
@image html Solution_1DLinearAdvDiscPETSc.png

Expected screen output:
\include 1D/LinearAdvection/DiscontinuousWaves_PETSc/output.log


\page linear_diff_sine_petsc 1D Linear Diffusion - Sine Wave

Location: \b hypar/Examples/1D/LinearDiffusion/SineWave_PETSc
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Linear Diffusion Equation (linearadr.h)

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_) 
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Numerical Method:
  + Spatial discretization (parabolic): 2nd order (Interp2PrimSecondOrder()),
                                        conservative (ParabolicFunctionCons1Stage())
  + Time integration: PETSc (SolvePETSc()) - <B>Crank-Nicholson</B> (TSCN - https://petsc.org/release/docs/manualpages/TS/TSCN.html)

Input files required:
---------------------

<B>.petscrc</B>
\include 1D/LinearDiffusion/SineWave_PETSc/petscrc

\b solver.inp
\include 1D/LinearDiffusion/SineWave_PETSc/solver.inp

\b boundary.inp
\include 1D/LinearDiffusion/SineWave_PETSc/boundary.inp

\b physics.inp
\include 1D/LinearDiffusion/SineWave_PETSc/physics.inp

To generate \b initial.inp (initial solution) and 
\b exact.inp (exact solution), compile and run the 
following code in the run directory. 
\include 1D/LinearDiffusion/SineWave_PETSc/aux/exact.c

Output:
-------
After running the code, there should be 6 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00005.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=10\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

Solutions at t=0,2,4,6,8,10: The following figure is obtained 
by plotting the solution files.
@image html Solution_1DLinearDiffSinePETSc.png

Since the exact solution is available at the final time 
, the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 1D/LinearDiffusion/SineWave_PETSc/errors.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 1D/LinearDiffusion/SineWave_PETSc/output.log


\page linear_diff_sine2d_petsc 2D Linear Diffusion - Sine Wave

Location: \b hypar/Examples/2D/LinearDiffusion/SineWave_PETSc
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Linear Diffusion Equation (linearadr.h)

Domain: \f$0 \le x,y < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \sin\left(2\pi x\right)\sin\left(2\pi y\right)\f$\n
Exact solution: \f$u\left(x,y,t\right) = \exp\left[-\pi^2 \left(4\nu_x + 4\nu_y\right) t\right] u0\left(x,y\right)\f$.

Numerical Method:
 + Spatial discretization (parabolic): 2nd order (Interp2PrimSecondOrder()), 
                                       conservative (ParabolicFunctionCons1Stage())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html) -
     Although the ARK methods are semi-implicit (IMEX), here they are used in the "fully implicit" mode, i.e., the implicit 
     method is used to solve the complete equation (Note the flag \b -ts_arkimex_fully_implicit in <B>.petscrc</B>).
   - Specific method: <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX4.html)


Input files required:
---------------------

<B>.petscrc</B>
\include 2D/LinearDiffusion/SineWave_PETSc/petscrc

\b solver.inp
\include 2D/LinearDiffusion/SineWave_PETSc/solver.inp

\b boundary.inp
\include 2D/LinearDiffusion/SineWave_PETSc/boundary.inp

\b physics.inp (specifies \f$\nu_x\f$ and \f$\nu_y\f$)
\include 2D/LinearDiffusion/SineWave_PETSc/physics.inp

To generate \b initial.inp (initial solution) and 
\b exact.inp (exact solution), compile and run the 
following code in the run directory. 
\include 2D/LinearDiffusion/SineWave_PETSc/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 2 output
files: \b op_00000.dat (initial solution at t=0) and 
\b op_00001.dat (final solution at t=10). Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following plots show the initial and final solutions:
@image html Solution_2DLinearDiffSinePETSc.png

Since the exact solution is available at the final time 
, the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/LinearDiffusion/SineWave_PETSc/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 2D/LinearDiffusion/SineWave_PETSc/output.log


\page euler2d_rtb_petsc_imex 2D Euler Equations (with gravitational force) - Rising Thermal Bubble

Location: \b hypar/Examples/2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.2).

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See the following reference:
+ Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
  Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
  Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
  http://dx.doi.org/10.1137/15M1044369.

Domain: \f$0 \le x,y \le 1000\,m\f$, 
        "slip-wall" (#_SLIP_WALL_) boundary conditions on all sides.

Initial solution: See references above.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes2D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes2D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes2D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes2D::p0)
  + Hydrostatic balance type 2 (#NavierStokes2D::HB)

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX4.html)


Input files required:
---------------------

<B>.petscrc</B>
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/petscrc

\b solver.inp
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/aux/init.c

Output:
-------
Note that \b iproc is set to 

      1 1

in \b solver.inp (i.e., 1 processor along \a x, and 1
processor along \a y). Thus, this example should be run
with 1 MPI rank (or change \b iproc).

After running the code, there should be a output
file \b op.bin (final solution),  
since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, the file is written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/aux/PostProcess.c

The following plot shows the potential temperature perturbation
contours at the final time t=200s. 
@image html Solution_2DNavStokRTBPETSc.png
It was plotted by using the 
above postprocessing code to write the solution in text format
and using the following MATLAB script:
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/PlotSolution.m

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include 2D/NavierStokes2D/RisingThermalBubble_PETSc_IMEX/output.log


\page navstok2d_flatplate_petsc 2D Navier-Stokes Equations -  Laminar Flow over Flat Plate

Location: \b hypar/Examples/2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h)

Reference: 
+ Hirsch, "Numerical Computation of Internal & External Flows",
  Volume 1 (Fundamentals of Computational Fluid Dynamics), 2nd 
  Edition, Elsevier, Section 12.3.2 (page 618-625).

Domain: \f$-0.25 \le x \le 1\f$, \f$0 \le y \le 0.25\f$

Boundary conditions:
+ Symmetry BC on \f$y=0, -0.25 \le x < 0\f$ (imposed through "slip-wall" 
  #_SLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0
  wall velocity).
+ Subsonic outflow on \f$y=0.25\f$ (#_SUBSONIC_OUTFLOW_) with pressure
  \f$p=1/\gamma\f$.
+ Subsonic inflow on \f$x=0\f$ (#_SUBSONIC_INFLOW_) with density \f$\rho=1\f$,
  and velocity \f$(u,v) = (0.3,0)\f$.
+ Subsonic outflow on \f$x=1\f$ (#_SUBSONIC_OUTFLOW_) with pressure
  \f$p=1/\gamma\f$.

Initial solution: Uniform flow with \f$\rho=1, u=0.3, v=0, p=1/\gamma\f$.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 100,000\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the plate length .
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.3\f$ (freestream Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html) -
     Although the ARK methods are semi-implicit (IMEX), here they are used in the "fully implicit" mode, i.e., the implicit 
     method is used to solve the complete equation (Note the flag \b -ts_arkimex_fully_implicit in <B>.petscrc</B>).
   - Specific method: <B>ARK2e</B> (TSARKIMEX2E - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX2E.html)

This is a steady state problem - the solution converges to a steady laminar flow over a 
flat plate, and the residuals decrease with time. The skin friction coefficient is given by
\f{equation}{
  c_f = \frac {0.664} {\sqrt{Re_x}}, Re_x = \frac {\rho u x} {\mu}
\f}
(Blasius solution).

Input files required:
---------------------

<B>.petscrc</B>
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/petscrc

\b solver.inp
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 3

in \b solver.inp (i.e., 2 processor along \a x, and 3
processor along \a y). Thus, this example should be run
with 6 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file (works only for 2D and 3D).
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The following plot was obtained by converting the binary file to the 
Tecplot format, and using VisIt to plot it (it shows the density 
and velocity):
@image html Solution_2DNavStokFlatPlate.png
The following plot shows a magnified view of the boundary layer:
@image html Solution_2DNavStokFlatPlateMagnified.png

The following file computes the skin friction as a function of the 
Reynolds number and writes to to a text file \b SkinFriction.dat with 4 
columns: Reynolds number, computed skin friction coefficient, exact skin 
friction coefficient (\f$0.664/\sqrt{Re_x}\f$), and normal velocity gradient 
on the plate surface (\f$\left.\partial u/\partial y\right|_{y=0}\f$). 
Compile and run it in the run directory.
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/aux/SkinFriction.c
The following figure showing the exact and computed skin friction coefficients
was obtained by plotting \b SkinFriction.dat:
@image html Solution_2DNavStokFlatPlateSkinFrictionPETSc.png

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include 2D/NavierStokes2D/FlatPlateLaminar_PETSc_Implicit/output.log

\page ns2d_ldsc_petsc_imex 2D Navier-Stokes Equations -  Lid-Driven Square Cavity

Location: \b hypar/Examples/2D/NavierStokes2D/LidDrivenCavity_PETSc_IMEX
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Navier-Stokes Equations (navierstokes2d.h)

Reference: 
+ Erturk, E., Corke, T.C., and Gokcol, C., "Numerical Solutions of
  2-D Steady Incompressible Driven Cavity Flow at High Reynolds Numbers",
  International Journal for Numerical Methods in Fluids, 48, 2005,
  http://dx.doi.org/10.1002/fld.953.
+ Ghia, U., Ghia, K.N., Shin, C.T., "High-Re Solutions for Incompressible
  Flow using the Navier-Stokes Equations and a Multigrid Method", Journal
  of Computational Physics, 48, 1982, http://dx.doi.org/10.1016/0021-9991(82)90058-4.

Note that this is an incompressible problem being solved here using the compressible
Navier-Stokes equations in terms of non-dimensional flow variables. The density and 
pressure are taken such that the speed of sound is 1.0, and the flow velocities 
specified in the initial and boundary conditions correspond to a characteristic 
Mach number of 0.1 (thus, they are 0.1 times the values in the above reference).

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See:
+ Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
  Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
  Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
  http://dx.doi.org/10.1137/15M1044369.

Domain: \f$0 \le x, y \le 1\f$

Boundary conditions:
+ No-slip wall BC on \f$x=0,1, 0 \le y \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ Moving no-slip wall BC on \f$y=1, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with specified 
  wall velocity of 0.1 in the x-direction).

Initial solution: \f$\rho=1, p=1/\gamma\f$. The velocities are specified according to 
the references above, but scaled by a factor of 0.1 to ensure that the characteristic
Mach number is 0.1.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 3200\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the cavity length and width.
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.1\f$ (characteristic Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX4.html)

Input files required:
---------------------

<B>.petscrc</B>
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_IMEX/petscrc

\b solver.inp
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_IMEX/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_IMEX/boundary.inp

\b physics.inp (\b Note: this file specifies \f$Re = 3200\f$,
change \a Re here for other Reynolds numbers.)
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_IMEX/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_IMEX/aux/init.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processor along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.dat, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
Since #HyPar::op_file_format is set to \a tecplot2d in \b solver.inp,
this file is in the ASCII Tecplot format and can be viewed in any
software that supports this format (e.g. VisIt).

The following plot shows the streamlines (colored by the velocity magnitude):
@image html Solution_2DNavStokLDSC_Re3200_PETSc_IMEX.png

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_IMEX/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output (for Reynolds number 3200):
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_IMEX/output.log

\page ns3d_bubble_petsc 3D Navier-Stokes Equations - Rising Thermal Bubble

Location: \b hypar/Examples/3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: \f$0 \le x,y,z < 1000\,{\rm m}\f$, \a "slip-wall" (#_SLIP_WALL_) boundaries 
        everywhere, with zero wall velocity.

Reference:
  + Kelly, J. F., Giraldo, F. X., "Continuous and discontinuous Galerkin methods for a scalable
  three-dimensional nonhydrostatic atmospheric model: Limited-area mode", J. Comput. Phys., 231,
  2012, pp. 7988-8008 (see section 5.1.2).
  + Giraldo, F. X., Kelly, J. F., Constantinescu, E. M., "Implicit-Explicit Formulations of a
  Three-Dimensional Nonhydrostatic Unified Model of the Atmosphere (NUMA)", SIAM J. Sci. Comput., 
  35 (5), 2013, pp. B1162-B1194 (see section 4.1).

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See the following reference:
+ Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
  Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
  Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
  http://dx.doi.org/10.1137/15M1044369.

Initial solution: A warm bubble in cool ambient atmosphere. Note that in this example, the
gravitational forces and rising of the bubble is along the \a y-axis.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes3D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes3D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes3D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes3D::p0)
  + Hydrostatic balance type 2 (#NavierStokes3D::HB)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>ARK 2e</B> (TSARKIMEX2E - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX2E.html)

Input files required:
---------------------

<B>.petscrc</B>
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/petscrc

\b solver.inp
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/physics.inp

\b weno.inp (optional)
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4 4

in \b solver.inp (i.e., 4 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 64 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=200\,{\rm s}\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The binary solution file contains the conserved variables
(\f$\rho,\rho u,\rho v,\rho w,e\f$). The Python script
plotSolution.py calculates the primitive
and reference variables 
\f$\rho,u,v,w,P,\theta,\rho_0,P_0,\pi,\theta_0\f$ (where the
subscript \f$0\f$ indicates the hydrostatic mean value) and plots 
3D figures for \f$\theta\f$ (potential temperature).

The following figure shows the potential temperature at the final time. Colored
contours are plotted in the X-Y plane, and contour lines are plotted in the Z-Y
plane:
@image html Solution_3DNavStok_Bubble3D_PETSc.png

The code <B>Examples/3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/aux/PostProcess.c</B>
can be used to convert the binary solution files to Tecplot or plain text files with the 
primitive and reference variables.

The file \b Extras/ExtractSlice.c can be used to extract a slice perpendicular to any dimension
at a specified location along that dimension. The extract slice is written out in the same
binary format as the original solutions files (with the same names op_xxxxx.bin) in a 
subdirectory called \b slices (\b Note: make the subdirectory called \a slices before running
this code). The following code (<B>Examples/3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/aux/PostProcessSlice.c</B>)
can then be used (in the \b slices subdirectory) to convert the binary slice solution file (with conserved variables 
\f$\rho,\rho u,\rho v,\rho w,e\f$) to Tecplot or plain text files with the primitive
and reference variables \f$\rho,u,v,w,P,\theta,\rho_0,P_0,\pi,\theta_0\f$.
\b Note that it needs the relevant \b solver.inp and \b physics.inp that can be created
as follows:
+ Copy the original solver.inp and physics.inp to the slices subdirectory.
+ In solver.inp, set \b ndims as \b 2, and remove the component of \b size and \b iproc 
  corresponding to the dimension being eliminated while extracting the slices (in this case,
  it is \a z or the 3rd component).
+ In physics.inp, remove the component of \b gravity corresponding to the dimension perpendicular
  to the slice.

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX/output.log


\page euler2d_vortex_petsc 2D Euler Equations - Isentropic Vortex Convection

Location: \b hypar/Examples/2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.5,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: PETSc (SolvePETSc()) - <B>Crank-Nicholson</B> (TSCN - https://petsc.org/release/docs/manualpages/TS/TSCN.html)

Input files required:
---------------------

<B>.petscrc</B>
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/petscrc

\b solver.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/physics.inp

\b lusolver.inp (optional)
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory.
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=20\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following plot shows the density contours at the final time t=1, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokVortexPETSc.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.
\b Note that the conservation error depends on the accuracy with which the 
implicit systems are solved. Since the systems are solved with relaxed tolerances
here (see \a ksp_atol, \a ksp_rtol, \a snes_atol, \a snes_rtol in <B>.petscrc</B>),
the conservation error is zero within round-off.

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include 2D/NavierStokes2D/InviscidVortexConvection_PETSc_Implicit/output.log


\page euler2d_low_mach_vortex_petsc 2D Euler Equations - Low-Mach Isentropic Vortex Convection

Location: \b hypar/Examples/2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
           Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
           Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
           http://dx.doi.org/10.1137/15M1044369.

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See the above reference.

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.1,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX4.html)

Input files required:
---------------------

<B>.petscrc</B>
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/petscrc

\b solver.inp
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/physics.inp

To generate \b initial.inp (initial solution) and \b exact.inp (exact solution), 
compile and run the following code in the run 
directory:
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      2 2

in \b solver.inp (i.e., 2 processors along \a x, and 2
processors along \a y). Thus, this example should be run
with 4 MPI ranks (or change \b iproc).

After running the code, there should be 11 output
files \b op_00000.dat, \b op_00001.dat, ... \b op_00010.dat; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=100\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following plot shows the density contours at the final time t=1, 
obtained from plotting \b op_00010.dat:
@image html Solution_2DNavStokLowMachVortexPETSc.gif

Since the exact solution is available at the final time 
(\a exact.inp is a copy of \a initial.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation.dat:
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/conservation.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.
\b Note that the conservation error depends on the accuracy with which the 
implicit systems are solved (see \a ksp_atol, \a ksp_rtol, \a snes_atol, 
\a snes_rtol in <B>.petscrc</B>).

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include 2D/NavierStokes2D/LowMachVortexConvection_PETSc_IMEX/output.log


\page euler2d_igwave_petsc 2D Euler Equations (with gravitational force) - Inertia-Gravity Waves

Location: \b hypar/Examples/2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + W. C. Skamarock and J. B. Klemp, "Efficiency and accuracy of 
    the Klemp-Wilhelmson timesplitting technique", Monthly Weather 
    Review, 122 (1994), pp. 2623–2630.
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.1).

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See the following reference:
+ Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
  Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
  Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
  http://dx.doi.org/10.1137/15M1044369.

Domain: \f$0 \le x \le 300,000\,m, 0 \le y \le 10,000\,m\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions along \f$x\f$, "slip-wall" (#_SLIP_WALL_) boundary conditions along \f$y\f$.

Initial solution: See references above.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes2D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes2D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes2D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes2D::p0)
  + Hydrostatic balance type 3 (#NavierStokes2D::HB)
  + Brunt-Vaisala frequency 0.01 (#NavierStokes2D::N_bv)

Numerical method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX4.html)


Input files required:
---------------------

<B>.petscrc</B>
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX/petscrc

\b solver.inp
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX/aux/init.c

Output:
-------
Note that \b iproc is set to 

      12 1

in \b solver.inp (i.e., 12 processors along \a x, and 1
processor along \a y). Thus, this example should be run
with 12 MPI ranks (or change \b iproc).

After running the code, there should be 16 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00015.bin; 
the first one is the solution at \f$t=0s\f$ and the final one
is the solution at \f$t=3000s\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX/aux/PostProcess.c

The following plot shows the potential temperature perturbation
contours at the final time t=3000. It was plotted using VisIt
(https://wci.llnl.gov/simulation/computer-codes/visit/) with 
tecplot2d format chosen in the above postprocessing code.
@image html Solution_2DNavStokIGWavePETSc.png

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output:
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX/output.log



\page ib_examples Immersed Boundaries Examples

The following are some examples are use the immersed boundary methodology to solve the PDE in the presence of
various geometries. To use the immersed boundary implementation in HyPar, an STL (https://en.wikipedia.org/wiki/STL_%28file_format%29) 
representation of the immersed body is necessary. Note:
+ The immersed boundary method is \b only implemented for 3-dimensional simulations (#HyPar::ndims = 3).
+ It can be used with only those physical models that define an immersed boundary implementation (#HyPar::IBFunction()), for
  example, the 3D Navier-Stokes (#NavierStokes3D) equations (NavierStokes3DIBAdiabatic(), NavierStokes3DIBIsothermal()).

3D Navier-Stokes Equations:
---------------------------

\subpage ns3d_cylinder_steady_incompressible_viscous \n
\subpage ns3d_cylinder_unsteady_incompressible_viscous

\subpage ns3d_shock_cylinder_interaction

\subpage ns3d_sphere_steady_incompressible_viscous_adiabatic \n
\subpage ns3d_sphere_steady_incompressible_viscous_isothermal \n
\subpage ns3d_sphere_steady_incompressible_viscous_isothermal_warm \n
\subpage ns3d_sphere_unsteady_compressible_viscous_adiabatic \n
\subpage ns3d_sphere_unsteady_compressible_viscous_isothermal

\page ns3d_cylinder_steady_incompressible_viscous Steady, incompressible, viscous flow around a cylinder

Location: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the cylinder defined by [-2,6] X [-2,2],
        and a stretched grid beyond this zone.
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b hypar/Examples/STLGeometries/cylinder.stl)

The following images shows the grid and the cylinder:
@image html Domain3D_Cylinder.png
@image html Domain2D_Cylinder.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
  + Taneda, S., "Experimental Investigation of the Wakes behind Cylinders and Plates at Low 
    Reynolds Numbers," Journal of the Physical Society of Japan, Vol. 11, 302–307, 1956. 
  + Dennis, S. C. R., Chang, G.-Z., "Numerical solutions for steady flow past a circular
    cylinder at Reynolds numbers up to 100", Journal of Fluid Mechanics, 42 (3), 1970,
    pp. 471-489.

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 10,15,20\f$ (#NavierStokes3D::Re) (\b Note: 
    since the diameter of the cylinder is 2.0, the cylinder-diameter-based Reynolds number is 
    \f$Re_D = 2Re = 20,30,40\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/

\b solver.inp
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 10 (corresponding to \f$Re_D\f$ of 20). To try other Reynolds 
numbers, change it here.
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/aux/init.c

Output:
-------

Note that \b iproc is set to 

      2 2 1

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 4 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The file \b Extras/ExtractSlice.c can be used to extract a slice perpendicular to any dimension
at a specified location along that dimension. The extracted slice is written out in the same
binary format as the original solutions files (with the same names op_xxxxx.bin) in a 
subdirectory called \b slices (\b Note: make the subdirectory called \a slices before running
this code). The codes Extras/BinaryToTecplot.c and Extras/BinaryToText.c
can then be used (in the \b slices subdirectory) to convert the binary slice solution file
to Tecplot or plain text files.
\b Note that it needs the relevant \b solver.inp that can be created as follows:
+ Copy the original solver.inp to the slices subdirectory.
+ In solver.inp, set \b ndims as \b 2, and remove the component of \b size and \b iproc 
  corresponding to the dimension being eliminated while extracting the slices (in this case,
  it is \a z or the 3rd component).

The following figure shows the flow (density and streamlines) at \f$Re_D=20\f$:
@image html Solution_3DNavStokCylinder_ReD020.png
The parameter \a Re can be changed to 20 in \a physics.inp to run this simulation for
\f$Re_D=40\f$, and following figure shows the solution:
@image html Solution_3DNavStokCylinder_ReD040.png

The following plot shows the wake length (\f$L/D\f$) as a function of the Reynolds
number (\f$Re_D\f$) for the computed solutions and experimental results reported
in the reference above:
@image html Solution_3DNavStokCylinder_WL.png

In addition to the main solution, the code also writes out a file with the aerodynamic
forces on the immersed body. This file is called \a surface.dat (if #HyPar::op_overwrite
is "yes") or \a surface_nnnnn.dat (if #HyPar::op_overwrite is "no", "nnnnn" is a numerical
index) (in this example, the file \b surface.dat is written out). This is an ASCII file in 
the Tecplot format, where the immersed body and the forces on it are represented using the 
"FETRIANGLE" type. The following image shows the surface pressure on the cylinder (front-view):
@image html IBSurface_3DNavStokCylinder.png
Since this is a 2D simulation, the value of the surface pressure on the end-surfaces of the 
cylinder (zmin and zmax) are not physically relevant.

Expected screen output (for \f$Re_D = 20\f$):
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible/output.log

  
\page ns3d_cylinder_unsteady_incompressible_viscous Unsteady, incompressible, viscous flow around a cylinder (vortex shedding)

Location: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the cylinder defined by [-4,12] X [-2,2],
        and a stretched grid beyond this zone.
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b hypar/Examples/STLGeometries/cylinder.stl)

The following images shows the grid and the cylinder:
@image html Domain3D_Cylinder2.png
@image html Domain2D_Cylinder2.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
  + Taneda, S., "Experimental Investigation of the Wakes behind Cylinders and Plates at Low 
    Reynolds Numbers," Journal of the Physical Society of Japan, Vol. 11, 302–307, 1956. 

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 50\f$ (#NavierStokes3D::Re) (\b Note: 
    since the diameter of the cylinder is 2.0, the cylinder-diameter-based Reynolds number is 
    \f$Re_D = 2Re = 100\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/

\b solver.inp
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/aux/init.c

Output:
-------

Note that \b iproc is set to 

      4 2 1

in \b solver.inp (i.e., 4 processors along \a x, 2
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 8 MPI ranks (or change \b iproc).

After running the code, there should be 401 solution files \b op_00000.bin, 
\b op_00001.bin, ..., \b op_00400.bin. Since #HyPar::op_overwrite is set to 
\a no in \b solver.inp, separate files are written for solutions at each
output time.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The file \b Extras/ExtractSlice.c can be used to extract a slice perpendicular to any dimension
at a specified location along that dimension. The extracted slice is written out in the same
binary format as the original solutions files (with the same names op_xxxxx.bin) in a 
subdirectory called \b slices (\b Note: make the subdirectory called \a slices before running
this code). The codes Extras/BinaryToTecplot.c and Extras/BinaryToText.c
can then be used (in the \b slices subdirectory) to convert the binary slice solution file
to Tecplot or plain text files.
\b Note that it needs the relevant \b solver.inp that can be created as follows:
+ Copy the original solver.inp to the slices subdirectory.
+ In solver.inp, set \b ndims as \b 2, and remove the component of \b size and \b iproc 
  corresponding to the dimension being eliminated while extracting the slices (in this case,
  it is \a z or the 3rd component).

The following animation shows the vorticity magnitude and has been created by
plotting op_00301.bin through op_00400.bin, and shows the vortex shedding:
@image html Solution_3DNavStokCylinder_Shedding.gif

Expected screen output:
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible/output.log

  
\page ns3d_shock_cylinder_interaction Inviscid Shock-Cylinder Interaction 

Location: \b hypar/Examples/3D/NavierStokes3D/2D_Shock_Cylinder_Interaction

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h - by default
                     #NavierStokes3D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 3D Euler
                     equations are solved.)

Domain: \f$-2.5 \le x \le 7.5\f$, \f$-5 \le y \le 5\f$
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b hypar/Examples/STLGeometries/cylinder.stl)

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_ (with post-shock flow conditions)
  + xmax: Supersonic inflow #_SUPERSONIC_OUTFLOW_ (with pre-shock flow conditions)
  + ymin and ymax: Slip walls #_SLIP_WALL_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
+ O. Boiron, G. Chiavassa, R. Donat, A high-resolution penalization method for large Mach number 
  flows in the presence of obstacles, Computers & Fluids, 38 (2009), pp. 703-714, Section 5.1

Initial solution: see reference

Other parameters:
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/2D_Shock_Cylinder_Interaction/

\b solver.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction/aux/init.c

Output:
-------

Note that \b iproc is set to 

      2 2 1

in \b solver.inp (i.e., 2 processors along \a x, 2
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 4 MPI ranks (or change \b iproc).

After running the code, there should be 81 solution files \b op_00000.bin 
(initial solution), 
\b op_00001.bin, ..., \b op_00080.bin (solution at t=2). 
Since #HyPar::op_overwrite is set to 
\a no in \b solver.inp, separate files are written for solutions at each
output time.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The file \b Extras/ExtractSlice.c can be used to extract a slice perpendicular to any dimension
at a specified location along that dimension. The extracted slice is written out in the same
binary format as the original solutions files (with the same names op_xxxxx.bin) in a 
subdirectory called \b slices (\b Note: make the subdirectory called \a slices before running
this code). The codes Extras/BinaryToTecplot.c and Extras/BinaryToText.c
can then be used (in the \b slices subdirectory) to convert the binary slice solution file
to Tecplot or plain text files.
\b Note that it needs the relevant \b solver.inp that can be created as follows:
+ Copy the original solver.inp to the slices subdirectory.
+ In solver.inp, set \b ndims as \b 2, and remove the component of \b size and \b iproc 
  corresponding to the dimension being eliminated while extracting the slices (in this case,
  it is \a z or the 3rd component).

The following plot shows the density contours for the final solution (t=2):
@image html Solution_3DNavStokShockCyl_Density.png
and the following is the numerical Schlieren image (contour plot of \f$\|\nabla\rho\|_2\f$):
@image html Solution_3DNavStokShockCyl_Schlieren.png
@image html Solution_3DNavStokShockCyl_Schlieren2.png

Expected screen output:
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction/output.log

\page ns3d_sphere_steady_incompressible_viscous_adiabatic Steady, incompressible, viscous flow around an adiabatic sphere

Location: \b hypar/Examples/3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Adiabatic

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the sphere defined by [-2,6] X [-2,2] X [-2,2],
        and a stretched grid beyond this zone.

Geometry: A sphere of radius 0.5 centered at (0,0)
          (\b hypar/Examples/STLGeometries/sphere.stl)

The following image shows the sphere:
@image html Surface3D_Sphere.png

The following images shows the grid and the sphere:
@image html Domain3D_Sphere1.png
@image html Domain3D_Sphere2.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
  + Taneda, S., “Experimental Investigation of Wake behind a Sphere at Low Reynolds Numbers,” 
    Journal of the Physical Society of Japan, Vol. 11, 1956.
  + Johnson, T.A. and Patel, V.C., “Flow Past a Sphere up to a Reynolds Number of 300,” 
    Journal of Fluid Mechanics, Vol. 378, 1999.

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 100\f$ (#NavierStokes3D::Re) 
    (\b Note: since the diameter of the sphere is 1.0, the diameter-based Reynolds number 
    is the same as the specified Reynolds number \f$Re_D = Re = 100\f$).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Adiabatic/

\b solver.inp
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Adiabatic/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Adiabatic/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 100. To try other Reynolds numbers, change it here.
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Adiabatic/physics.inp

\b sphere.stl : the filename "sphere.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/sphere.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Adiabatic/aux/init.c

Output:
-------

Note that \b iproc is set to 

      8 4 4

in \b solver.inp (i.e., 8 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 128 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The following figure shows the flow (pressure and streamlines) at \f$Re_D=100\f$:
@image html Solution_3DNavStokSphere_ReD100.png

In addition to the main solution, the code also writes out a file with the aerodynamic
forces on the immersed body. This file is called \a surface.dat (if #HyPar::op_overwrite
is "yes") or \a surface_nnnnn.dat (if #HyPar::op_overwrite is "no", "nnnnn" is a numerical
index) (in this example, the file \b surface.dat is written out). This is an ASCII file in 
the Tecplot format, where the immersed body and the forces on it are represented using the 
"FETRIANGLE" type. The following image shows the surface pressure and temperature
on the sphere (front-view):
@image html IBSurface_3DNavStokSphereAdiabatic_Pressure.png
@image html IBSurface_3DNavStokSphereAdiabatic_Temperature.png

Expected screen output:
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Adiabatic/output.log

\page ns3d_sphere_steady_incompressible_viscous_isothermal Steady, incompressible, viscous flow around an isothermal sphere

Location: \b hypar/Examples/3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the sphere defined by [-2,6] X [-2,2] X [-2,2],
        and a stretched grid beyond this zone.

Geometry: A sphere of radius 0.5 centered at (0,0)
          (\b hypar/Examples/STLGeometries/sphere.stl)

The following image shows the sphere:
@image html Surface3D_Sphere.png

The following images shows the grid and the sphere:
@image html Domain3D_Sphere1.png
@image html Domain3D_Sphere2.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + The immersed body wall is specified as isothermal (#_IB_ISOTHERMAL_)
    with specified wall temperature \f$T_{\rm wall} = 1/\gamma\f$.

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 100\f$ (#NavierStokes3D::Re) 
    (\b Note: since the diameter of the sphere is 1.0, the diameter-based Reynolds number 
    is the same as the specified Reynolds number \f$Re_D = Re = 100\f$).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal/

\b solver.inp
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 100. To try other Reynolds numbers, change it here. This file
also specifies the sphere wall temperature as the same as the
freestream temperature.
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal/physics.inp

\b sphere.stl : the filename "sphere.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/sphere.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal/aux/init.c

Output:
-------

Note that \b iproc is set to 

      8 4 4

in \b solver.inp (i.e., 8 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 128 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

The following figure shows the flow (pressure and streamlines) at \f$Re_D=100\f$:
@image html Solution_3DNavStokSphereIsothermal_ReD100.png

In addition to the main solution, the code also writes out a file with the aerodynamic
forces on the immersed body. This file is called \a surface.dat (if #HyPar::op_overwrite
is "yes") or \a surface_nnnnn.dat (if #HyPar::op_overwrite is "no", "nnnnn" is a numerical
index) (in this example, the file \b surface.dat is written out). This is an ASCII file in 
the Tecplot format, where the immersed body and the forces on it are represented using the 
"FETRIANGLE" type. The following image shows the surface pressure and temperature
on the sphere (front-view):
@image html IBSurface_3DNavStokSphereIsothermal_Pressure.png
@image html IBSurface_3DNavStokSphereIsothermal_Temperature.png

Expected screen output:
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal/output.log

\page ns3d_sphere_steady_incompressible_viscous_isothermal_warm Steady, incompressible, viscous flow around an isothermal warm sphere

Location: \b hypar/Examples/3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal_WarmSphere

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the sphere defined by [-2,6] X [-2,2] X [-2,2],
        and a stretched grid beyond this zone.

Geometry: A sphere of radius 0.5 centered at (0,0)
          (\b hypar/Examples/STLGeometries/sphere.stl)

The following image shows the sphere:
@image html Surface3D_Sphere.png

The following images shows the grid and the sphere:
@image html Domain3D_Sphere1.png
@image html Domain3D_Sphere2.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + The immersed body wall is specified as isothermal (#_IB_ISOTHERMAL_)
    with specified wall temperature \f$T_{\rm wall} = 2/\gamma\f$ 
    (twice the freestream temperature).

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 100\f$ (#NavierStokes3D::Re) 
    (\b Note: since the diameter of the sphere is 1.0, the diameter-based Reynolds number 
    is the same as the specified Reynolds number \f$Re_D = Re = 100\f$).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal_WarmSphere/

\b solver.inp
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal_WarmSphere/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal_WarmSphere/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 100. To try other Reynolds numbers, change it here. This file
also specifies the sphere wall temperature that is twice the freestream
temperature.
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal_WarmSphere/physics.inp

\b sphere.stl : the filename "sphere.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/sphere.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal_WarmSphere/aux/init.c

Output:
-------

Note that \b iproc is set to 

      8 4 4

in \b solver.inp (i.e., 8 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 128 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

In addition to the main solution, the code also writes out a file with the aerodynamic
forces on the immersed body. This file is called \a surface.dat (if #HyPar::op_overwrite
is "yes") or \a surface_nnnnn.dat (if #HyPar::op_overwrite is "no", "nnnnn" is a numerical
index) (in this example, the file \b surface.dat is written out). This is an ASCII file in 
the Tecplot format, where the immersed body and the forces on it are represented using the 
"FETRIANGLE" type. 

The following figure shows the flow (pressure and velocity vectors) at \f$Re_D=100\f$ 
(A 2D x-y slice through the middle of the domain and the sphere surface are shown):
@image html Solution_3DNavStokSphereIsothermal_ReD100_WarmSphere.png

The following figure shows the temperature:
@image html Solution_3DNavStokSphereIsothermal_ReD100_T_WarmSphere.png

Expected screen output:
\include 3D/NavierStokes3D/Sphere/Steady_Viscous_Incompressible_Isothermal_WarmSphere/output.log

\page ns3d_sphere_unsteady_compressible_viscous_adiabatic Unsteady, compressible, viscous flow around an adiabatic sphere

Location: \b hypar/Examples/3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Adiabatic

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the sphere defined by [-2,6] X [-2,2] X [-2,2],
        and a stretched grid beyond this zone.

Geometry: A sphere of radius 0.5 centered at (0,0)
          (\b hypar/Examples/STLGeometries/sphere.stl)

The following image shows the sphere:
@image html Surface3D_Sphere.png

The following images shows the grid and the sphere:
@image html Domain3D_Sphere1.png
@image html Domain3D_Sphere2.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Initial solution: \f$\rho=1, u=0.8, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.8\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 1000000\f$ (#NavierStokes3D::Re) 
    (\b Note: since the diameter of the sphere is 1.0, the diameter-based Reynolds number 
    is the same as the specified Reynolds number \f$Re_D = Re = 1000000\f$).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Adiabatic/

\b solver.inp
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Adiabatic/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Adiabatic/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 1,000,000. To try other Reynolds numbers, change it here.
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Adiabatic/physics.inp

\b sphere.stl : the filename "sphere.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/sphere.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Adiabatic/aux/init.c

Output:
-------

Note that \b iproc is set to 

      8 4 4

in \b solver.inp (i.e., 8 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 128 MPI ranks (or change \b iproc).

After running the code, there should be 101 output files
\b op_nnnnn.bin, since #HyPar::op_overwrite is set to \a no in \b solver.inp.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

In addition to the main solution, the code also writes out a file with the aerodynamic
forces on the immersed body. This file is called \a surface_nnnnn.dat (if #HyPar::op_overwrite
is "no") or \a surface.dat (if #HyPar::op_overwrite is "yes") (in this example, the files 
\b surface_nnnnn.dat are written out). This is an ASCII file in 
the Tecplot format, where the immersed body and the forces on it are represented using the 
"FETRIANGLE" type. 

The following visualizations were generated using these solution files.
  + Pressure and velocity vector on a 2D x-y plane:
@image html Solution_3DNavStokSphere_ReD1mil_2Dflow.gif
  + Streamlines colored by velocity magnitude and sphere
    surface colored by temperature:
@image html Solution_3DNavStokSphere_ReD1mil_3Dflow.gif

Expected screen output:
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Adiabatic/output.log

\page ns3d_sphere_unsteady_compressible_viscous_isothermal Unsteady, compressible, viscous flow around an isothermal sphere

Location: \b hypar/Examples/3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Isothermal

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: The domain consists of a fine uniform grid around the sphere defined by [-2,6] X [-2,2] X [-2,2],
        and a stretched grid beyond this zone.

Geometry: A sphere of radius 0.5 centered at (0,0)
          (\b hypar/Examples/STLGeometries/sphere.stl)

The following image shows the sphere:
@image html Surface3D_Sphere.png

The following images shows the grid and the sphere:
@image html Domain3D_Sphere1.png
@image html Domain3D_Sphere2.png

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + The immersed body wall is specified as isothermal (#_IB_ISOTHERMAL_) 
    with surface temperature \f$T_{\rm wall}=1/\gamma\f$.

Initial solution: \f$\rho=1, u=0.8, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.8\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 1000000\f$ (#NavierStokes3D::Re) 
    (\b Note: since the diameter of the sphere is 1.0, the diameter-based Reynolds number 
    is the same as the specified Reynolds number \f$Re_D = Re = 1000000\f$).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Isothermal/

\b solver.inp
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Isothermal/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Isothermal/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 1,000,000. To try other Reynolds numbers, change it here.
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Isothermal/physics.inp

\b sphere.stl : the filename "sphere.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/sphere.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Isothermal/aux/init.c

Output:
-------

Note that \b iproc is set to 

      8 4 4

in \b solver.inp (i.e., 8 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 128 MPI ranks (or change \b iproc).

After running the code, there should be 101 output files
\b op_nnnnn.bin, since #HyPar::op_overwrite is set to \a no in \b solver.inp.
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following two codes
are available to convert the binary output file:
+ \b hypar/Extras/BinaryToTecplot.c - convert binary output file to 
  Tecplot file.
+ \b hypar/Extras/BinaryToText.c - convert binary output file to
  an ASCII text file (to visualize in, for example, MATLAB).

In addition to the main solution, the code also writes out a file with the aerodynamic
forces on the immersed body. This file is called \a surface_nnnnn.dat (if #HyPar::op_overwrite
is "no") or \a surface.dat (if #HyPar::op_overwrite is "yes") (in this example, the files 
\b surface_nnnnn.dat are written out). This is an ASCII file in 
the Tecplot format, where the immersed body and the forces on it are represented using the 
"FETRIANGLE" type. 

The following visualizations were generated using these solution files.
  + Pressure and velocity vector on a 2D x-y plane and sphere surface:
@image html Solution_3DNavStokSphere_ReD1mil_P_and_V.gif
  + Temperature on a 2D x-y plane and sphere surface:
@image html Solution_3DNavStokSphere_ReD1mil_T.gif

Expected screen output:
\include 3D/NavierStokes3D/Sphere/Unsteady_Viscous_Compressible_Isothermal/output.log

\page md_examples Multidomain Examples
HyPar has the capability to run ensemble simulations on multiple domains that 
of varying sizes and initial solution using the same discretization methods, 
time step sizes, final time, physics parameters, and boundary conditions. The following examples
demonstrate this functionality: 

\subpage md_linadv_1d_sine\n
\subpage md_euler2d_vortex

Other examples can also be run with this feature; the modifications needed are as follows:
  + Need to have an input file \b simulation.inp that specifies the number of
    simulations.
  + Modify the code to generate the intial solution (and exact solution, if available)
    to write out initial solution files for each simulation with the filenames 
    properly indexed, i.e., \b initial_<n>.inp, where n = 0, 1, ..., 
    number of simulations-1.
  + Modify \a size and \a iproc in \b solver.inp to specify the grid dimensions
    and MPI rank distribution for each simulation.
  + Note that multiple output files will be written, one set for each simulation.

\page md_linadv_1d_sine 1D Linear Advection - Sine Wave - Ensemble Simulation

Location: \b hypar/Examples/1D/LinearAdvection/SineWave_Ensemble
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Linear Advection Equation (linearadr.h)

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Ensemble Runs: 4 simulations on grids with 40, 80, 160, and 320 points.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b simulation.inp
\include 1D/LinearAdvection/SineWave_Ensemble/simulation.inp

\b solver.inp
\include 1D/LinearAdvection/SineWave_Ensemble/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/SineWave_Ensemble/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/SineWave_Ensemble/physics.inp

To generate \b initial_<n>.inp (initial solution) and
\b exact_<n>.inp (exact solution), compile and run the 
following code in the run directory. 
\include 1D/LinearAdvection/SineWave_Ensemble/aux/init.c
The index \a n is the simulation index (0, ..., 3 in 
this case).

Output:
-------
After running the code, there should be 24 output
files \b op_<n>_00000.dat, ... \b op_<n>_00005.dat, with
6 files for each of the 4 simulations in the set;
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. All the files are ASCII 
text (#HyPar::op_file_format is set to \a text in \b solver.inp).
In these files, the first column is grid index, the second column 
is x-coordinate, and the third column is the solution.

Solutions for each simulation: The following figure is obtained 
by plotting \a op_<n>_00005.dat (n=0,1,2,3): 
@image html Solution_1DLinearAdvSine_Ensemble.png

Since the exact solution is available at the final time,
the errors are calculated and reported on screen (see below)
as well as \b errors_<n>.dat. The errors for each of the 
simulation are as follows:
\include 1D/LinearAdvection/SineWave_Ensemble/errors_0.dat
\include 1D/LinearAdvection/SineWave_Ensemble/errors_1.dat
\include 1D/LinearAdvection/SineWave_Ensemble/errors_2.dat
\include 1D/LinearAdvection/SineWave_Ensemble/errors_3.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

\b Note: the errors above show 5th order convergence with grid size, 
as expected.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation_<n>.dat:
\include 1D/LinearAdvection/SineWave_Ensemble/conservation_0.dat
\include 1D/LinearAdvection/SineWave_Ensemble/conservation_1.dat
\include 1D/LinearAdvection/SineWave_Ensemble/conservation_2.dat
\include 1D/LinearAdvection/SineWave_Ensemble/conservation_3.dat
The numbers are: number of grid points (#HyPar::dim_global),
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 1D/LinearAdvection/SineWave_Ensemble/output.log

\page md_euler2d_vortex 2D Euler Equations - Isentropic Vortex Convection - Ensemble Simulation

Location: \b hypar/Examples/2D/NavierStokes2D/InviscidVortexConvection_Ensemble
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.1,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Ensemble Runs: 4 simulations on grids with (128,16), (64,32), (32,64), and (16,128) points.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b simulation.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/simulation.inp

\b solver.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/physics.inp

To generate \b initial_<n>.inp (initial solution) and \b exact_<n>.inp
(exact solution), compile and run the following code in the run 
directory.
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/aux/exact.c
The index \a n is the simulation index (0, ..., 3 in 
this case).

Output:
-------
Note that \b iproc for each simulation are specified as:

      8 1
      4 2
      2 4
      1 8

in \b solver.inp Thus, this example should be run
with 8 MPI ranks (or change \b iproc).

After running the code, there should be 44 output
files \b op_<n>_00000.dat, \b op_<n>_00001.dat, ... \b op_<n>_00010.dat,
with 11 files for each of the 4 simulations in the set; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=20\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, separate files are written
for solutions at each output time. 
  
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation shows the density contours as the vortex
convects over the domain, for each of the 4 simulations:
@image html Solution_md_2DNavStokVortex_0.gif
@image html Solution_md_2DNavStokVortex_1.gif
@image html Solution_md_2DNavStokVortex_2.gif
@image html Solution_md_2DNavStokVortex_3.gif

Since the exact solution is available at the final time, the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors_<n>.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/errors_0.dat
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/errors_1.dat
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/errors_2.dat
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/errors_3.dat
The numbers for each simulation are: number of grid points in each dimension 
(#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation error and prints it to screen, as well
as the file \b conservation_<n>.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/conservation_0.dat
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/conservation_1.dat
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/conservation_2.dat
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/conservation_3.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
and conservation error (#HyPar::ConservationError) of each component.

Expected screen output:
\include 2D/NavierStokes2D/InviscidVortexConvection_Ensemble/output.log

\page sg_examples Sparse Grids Examples
The following are some examples are simulated using the <B>sparse grids</B>
method. Familiarity with the sparse grids approach (specifically the combination
technique approach) is assumed.

\subpage sg_linear_adv_sinewave \n
\subpage sg_linear_adv_2d_sine_varyingadv

\subpage sg_euler2d_density_sinewave\n
\subpage sg_euler2d_vortconv\n
\subpage sg_euler2d_bubble

Any other simulation can also be run with the sparse grids method (as long as
the number of spatial dimensions is greater than 1); only the following is needed: 
  + Input file \b sparse_grids.inp that specifies parameters related to the sparse
    grids method.

As for the output, instead of \b op*.*, the full grid solution will be written
to \b op_fg*.*, and it will be in the same format/structure as the op*.* files.
If specified in the inputs, the sparse grid solution files will also be 
written to \b op_sg_<n>*.*.

\b Note: Sparse grids work well with smooth, linear simulations. Adapting them
for hyperbolic simulations with shocks and discontinuities is an area of active
research.

\page sg_linear_adv_sinewave 2D Linear Advection - Sine Wave

Location: \b hypar/Examples/2D/LinearAdvection/SineWave_SparseGrids
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Linear Advection Equation (linearadr.h)

Domain: \f$0 \le x,y < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \sin\left(2\pi x\right)\cos\left(2\pi y\right)\f$\n
Exact solution: \f$u\left(x,y,t\right) = u_0\left(x-a_xt,y-a_yt\right)\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b sparse_grids.inp
\include 2D/LinearAdvection/SineWave_SparseGrids/sparse_grids.inp

\b Note: The remaining files are the same as what would be 
required for a conventional (non-sparse-grids) simulation.

\b solver.inp
\include 2D/LinearAdvection/SineWave_SparseGrids/solver.inp

\b boundary.inp
\include 2D/LinearAdvection/SineWave_SparseGrids/boundary.inp

\b physics.inp (specifies \f$a_x\f$ and \f$a_y\f$)
\include 2D/LinearAdvection/SineWave_SparseGrids/physics.inp

To generate \b initial.inp and \b exact.inp, compile and run 
the following code in the run directory. (It is the same file 
as the one used in running a conventional non-sparse-grids 
simulation).
\include 2D/LinearAdvection/SineWave_SparseGrids/aux/exact.c

Output:
-------
Note that \b iproc does \b not need to be set for simulations
using sparse grids. HyPar will automatically calculate the load
balanced processor distribution for each sparse grid. If too
many processors are specified, then it will return an error.

After running the code, there should be the following output
files:

  + op_fg_00000.dat, ..., op_fg_00010.dat: these contain the full
    grid solution at \f$t=0, ..., 1\f$.
  + op_sg_<n>_00000.dat, ..., op_sg_<n>_00010.dat: these contain
    the solution on each of the sparse grids in the combination
    technique. These are written out because \b write_sg_solutions
    is set to \b yes in \b sparse_grids.inp 
    (SparseGridsSimulation::m_write_sg_solutions).
    
    
Since #HyPar::op_overwrite is set to \a no in \b solver.inp, 
separate files are written for solutions at each output time. 

#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation was generated from the full grid solution files
that were computed by using the combination technique on the sparse grids:
@image html Solution_SG_full_2DLinearAdvSineWave.gif

The following animations show the solution on some of the sparse grids
in the combination technique. The simulation was actually carried out on 
these grids. Note the different grid sizes.
@image html Solution_SG_sg_2DLinearAdvSineWave_00.gif
@image html Solution_SG_sg_2DLinearAdvSineWave_02.gif
@image html Solution_SG_sg_2DLinearAdvSineWave_04.gif
@image html Solution_SG_sg_2DLinearAdvSineWave_06.gif
@image html Solution_SG_sg_2DLinearAdvSineWave_08.gif
@image html Solution_SG_sg_2DLinearAdvSineWave_10.gif

Since the exact solution is available at the final time, the numerical 
errors are calculated for the recombined full grid solution and reported 
on screen (see below) as well as \b errors_fg.dat:
\include 2D/LinearAdvection/SineWave_SparseGrids/errors_fg.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since \b write_sg_errors is set to to \b yes in \b sparse_grids.inp,
(SparseGridsSimulation::m_print_sg_errors), the errors are also
calculated and reported (on screen and in the files \b errors_<n>.dat )
for each of the sparse grids. Following are the errors corresponding
to the sparse grids whose solutions are shown above.
\include 2D/LinearAdvection/SineWave_SparseGrids/errors_00.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/errors_02.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/errors_04.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/errors_06.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/errors_08.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/errors_10.dat

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation errors for each of the sparse grids
and prints it to screen, as well as the files \b conservation_<n>.dat:
\include 2D/LinearAdvection/SineWave_SparseGrids/conservation_00.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/conservation_02.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/conservation_04.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/conservation_06.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/conservation_08.dat
\include 2D/LinearAdvection/SineWave_SparseGrids/conservation_10.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt), and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 2D/LinearAdvection/SineWave_SparseGrids/output.log

\page sg_linear_adv_2d_sine_varyingadv 2D Linear Advection - Sine Wave with Spatially-Varying Advection Speed

Location: \b hypar/Examples/2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Linear Advection Equation (linearadr.h)

Domain: \f$0 \le x,y < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions on all boundaries.

Initial solution: \f$u\left(x,y,0\right) = u_0\left(x,y\right)= \cos\left(4\pi y\right)\f$\n

Advection speed: \f$a_x\left(x,y\right) = \sin\left(4\pi y\right)\f$, \f$a_y\left(x,y\right) = -\cos\left(4\pi x\right)\f$

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b sparse_grids.inp
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/sparse_grids.inp

\b Note: The remaining files are the same as what would be 
required for a conventional (non-sparse-grids) simulation.

\b solver.inp
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/solver.inp

\b boundary.inp
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/boundary.inp

\b physics.inp (specifies filename for advection field)
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/physics.inp
\b Note: Do not include the ".inp" extension in the filename above.

To generate \b initial.inp and \b advection.inp, compile and run the 
following code in the run directory. 
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/aux/init.c

Output:
-------
Note that \b iproc does \b not need to be set for simulations
using sparse grids. HyPar will automatically calculate the load
balanced processor distribution for each sparse grid. If too
many processors are specified, then it will return an error.

After running the code, there should be the following output
files:

  + op_fg_00000.dat, ..., op_fg_00020.dat: these contain the full
    grid solution at \f$t=0, ..., 4\f$.
  + op_sg_<n>_00000.dat, ..., op_sg_<n>_00020.dat: these contain
    the solution on each of the sparse grids in the combination
    technique. These are written out because \b write_sg_solutions
    is set to \b yes in \b sparse_grids.inp 
    (SparseGridsSimulation::m_write_sg_solutions).
    
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the final column is the 
solution.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation was generated from the full grid solution files
that were computed by using the combination technique on the sparse grids:
@image html Solution_SG_full_2DLinearAdvSine_VaryingAdv.gif

The following animations show the solution on some of the sparse grids
in the combination technique. The simulation was actually carried out on 
these grids. Note the different grid sizes.
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_0.gif
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_2.gif
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_4.gif
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_6.gif
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_8.gif

In addition to the usual output files, the linear advection physics 
module writes out the following files:
+ \b advection_field_<n>_00000.dat, ..., \b advection_field_<n>_00020.dat: 
  These files
  share the same format as the solution output files \b op_*.dat 
  and contains the advection field \f$a_x\left(x,y\right), a_y\left(x,y\right)\f$.

The following figure shows the vector plot of the velocity field on
some of the sparse grids:
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_0.png
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_2.png
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_4.png
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_6.png
@image html Solution_SG_sg_2DLinearAdvSine_VaryingAdv_8.png

Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation errors for each of the sparse grids
and prints it to screen, as well as the files \b conservation_<n>.dat:
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/conservation_0.dat
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/conservation_2.dat
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/conservation_4.dat
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/conservation_6.dat
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/conservation_8.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global),
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt), and conservation error (#HyPar::ConservationError).

Expected screen output:
\include 2D/LinearAdvection/SineWave_NonConstantAdvection_SparseGrids/output.log

\page sg_euler2d_density_sinewave 2D Euler Equations - Density Sine Wave

Location: \b hypar/Examples/2D/NavierStokes2D/DensitySineWave_SparseGrids
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Domain: \f$0 \le x,y \le 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 1,\ v_\infty = 1,\ p_\infty = \frac{1}{\gamma}
\f}
and the density wave is
\f{align}{
\rho = \rho_\infty + \frac{1}{10} \sin\left(2\pi x\right)\cos\left(2\pi y\right)
\f}

Numerical method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b sparse_grids.inp
\include 2D/NavierStokes2D/DensitySineWave_SparseGrids/sparse_grids.inp

\b Note: The remaining files are the same as what would be 
required for a conventional (non-sparse-grids) simulation.

\b solver.inp
\include 2D/NavierStokes2D/DensitySineWave_SparseGrids/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/DensitySineWave_SparseGrids/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/DensitySineWave_SparseGrids/physics.inp

\b lusolver.inp
\include 2D/NavierStokes2D/DensitySineWave_SparseGrids/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory. (It is the same file as the one used in running a
conventional non-sparse-grids simulation).
\include 2D/NavierStokes2D/DensitySineWave_SparseGrids/aux/exact.c

Output:
-------
Note that \b iproc does \b not need to be set for simulations
using sparse grids. HyPar will automatically calculate the load
balanced processor distribution for each sparse grid. If too
many processors are specified, then it will return an error.

After running the code, there should be the following output
files:
  + op_fg_00000.dat, ..., op_fg_00010.dat: these contain the full
    grid solution at \f$t=0, ..., 1\f$.

Since \b write_sg_solutions is set to \b no in \b sparse_grids.inp
(SparseGridsSimulation::m_write_sg_solutions), sparse grid solution
files are not written out.
    
Since #HyPar::op_overwrite is set to \a no in \b solver.inp, 
separate files are written for solutions at each output time. 
    
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation shows the density contours as the vortex
convects over the domain:, 
@image html Solution_SG_2DNavStokDensityWave.gif

Since the exact solution is available at the final time, the numerical 
errors are calculated for the recombined full grid solution and reported 
on screen (see below) as well as \b errors_fg.dat:
\include 2D/NavierStokes2D/DensitySineWave_SparseGrids/errors_fg.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since \b write_sg_errors is set to to \b no in \b sparse_grids.inp
(SparseGridsSimulation::m_print_sg_errors),
the errors for each of the sparse grids are not computed or reported.

Expected screen output:
\include 2D/NavierStokes2D/DensitySineWave_SparseGrids/output.log

\page sg_euler2d_vortconv 2D Euler Equations - Isentropic Vortex Convection

Location: \b hypar/Examples/2D/NavierStokes2D/InviscidVortexConvection_SparseGrids
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.1,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Input files required:
---------------------

\b sparse_grids.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_SparseGrids/sparse_grids.inp

\b Note: The remaining files are the same as what would be 
required for a conventional (non-sparse-grids) simulation.

\b solver.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_SparseGrids/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_SparseGrids/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_SparseGrids/physics.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory. (It is the same file as the one used in running a
conventional non-sparse-grids simulation).
\include 2D/NavierStokes2D/InviscidVortexConvection_SparseGrids/aux/exact.c

Output:
-------
Note that \b iproc does \b not need to be set for simulations
using sparse grids. HyPar will automatically calculate the load
balanced processor distribution for each sparse grid. If too
many processors are specified, then it will return an error.

After running the code, there should be the following output
files:

  + op_fg_00000.dat, ..., op_fg_00010.dat: these contain the full
    grid solution at \f$t=0, ..., 20\f$.
  + op_sg_<n>_00000.dat, ..., op_sg_<n>_00010.dat: these contain
    the solution on each of the sparse grids in the combination
    technique. These are written out because \b write_sg_solutions
    is set to \b yes in \b sparse_grids.inp 
    (SparseGridsSimulation::m_write_sg_solutions).
    
Since #HyPar::op_overwrite is set to \a no in \b solver.inp, 
separate files are written for solutions at each output time. 
    
#HyPar::op_file_format is set to \a tecplot2d in \b solver.inp, and
thus, all the files are in a format that Tecplot (http://www.tecplot.com/)
or other visualization software supporting the Tecplot format 
(e.g. VisIt - https://wci.llnl.gov/simulation/computer-codes/visit/)
can read. In these files, the first two lines are the Tecplot headers, 
after which the data is written out as: the first two columns are grid indices, 
the next two columns are x and y coordinates, and the remaining columns are the 
solution components.  #HyPar::op_file_format can be set to \a text to get the solution
files in plain text format (which can be read in and visualized in
MATLAB for example).

The following animation shows the density contours as the vortex
convects over the domain:, 
@image html Solution_SG_2DNavStokVortex.gif

The following animations show the solution on some of the sparse grids
in the combination technique. The simulation was actually carried out on 
these grids. Note the different grid sizes.
@image html Solution_SG_2DNavStokVortex_00.gif
@image html Solution_SG_2DNavStokVortex_02.gif
@image html Solution_SG_2DNavStokVortex_04.gif
@image html Solution_SG_2DNavStokVortex_06.gif
@image html Solution_SG_2DNavStokVortex_08.gif
@image html Solution_SG_2DNavStokVortex_10.gif

\b Note: Some of the coarse grid solutions look like nonsense, but they are 
still important constituents of the combination technique and contribute 
towards the full grid solution.

Since the exact solution is available at the final time, the numerical 
errors are calculated for the recombined full grid solution and reported 
on screen (see below) as well as \b errors_fg.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection_SparseGrids/errors_fg.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Since \b write_sg_errors is set to to \b no in \b sparse_grids.inp
(SparseGridsSimulation::m_print_sg_errors),
the errors for each of the sparse grids are not computed or reported.


Since #HyPar::ConservationCheck is set to \a yes in \b solver.inp,
the code checks for conservation errors for each of the sparse grids
and prints it to screen, but does not write it to files because 
\b write_sg_errors is set to to \b no in \b sparse_grids.inp.

Expected screen output:
\include 2D/NavierStokes2D/InviscidVortexConvection_SparseGrids/output.log

\page sg_euler2d_bubble 2D Euler Equations (with gravitational force) - Rising Thermal Bubble

Location: \b hypar/Examples/2D/NavierStokes2D/RisingThermalBubble_SparseGrids
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference:
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.2).

Domain: \f$0 \le x,y \le 1000\,m\f$, 
        "slip-wall" (#_SLIP_WALL_) boundary conditions on all sides.

Initial solution: See references above.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes2D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes2D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes2D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes2D::p0)
  + Hydrostatic balance type 2 (#NavierStokes2D::HB)

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Input files required:
---------------------

\b sparse_grids.inp
\include 2D/NavierStokes2D/RisingThermalBubble_SparseGrids/sparse_grids.inp

\b Note: The remaining files are the same as what would be 
required for a conventional (non-sparse-grids) simulation.

\b solver.inp
\include 2D/NavierStokes2D/RisingThermalBubble_SparseGrids/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/RisingThermalBubble_SparseGrids/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/RisingThermalBubble_SparseGrids/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/RisingThermalBubble_SparseGrids/aux/init.c

Output:
-------
Note that \b iproc does \b not need to be set for simulations
using sparse grids. HyPar will automatically calculate the load
balanced processor distribution for each sparse grid. If too
many processors are specified, then it will return an error.

After running the code, there should be the following output
files:
  + op_fg_00000.dat, ..., op_fg_00030.dat: these contain the full
    grid solution at \f$t=0, ..., 300\f$.

Since \b write_sg_solutions is set to \b no in \b sparse_grids.inp
(SparseGridsSimulation::m_write_sg_solutions), sparse grid solution
files are not written out.
    
#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/RisingThermalBubble_SparseGrids/aux/PostProcess.c

The following animation shows the potential temperature perturbation. 
It was plotted using VisIt
(https://wci.llnl.gov/simulation/computer-codes/visit/) with 
tecplot2d format chosen in the above postprocessing code.
@image html Solution_SG_2DNavStokRTB.gif

Expected screen output:
\include 2D/NavierStokes2D/RisingThermalBubble_SparseGrids/output.log


\page librom_examples libROM Examples

The following are some examples that use the reduced-order modeling capabilities
implemented in the \b libROM library (https://www.librom.net/). To run them,
HyPar needs to be compiled \b with \b libROM. Familiarity with using libROM and
the mathematical concepts implemented in it is assumed.

\b Note:
+ In general, any other example or simulation can use the capabilities of libROM
  by providing the libROM-specific input file (librom.inp).
+ If HyPar is not compiled with libROM, these examples will run as a regular 
  HyPar simulation.

These examples may come with a Python script to plot the solutions from the binary
solution file. This script will need the environment variable 
\b HYPAR_DIR (/path/to/hypar) to be defined, so make sure this variable exists.

Dynamic Mode Decomposition
==========================

Training Examples
-----------------

The following examples demonstrate training a #DMDROMObject using the solutions
from a PDE simulation. The same procedure (i.e. providing the libROM input file 
\b librom.inp with the appropriate inputs) can be used with any other simulation.\n
\n
\subpage linear_adv_sine_librom_dmd_train \n
\subpage linear_adv_disc_librom_dmd_train \n
\subpage sod_shock_tube_librom_dmd_train \n
\subpage shu_osher_librom_dmd_train \n
\n
\subpage euler2d_vortex_librom_dmd_train \n
\subpage euler2d_riemann4_librom_dmd_train \n
\subpage euler2d_rtb_librom_dmd_train \n
\subpage euler2d_igwave_petsc_librom_dmd_train (Requires PETSc)\n
\subpage navstok2d_ldsc_librom_dmd_train \n
\subpage vlasov_1d1v_twostreaminstability_librom_dmd_train \n
\n
\subpage ns3d_cylinder_steady_incompressible_viscous_librom_dmd_train \n
\subpage ns3d_cylinder_unsteady_incompressible_viscous_librom_dmd_train \n
\subpage ns3d_shock_cylinder_interaction_librom_dmd_train \n
\subpage ns3d_bubble_petsc_librom_dmd_train (Requires PETSc)\n

Prediction Examples
-------------------

The following examples demonstrate using a trained #DMDROMObject to compute the
solution instead of the discretizing the PDE in space and time. The same steps
(i.e. providing the libROM input file \b librom.inp and trained DMD objects
consistent with the computational domain) can be used for any other simulation.\n
\n
\subpage linear_adv_sine_librom_dmd_predict \n
\subpage linear_adv_disc_librom_dmd_predict \n
\subpage sod_shock_tube_librom_dmd_predict \n
\n
\subpage euler2d_vortex_librom_dmd_predict \n
\subpage navstok2d_ldsc_librom_dmd_predict \n
\subpage vlasov_1d1v_twostreaminstability_librom_dmd_predict \n
\n
\subpage ns3d_cylinder_steady_incompressible_viscous_librom_dmd_predict \n
\subpage ns3d_shock_cylinder_interaction_librom_dmd_predict \n

\page linear_adv_sine_librom_dmd_train 1D Linear Advection - Sine Wave

See \ref linear_adv_sine to familiarize yourself with this case.

Location: \b hypar/Examples/1D/LinearAdvection/SineWave_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Linear Advection Equation (linearadr.h)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Numerical method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: 3rd order SSP RK3 (TimeRK(), _RK_SSP3_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 2 (DMDROMObject::m_rdim)

Input files required:
---------------------

\b librom.inp
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/librom.inp

\b solver.inp
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/physics.inp

\b lusolver.inp (optional)
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/lusolver.inp

\b weno.inp (optional)
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/aux/init.c

Output:
-------
After running the code, there should be the following output
files:

+ 11 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
these are the \b HyPar solutions\b. 

+ 11 output files \b op_rom_00000.bin, \b op_rom_00001.bin, ... \b op_rom_00010.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=1\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a text, and GNUPlot 
or something similar can be used to plot the resulting text files.

The following plot shows the final solution - FOM (full-order model) refers to
the HyPar solution, ROM (reduced-order model) refers to the DMD solution.
@image html Solution_1DLinearAdvSinelibROM_DMD.png

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up), and total wall time.

By default, the code will write the trained DMD object(s) to files in a 
subdirectory (#DMDROMObject::m_dirname - default value is "DMD"). If the
subdirectory does not exist, the code \b may \b not report an error
(or give some HDF5 file-writing error); the DMD
objects will not be written! If the subdirectory exists, several files
will exist after the simulation is complete - they are in a format that
is readable by libROM.

Expected screen output:
\include 1D/LinearAdvection/SineWave_libROM_DMD_Train/out.log


\page linear_adv_disc_librom_dmd_train 1D Linear Advection - Discontinuous Waves (Time Windowed DMD)

See \ref linear_adv_disc to familiarize yourself with this case.

Location: \b hypar/Examples/1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Linear Advection Equation (linearadr.h)

Reduced Order Modeling: This example trains time-windowed DMD objects and 
then predicts the solution using the DMDs at the same times that the
actual HyPar solution is written at.

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$-1 \le x \le 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution:
  \f{equation}{
    u\left(x,0\right) = \left\{\begin{array}{lc} 
                          \exp\left(-\log\left(2\right)\frac{\left(x+7\right)^2}{0.0009}\right) & -0.8\le x \le -0.6 \\
                          1 & -0.4\le x \le -0.2 \\
                          1 - \left|10\left(x-0.1\right)\right| & 0\le x \le 0.2 \\
                          \sqrt{1-100\left(x-0.5\right)^2} & 0.4\le x \le 0.6 \\
                          0 & {\rm otherwise}
                        \end{array}\right.
  \f}

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Number of samples per time window: 100 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

\b librom.inp
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/librom.inp

\b solver.inp
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/physics.inp

\b lusolver.inp (optional)
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/lusolver.inp

\b weno.inp (optional)
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/aux/init.c

Output:
-------
After running the code, there should be the following output
files:

+ 51 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00050.bin; 
these are the \b HyPar solutions\b. 

+ 51 output files \b op_rom_00000.bin, \b op_rom_00001.bin, ... \b op_rom_00050.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=2\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a text, and GNUPlot 
or something similar can be used to plot the resulting text files.

The following plot shows the final solution - FOM (full-order model) refers to
the HyPar solution, ROM (reduced-order model) refers to the DMD solution.
@image html Solution_1DLinearAdvDisclibROM_DMD.png

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up), and total wall time.

By default, the code will write the trained DMD object(s) to files in a 
subdirectory (#DMDROMObject::m_dirname - default value is "DMD"). If the
subdirectory does not exist, the code \b may \b not report an error
(or give some HDF5 file-writing error); the DMD
objects will not be written! If the subdirectory exists, several files
will exist after the simulation is complete - they are in a format that
is readable by libROM.

Expected screen output:
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Train/out.log

\page sod_shock_tube_librom_dmd_train 1D Euler Equations - Sod Shock Tube (Time Windowed DMD)

See \ref sod_shock_tube to familiarize yourself with this case.

Description: 
-------------------

Location: \b hypar/Examples/1D/Euler1D/SodShockTube_libROM_DMD_Train 
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Euler equations (euler1d.h)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

References: 
  + G.A. Sod, "A survey of several finite difference methods 
    for systems of nonlinear hyperbolic conservation laws," 
    J. Comput. Phys., 27, 1 (1978).
  + C. B. Laney, "Computational Gasdynamics", Cambridge 
    University Press, 1998.

Domain: \f$0 \le x \le 1.0\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 1, u = 0, p = 1\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.125, u = 0, p = 0.1\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 8 (DMDROMObject::m_rdim)
 + Number of samples per time window: 40 (DMDROMObject::m_num_window_samples)

Input files required:
--------------------

\b librom.inp
\include 1D/Euler1D/SodShockTube_libROM_DMD_Train/librom.inp

\b solver.inp:
\include 1D/Euler1D/SodShockTube_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 1D/Euler1D/SodShockTube_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 1D/Euler1D/SodShockTube_libROM_DMD_Train/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/SodShockTube_libROM_DMD_Train/aux/init.c

Output:
-------
After running the code, there should be the following output
files:

+ 9 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00008.bin; 
these are the \b HyPar solutions\b. 

+ 9 output files \b op_rom_00000.bin, \b op_rom_00001.bin, ... \b op_rom_00008.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=2\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. It will plot the 3 conserved variables - density, momentum, and energy. 
Alternatively, #HyPar::op_file_format can be set to \a text, and GNUPlot 
or something similar can be used to plot the resulting text files.

The following plot shows the final solution (density) - FOM (full-order model) refers to
the HyPar solution, ROM (reduced-order model) refers to the DMD solution.
@image html Solution_1DSodShockTubeROM_DMD.png

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 1D/Euler1D/SodShockTube_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up), and total wall time.

By default, the code will write the trained DMD object(s) to files in a 
subdirectory (#DMDROMObject::m_dirname - default value is "DMD"). If the
subdirectory does not exist, the code \b may \b not report an error
(or give some HDF5 file-writing error); the DMD
objects will not be written! If the subdirectory exists, several files
will exist after the simulation is complete - they are in a format that
is readable by libROM.

Expected screen output:
\include 1D/Euler1D/SodShockTube_libROM_DMD_Train/out.log

\page euler2d_vortex_librom_dmd_train 2D Euler Equations - Isentropic Vortex Convection

See \ref euler2d_vortex to familiarize yourself with this case.

Location: \b hypar/Examples/2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.1,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)

Input files required:
---------------------

\b librom.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/librom.inp

\b solver.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/weno.inp

\b lusolver.inp (optional)
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/lusolver.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory.
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 21 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00020.bin; 
these are the \b HyPar solutions\b. 

+ 21 output files \b op_rom_00000.bin, \b op_rom_00001.bin, ... \b op_rom_00020.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=20\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a tecplot2d, and Tecplot/VisIt
or something similar can be used to plot the resulting files.

The following plot shows evolution of the density - FOM (full-order model) refers to
the HyPar solution, ROM (reduced-order model) refers to the DMD solution, and "Diff" is
their difference.
@image html Solution_2DNavStokVortex_libROM_DMD.gif

\b Wall \b clock \b times:
- PDE solution: 5.85 seconds
- DMD training time: 5.25 seconds
- DMD prediction/query time: 0.28 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

By default, the code will write the trained DMD object(s) to files in a 
subdirectory (#DMDROMObject::m_dirname - default value is "DMD"). If the
subdirectory does not exist, the code \b may \b not report an error
(or give some HDF5 file-writing error); the DMD
objects will not be written! If the subdirectory exists, several files
will exist after the simulation is complete - they are in a format that
is readable by libROM.

Expected screen output:
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Train/out.log

\page euler2d_riemann4_librom_dmd_train 2D Euler Equations - Riemann Problem Case 4 (Time Windowed DMD)

See \ref euler2d_riemann4 to familiarize yourself with this case.

Location: \b hypar/Examples/2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Reference:
  + P. Lax and X.-D. Liu, "Solution of two-dimensional Riemann
    problems of gas dynamics by positive schemes," SIAM J Sci 
    Comp 19 (1998), 319–340.

Domain: \f$-0.5 \le x,y \le 0.5\f$, \a "extrapolate" (#_EXTRAPOLATE_)
        boundary conditions.

Initial solution: see \b Case \b 4 in the reference.

Numerical method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Number of samples per time window: 50 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

\b librom.inp
\include 2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train/librom.inp

\b solver.inp
\include 2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory.
\include 2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 11 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
these are the \b HyPar solutions\b. 

+ 11 output files \b op_rom_00000.bin, \b op_rom_00001.bin, ... \b op_rom_00010.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=0.25\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a tecplot2d, and Tecplot/VisIt
or something similar can be used to plot the resulting text files.

The following plot shows the final solution (density) - FOM (full-order model) refers to
the HyPar solution, ROM (reduced-order model) refers to the DMD solution, and Diff
is the difference between the two.
@image html Solution_2DNavStokRiemann4_libROMDMD.png

\b Wall \b clock \b times:
- PDE solution: 111.1 seconds
- DMD training time: 17.6 seconds
- DMD prediction/query time: 1.4 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

By default, the code will write the trained DMD object(s) to files in a 
subdirectory (#DMDROMObject::m_dirname - default value is "DMD"). If the
subdirectory does not exist, the code \b may \b not report an error
(or give some HDF5 file-writing error); the DMD
objects will not be written! If the subdirectory exists, several files
will exist after the simulation is complete - they are in a format that
is readable by libROM.

Expected screen output:
\include 2D/NavierStokes2D/Riemann2DCase4_libROM_DMD_Train/out.log

\page vlasov_1d1v_twostreaminstability_librom_dmd_train 2D (1D-1V) Vlasov Equation - Two-Stream Instability (Time-Windowed DMD)

See \ref vlasov_1d1v_twostreaminstability to familiarize yourself with this case.

Location: \b hypar/Examples/2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D (1D-1V) Vlasov Equation (vlasov.h)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Domain: 
  + \f$0 \le x < 2\pi\f$, \a "periodic" (#_PERIODIC_)
  + \f$-7 \le y < 7\f$, \a "dirichlet" (#_DIRICHLET_) (\f$f = 0\f$)


Initial solution: 
\f$f\left(x,v\right) = \frac{4}{\pi T}\left(1+\frac{1}{10}\cos\left(2k\pi\frac{x}{L}\right)\right)\left[\exp\left(-\frac{\left(v-2\right)^2}{2T}\right) + \exp\left(-\frac{\left(v+2\right)^2}{2T}\right)\right]\f$, \f$k=1,T=1,L=2\pi\f$.

Self-consistent electric field is computed by solving the Poisson equation
in a periodic domain using Fourier transforms. This examples *requires* HyPar
to be compiled with FFTW (http://www.fftw.org/).

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Number of samples per time window: 200 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

\b librom.inp
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Train/librom.inp

\b solver.inp
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Train/physics.inp

To generate \b initial.inp (initial solution), compile and run the 
following code in the run directory. 
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Train/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 26 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00025.bin; 
these are the \b HyPar solutions\b. 

+ 26 output files \b op_rom_00000.bin, \b op_rom_00001.bin, ... \b op_rom_00025.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=5\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a tecplot2d or \a text, 
and Tecplot/VisIt or something similar can be used to plot the resulting files.

The following plot shows the evolution of the distribution function - 
FOM (full-order model) refers to the HyPar solution, 
ROM (reduced-order model) refers to the DMD solution, and 
Diff is the difference between the two.
@image html Solution_1D1VVlasov_TwoStreamInstability_libROM_DMD.gif

\b Wall \b clock \b times:
- PDE solution: 11.34 seconds
- DMD training time: 2.30 seconds
- DMD prediction/query time: 0.34 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

By default, the code will write the trained DMD object(s) to files in a 
subdirectory (#DMDROMObject::m_dirname - default value is "DMD"). If the
subdirectory does not exist, the code \b may \b not report an error
(or give some HDF5 file-writing error); the DMD
objects will not be written! If the subdirectory exists, several files
will exist after the simulation is complete - they are in a format that
is readable by libROM.

Expected screen output:
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Train/out.log

\page navstok2d_ldsc_librom_dmd_train 2D Navier-Stokes Equations -  Lid-Driven Square Cavity (Time-Windowed DMD)

See \ref navstok2d_ldsc to familiarize yourself with this case.

Location: \b hypar/Examples/2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Navier-Stokes Equations (navierstokes2d.h)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Reference: 
+ Erturk, E., Corke, T.C., and Gokcol, C., "Numerical Solutions of
  2-D Steady Incompressible Driven Cavity Flow at High Reynolds Numbers",
  International Journal for Numerical Methods in Fluids, 48, 2005,
  http://dx.doi.org/10.1002/fld.953.
+ Ghia, U., Ghia, K.N., Shin, C.T., "High-Re Solutions for Incompressible
  Flow using the Navier-Stokes Equations and a Multigrid Method", Journal
  of Computational Physics, 48, 1982, http://dx.doi.org/10.1016/0021-9991(82)90058-4.

Note that this is an incompressible problem being solved here using the compressible
Navier-Stokes equations in terms of non-dimensional flow variables. The density and 
pressure are taken such that the speed of sound is 1.0, and the flow velocities 
specified in the initial and boundary conditions correspond to a characteristic 
Mach number of 0.1 (thus, they are 0.1 times the values in the above reference).

Domain: \f$0 \le x, y \le 1\f$

Boundary conditions:
+ No-slip wall BC on \f$x=0,1, 0 \le y \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ Moving no-slip wall BC on \f$y=1, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with specified 
  wall velocity of 0.1 in the x-direction).

Initial solution: \f$\rho=1, p=1/\gamma\f$. The velocities are specified according to 
the references above, but scaled by a factor of 0.1 to ensure that the characteristic
Mach number is 0.1.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 100, 1000, 3200\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the cavity length and width.
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.1\f$ (characteristic Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: 4th order, 4-stage Runge-Kutta (TimeRK(), #_RK_44_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Sampling frequency: 50 (libROMInterface::m_sampling_freq)
 + Number of samples per time window: 200 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

\b librom.inp
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Train/librom.inp

\b solver.inp
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Train/boundary.inp

\b physics.inp (\b Note: this file specifies \f$Re = 3200\f$,
change \a Re here for other Reynolds numbers.)
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Train/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Train/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processor along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 1 output file \b op.bin; 
this is the \b HyPar solutions\b. 

+ 1 output file \b op_rom.bin; 
this is the \b predicted solutions from the DMD object(s)\b.

All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotStreamlines.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a tecplot2d, and Tecplot/VisIt
or something similar can be used to plot the resulting text files.

The following plot shows the final solution (velocity streamlines colored by
the velocity magnitude) - FOM (full-order model) refers to
the HyPar solution, ROM (reduced-order model) refers to the DMD solution, and Diff
is the difference in the velocity magnitude between the two.
@image html Solution_2DNavStokLDSC_Re3200_libROM_DMD.png

\b Wall \b clock \b times:
- PDE solution: 554.6 seconds
- DMD training time: 58.6 seconds
- DMD prediction/query time: 0.3 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

By default, the code will write the trained DMD object(s) to files in a 
subdirectory (#DMDROMObject::m_dirname - default value is "DMD"). If the
subdirectory does not exist, the code \b may \b not report an error
(or give some HDF5 file-writing error); the DMD
objects will not be written! If the subdirectory exists, several files
will exist after the simulation is complete - they are in a format that
is readable by libROM.

Expected screen output:
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Train/out.log

\page ns3d_cylinder_steady_incompressible_viscous_librom_dmd_train 3D Navier-Stokes - Steady, incompressible, viscous flow around a cylinder (Time-Windowed DMD)

See \ref ns3d_cylinder_steady_incompressible_viscous to familiarize yourself with this case.

Location: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Domain: The domain consists of a fine uniform grid around the cylinder defined by [-2,6] X [-2,2],
        and a stretched grid beyond this zone.
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b hypar/Examples/STLGeometries/cylinder.stl)

See \ref ns3d_cylinder_steady_incompressible_viscous for pictures showing the domain and geometry.

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
  + Taneda, S., "Experimental Investigation of the Wakes behind Cylinders and Plates at Low 
    Reynolds Numbers," Journal of the Physical Society of Japan, Vol. 11, 302–307, 1956. 
  + Dennis, S. C. R., Chang, G.-Z., "Numerical solutions for steady flow past a circular
    cylinder at Reynolds numbers up to 100", Journal of Fluid Mechanics, 42 (3), 1970,
    pp. 471-489.

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 10\f$ (#NavierStokes3D::Re) (\b Note: 
    since the diameter of the cylinder is 2.0, the cylinder-diameter-based Reynolds number is 
    \f$Re_D = 2Re = 20\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Sampling frequency: 10 (libROMInterface::m_sampling_freq)
 + Number of samples per time window: 500 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train

\b librom.inp
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train/librom.inp

\b solver.inp
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 10 (corresponding to \f$Re_D\f$ of 20). To try other Reynolds 
numbers, change it here.
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train/aux/init.c

Output:
-------

Note that \b iproc is set to 

      8 4 1

in \b solver.inp (i.e., 8 processors along \a x, 4
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 32 MPI ranks (or change \b iproc).

Please see the original example in the "Immersed Boundaries Examples" section
for a full description of the output files and how to view them. The following only
contains libROM-specific comments.

After running the code, there should be the following output
files:

+ 1 output file \b op.bin; 
this is the \b HyPar solutions\b. 

+ 1 output file \b op_rom.bin; 
this is the \b predicted solutions from the DMD object(s)\b.

All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a tecplot3d, and Tecplot/VisIt
or something similar can be used to plot the resulting text files.

The following figure shows the flow at \f$Re_D=20\f$. The pressure is plotted in the overall
domain, and the wake is shown by plotting the x-velocity, where it is negative.
FOM (full-order model) refers to the HyPar solution, ROM (reduced-order model) refers to 
the DMD solution, and Diff is the difference between the two.
@image html Solution_3DNavStokCylinder_ReD020_libROM_DMD.png

\b Wall \b clock \b times:
- PDE solution: 3038 seconds
- DMD training time: 249 seconds
- DMD prediction/query time: 1.3 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Train/out.log

\page ns3d_cylinder_unsteady_incompressible_viscous_librom_dmd_train 3D Navier-Stokes - Unsteady, incompressible, viscous flow around a cylinder with vortex shedding (Time-Windowed DMD)

See \ref ns3d_cylinder_unsteady_incompressible_viscous to familiarize yourself with this case.

Location: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Domain: The domain consists of a fine uniform grid around the cylinder defined by [-4,12] X [-2,2],
        and a stretched grid beyond this zone.
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b hypar/Examples/STLGeometries/cylinder.stl)

See \ref ns3d_cylinder_unsteady_incompressible_viscous for pictures show the domain and geometry.

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
  + Taneda, S., Experimental Investigation of the Wakes behind Cylinders and Plates at Low 
    Reynolds Numbers, Journal of the Physical Society of Japan, Vol. 11, 302–307, 1956. 

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 50\f$ (#NavierStokes3D::Re) (\b Note: 
    since the diameter of the cylinder is 2.0, the cylinder-diameter-based Reynolds number is 
    \f$Re_D = 2Re = 100\f$.

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (NavierStokes3DParabolicFunction())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Sampling frequency: 2 (libROMInterface::m_sampling_freq)
 + Number of samples per time window: 1000 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train/

\b librom.inp
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train/librom.inp

\b solver.inp
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), the following code will generate the
domain with a clustered mesh and uniform horizontal flow. It will take
time to develop the time-periodic vortex shedding flow from this initial solution.
The following steps are recommended:
- Use the code below to generate a initial solution file (initial.inp) in a separate directory,
  and set the \a n_iter to 50000 in that directory. Change \a op_overwrite to "yes" to avoid
  writing too many solution files. Run HyPar.
- Use the code hypar/Extras/BinaryOPToInitialSolution.c to generate
  the initial solution file from the final solution of the above simulation.
  This code will write out a file called "solution.inp"; just rename it
  to "initial.inp" and place it in the directory where this simulation will
  be run.
- Now, using this initial solution, where the time-periodic vortex shedding has developed,
  run this simulation with frequent solution output to visualize the vortex shedding.

\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train/aux/init.c

Output:
-------

Note that \b iproc is set to 

      8 4 1

in \b solver.inp (i.e., 8 processors along \a x, 4
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 32 MPI ranks (or change \b iproc).

Please see the original example in the "Immersed Boundaries Examples" section
for a full description of the output files and how to view them. The following only
contains libROM-specific comments.

After running the code, there should be the following output
files:

+ 41 output file \b op_00000.bin, ..., op_00040.bin; 
these are the \b HyPar solutions\b. 

+ 41 output file \b op_rom_00000.bin, ..., op_rom_00040.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a tecplot3d, and Tecplot/VisIt
or something similar can be used to plot the resulting text files.

The following figure shows the vortex shedding. The pressure is plotted in the overall
domain, and the wake is shown by plotting the x-velocity, where it is negative.
FOM (full-order model) refers to the HyPar solution, ROM (reduced-order model) refers to 
the DMD solution, and Diff is the difference between the two.
@image html Solution_3DNavStokCylinder_ReD100_libROM_DMD.gif

\b Wall \b clock \b times:
- PDE solution: 2413 seconds
- DMD training time: 2113 seconds
- DMD prediction/query time: 9.8 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 3D/NavierStokes3D/2D_Cylinder/Unsteady_Viscous_Incompressible_libROM_DMD_Train/out.log

\page euler2d_rtb_librom_dmd_train 2D Euler Equations (with gravitational force) - Rising Thermal Bubble (Time-Windowed DMD)

See \ref euler2d_rtb to familiarize yourself with this case.

Location: \b hypar/Examples/2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Reference:
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.2).

Domain: \f$0 \le x,y \le 1000\,m\f$, 
        "slip-wall" (#_SLIP_WALL_) boundary conditions on all sides.

Initial solution: See references above.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes2D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes2D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes2D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes2D::p0)
  + Hydrostatic balance type 2 (#NavierStokes2D::HB)

Numerical method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSPRK3 (TimeRK(), #_RK_SSP3_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Sampling frequency: 7 (libROMInterface::m_sampling_freq)
 + Number of samples per time window: 1000 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

\b librom.inp
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/librom.inp

\b solver.inp
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/physics.inp

\b weno.inp (optional)
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/aux/init.c

Output:
-------
Note that \b iproc is set to 

      6 6

in \b solver.inp (i.e., 6 processors along \a x, and 6
processor along \a y). Thus, this example should be run
with 36 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 41 output file \b op_00000.bin, ..., \b op_00040.bin; 
this is the \b HyPar solutions\b. 

+ 41 output file \b op_rom_00000.bin, ..., \b op_rom_00040.bin; 
this is the \b predicted solutions from the DMD object(s)\b.

All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/aux/PostProcess.cpp

The following animation shows the potential temperature for the HyPar solutions.
It was plotted using VisIt
(https://wci.llnl.gov/simulation/computer-codes/visit/) with 
tecplot2d format chosen in the above postprocessing code.
@image html Solution_2DNavStokRTB_libROM_DMD_FOM.gif
The following animation shows the same for the libROM solution predicted
by the DMD object:
@image html Solution_2DNavStokRTB_libROM_DMD_ROM.gif

The provided Python script (\b plotDiff.py) computes and plots the diff
of the HyPar and DMD solutions for each of the conserved variables. The
following figure shows this diff at the final time for the density.
@image html Solution_2DNavStokRTB_libROM_DMD_diff_density.png

\b Wall \b clock \b times:
- PDE solution: 1419 seconds
- DMD training time: 1068 seconds
- DMD prediction/query time: 2.73 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 2D/NavierStokes2D/RisingThermalBubble_libROM_DMD_Train/out.log

\page shu_osher_librom_dmd_train 1D Euler Equations - Shu-Osher Problem (Time Windowed DMD)

See \ref shu_osher to familiarize yourself with this case.

Description: 
------------

Location: \b hypar/Examples/1D/Euler1D/ShuOsherProblem_libROM_DMD_Train 
          (This directory contains all the input files needed
          to run this case. If there is a \a Run.m, run it in
          MATLAB to quickly set up, run, and visualize the 
          example).

Governing equations: 1D Euler equations (euler1d.h)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

References: 
  + C.-W. Shu and S. Osher, "Efficient implementation of
    essentially non-oscillatory schemes ,II," J. Comput.
    Phys., 83 (1989), pp. 32–78

Domain: \f$-5 \le x \le 5\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ -5 \le x < -4\f$: \f$\rho = 27/7, u = 4\sqrt{35}/7, p = 31/3\f$
  + \f$ -4 \le x \le 5\f$: \f$\rho = 1 + 0.2\sin\left(5x\right), u = 0, p = 1\f$

Numerical Method:
 + Spatial discretization (hyperbolic): Characteristic-based 5th order WENO (Interp1PrimFifthOrderWENOChar())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Number of samples per time window: 60 (DMDROMObject::m_num_window_samples)

Input files required:
--------------------
\b librom.inp:
\include 1D/Euler1D/ShuOsherProblem_libROM_DMD_Train/librom.inp

\b solver.inp:
\include 1D/Euler1D/ShuOsherProblem_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 1D/Euler1D/ShuOsherProblem_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 1D/Euler1D/ShuOsherProblem_libROM_DMD_Train/physics.inp

\b weno.inp (optional)
\include 1D/Euler1D/ShuOsherProblem_libROM_DMD_Train/weno.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/ShuOsherProblem_libROM_DMD_Train/aux/init.c

Output:
-------
After running the code, there should be the following output
files:

+ 11 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
these are the \b HyPar solutions\b. 

+ 11 output files \b op_rom_00000.bin, \b op_rom_00001.bin, ... \b op_rom_00010.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=1.8\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. It will plot the 3 conserved variables - density, momentum, and energy. 
Alternatively, #HyPar::op_file_format can be set to \a text, and GNUPlot 
or something similar can be used to plot the resulting text files.

The following plot shows the final solution (density) - FOM (full-order model) refers to
the HyPar solution, ROM (reduced-order model) refers to the DMD solution.
@image html Solution_1DShuOsherProblem_libROM_DMD.png

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 1D/Euler1D/ShuOsherProblem_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points (#HyPar::dim_global), 
number of processors (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up), and total wall time.

By default, the code will write the trained DMD object(s) to files in a 
subdirectory (#DMDROMObject::m_dirname - default value is "DMD"). If the
subdirectory does not exist, the code \b may \b not report an error
(or give some HDF5 file-writing error); the DMD
objects will not be written! If the subdirectory exists, several files
will exist after the simulation is complete - they are in a format that
is readable by libROM.

Expected screen output:
\include 1D/Euler1D/ShuOsherProblem_libROM_DMD_Train/out.log

\page ns3d_shock_cylinder_interaction_librom_dmd_train 3D Navier-Stokes - Inviscid Shock-Cylinder Interaction (Time-Windowed DMD)

See \ref ns3d_shock_cylinder_interaction to familiarize yourself with this case.

Location: \b hypar/Examples/3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h - by default
                     #NavierStokes3D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 3D Euler
                     equations are solved.)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Domain: \f$-2.5 \le x \le 7.5\f$, \f$-5 \le y \le 5\f$
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b hypar/Examples/STLGeometries/cylinder.stl)

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_ (with post-shock flow conditions)
  + xmax: Supersonic inflow #_SUPERSONIC_OUTFLOW_ (with pre-shock flow conditions)
  + ymin and ymax: Slip walls #_SLIP_WALL_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
+ O. Boiron, G. Chiavassa, R. Donat, A high-resolution penalization method for large Mach number 
  flows in the presence of obstacles, Computers & Fluids, 38 (2009), pp. 703-714, Section 5.1

Initial solution: see reference

Other parameters:
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: RK4 (TimeRK(), #_RK_44_)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 32 (DMDROMObject::m_rdim)
 + Number of samples per time window: 100 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

These files are all located in: \b hypar/Examples/3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train/

\b librom.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train/librom.inp

\b solver.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train/aux/init.c

Output:
-------

Note that \b iproc is set to 

      4 4 1

in \b solver.inp (i.e., 4 processors along \a x, 4
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 16 MPI ranks (or change \b iproc).

Please see the original example in the "Immersed Boundaries Examples" section
for a full description of the output files and how to view them. The following only
contains libROM-specific comments.

After running the code, there should be the following output
files:

+ 11 output file \b op_00000.bin, ..., op_00010.bin; 
these are the \b HyPar solutions\b. 

+ 11 output file \b op_rom_00000.bin, ..., op_rom_00010.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a tecplot3d, and Tecplot/VisIt
or something similar can be used to plot the resulting text files.

The following plot shows the density contours for the final solution (t=2).
FOM (full-order model) refers to the HyPar solution, ROM (reduced-order model) refers to 
the DMD solution, and Diff is the difference between the two.
@image html Solution_3DNavStokShockCyl_Density_libROM_DMD.png

\b Wall \b clock \b times:
- PDE solution: 815 seconds
- DMD training time: 155 seconds
- DMD prediction/query time: 33 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Train/out.log

\page linear_adv_sine_librom_dmd_predict 1D Linear Advection - Sine Wave

See \ref linear_adv_sine to familiarize yourself with this case. This example uses 
a DMD object that has already been trained (see \ref linear_adv_sine_librom_dmd_train).

Location: \b hypar/Examples/1D/LinearAdvection/SineWave_libROM_DMD_Predict
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Linear Advection Equation (linearadr.h)

Reduced Order Modeling: This example predicts the solution from a trained 
DMD object. The code does \b not solve the PDE by discretizing in space and
integrating in time.

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$0 \le x < 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution: \f$u\left(x,0\right) = \sin\left(2\pi x\right)\f$

Reduced Order Modeling:
 + Mode: predict (libROMInterface::m_mode, #_ROM_MODE_PREDICT_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)

Note:
-----

In this mode, HyPar will run just like an usual PDE simulation, except that it will swap out
the numerical spatial discretization and time integration with the ROM-based prediction. The
input files and output files will be the same as a regular simulation with the following 
comments:
+ Numerical method inputs are ignored (eg. those that specify spatial discretization and
  time integration methods).
+ In solver.inp, the values for dt, n_iter, and file_op_iter is used only to compute the 
  simulation times at which to compute and write the solution. The time step size, dt, 
  need not respect any CFL criterion.
+ #HyPar::ConservationCheck is set to "no" since it is not possible to compute conservation
  loss for a general domain (because boundary fluxes are not being computed).

Input files required:
---------------------

\b librom.inp
\include 1D/LinearAdvection/SineWave_libROM_DMD_Predict/librom.inp

<B> DMD Object(s) </B>:\n
The trained DMD object(s) must be located in the directory specified in \b librom.inp
as \a dmd_dirname (#DMDROMObject::m_dirname). For this example, they were generated
using \ref linear_adv_sine_librom_dmd_train.

\b solver.inp
\include 1D/LinearAdvection/SineWave_libROM_DMD_Predict/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/SineWave_libROM_DMD_Predict/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/SineWave_libROM_DMD_Predict/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/SineWave_libROM_DMD_Predict/aux/init.c

Output:
-------
After running the code, there should be the following output
files:

+ 31 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00030.bin; 
these are the solutions as <B> predicted by the ROM </B>. 

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=3\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b Examples/Python/plotSolution_1DBinary.py)
can be used to generate plots from these binary files.
Alternatively, #HyPar::op_file_format can be set to \a text, and GNUPlot 
or something similar can be used to plot the resulting text files.

The animation shows the evolution of the solution.
@image html Solution_1DLinearAdvSinelibROM_DMD_Predict.gif

Expected screen output:
\include 1D/LinearAdvection/SineWave_libROM_DMD_Predict/out.log


\page linear_adv_disc_librom_dmd_predict 1D Linear Advection - Discontinuous Waves (Time Windowed DMD)

See \ref linear_adv_disc to familiarize yourself with this case. This example uses 
a DMD object that has already been trained (see \ref linear_adv_disc_librom_dmd_train).

Location: \b hypar/Examples/1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Predict
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Linear Advection Equation (linearadr.h)

Reduced Order Modeling: This example predicts the solution from trained time-windowed 
DMD objects. The code does \b not solve the PDE by discretizing in space and
integrating in time.

References:
  + Ghosh, D., Baeder, J. D., "Compact Reconstruction Schemes with 
    Weighted ENO Limiting for Hyperbolic Conservation Laws", 
    SIAM Journal on Scientific Computing, 34 (3), 2012, A1678–A1706

Domain: \f$-1 \le x \le 1\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions

Initial solution:
  \f{equation}{
    u\left(x,0\right) = \left\{\begin{array}{lc} 
                          \exp\left(-\log\left(2\right)\frac{\left(x+7\right)^2}{0.0009}\right) & -0.8\le x \le -0.6 \\
                          1 & -0.4\le x \le -0.2 \\
                          1 - \left|10\left(x-0.1\right)\right| & 0\le x \le 0.2 \\
                          \sqrt{1-100\left(x-0.5\right)^2} & 0.4\le x \le 0.6 \\
                          0 & {\rm otherwise}
                        \end{array}\right.
  \f}

Reduced Order Modeling:
 + Mode: predict (libROMInterface::m_mode, #_ROM_MODE_PREDICT_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)

Note:
-----

In this mode, HyPar will run just like an usual PDE simulation, except that it will swap out
the numerical spatial discretization and time integration with the ROM-based prediction. The
input files and output files will be the same as a regular simulation with the following 
comments:
+ Numerical method inputs are ignored (eg. those that specify spatial discretization and
  time integration methods).
+ In solver.inp, the values for dt, n_iter, and file_op_iter is used only to compute the 
  simulation times at which to compute and write the solution. The time step size, dt, 
  need not respect any CFL criterion.
+ #HyPar::ConservationCheck is set to "no" since it is not possible to compute conservation
  loss for a general domain (because boundary fluxes are not being computed).

Input files required:
---------------------

\b librom.inp
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Predict/librom.inp

<B> DMD Object(s) </B>:\n
The trained DMD object(s) must be located in the directory specified in \b librom.inp
as \a dmd_dirname (#DMDROMObject::m_dirname). For this example, they were generated
using \ref linear_adv_disc_librom_dmd_train.

\b solver.inp
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Predict/solver.inp

\b boundary.inp
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Predict/boundary.inp

\b physics.inp
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Predict/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory. \b Note: if the
final time is an integer multiple of the time period,
the file \b initial.inp can also be used as the exact
solution \b exact.inp (i.e. create a sym link called 
\a exact.inp pointing to \a initial.inp, or just copy
\a initial.inp to \a exact.inp).
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Predict/aux/init.c

Output:
-------
After running the code, there should be the following output
files:

+ 21 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00020.bin; 
these are the solutions as <B> predicted by the ROM </B>. 

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=3\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b Examples/Python/plotSolution_1DBinary.py)
can be used to generate plots from these binary files.
Alternatively, #HyPar::op_file_format can be set to \a text, and GNUPlot 
or something similar can be used to plot the resulting text files.

The animation shows the evolution of the solution.
@image html Solution_1DLinearAdvDisclibROM_DMD_Predict.gif

Expected screen output:
\include 1D/LinearAdvection/DiscontinuousWaves_libROM_DMD_Predict/out.log

\page sod_shock_tube_librom_dmd_predict 1D Euler Equations - Sod Shock Tube (Time Windowed DMD)

See \ref sod_shock_tube to familiarize yourself with this case. This example
uses a DMD object that has already been trained (see \ref sod_shock_tube_librom_dmd_train).

Location: \b hypar/Examples/1D/Euler1D/SodShockTube_libROM_DMD_Predict 
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 1D Euler equations (euler1d.h)

Reduced Order Modeling: This example predicts the solution from a trained 
DMD object. The code does \b not solve the PDE by discretizing in space and
integrating in time.

References: 
  + G.A. Sod, "A survey of several finite difference methods 
    for systems of nonlinear hyperbolic conservation laws," 
    J. Comput. Phys., 27, 1 (1978).
  + C. B. Laney, "Computational Gasdynamics", Cambridge 
    University Press, 1998.

Domain: \f$0 \le x \le 1.0\f$, \a "extrapolate" (#_EXTRAPOLATE_) 
        boundary conditions

Initial Solution:
  + \f$ 0 \le x < 0.5\f$: \f$\rho = 1, u = 0, p = 1\f$
  + \f$ 0.5 \le x \le 1\f$: \f$\rho = 0.125, u = 0, p = 0.1\f$

Reduced Order Modeling:
 + Mode: predict (libROMInterface::m_mode, #_ROM_MODE_PREDICT_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)

Note:
-----

In this mode, HyPar will run just like an usual PDE simulation, except that it will swap out
the numerical spatial discretization and time integration with the ROM-based prediction. The
input files and output files will be the same as a regular simulation with the following 
comments:
+ Numerical method inputs are ignored (eg. those that specify spatial discretization and
  time integration methods).
+ In solver.inp, the values for dt, n_iter, and file_op_iter is used only to compute the 
  simulation times at which to compute and write the solution. The time step size, dt, 
  need not respect any CFL criterion.
+ #HyPar::ConservationCheck is set to "no" since it is not possible to compute conservation
  loss for a general domain (because boundary fluxes are not being computed).

Input files required:
--------------------

\b librom.inp
\include 1D/Euler1D/SodShockTube_libROM_DMD_Predict/librom.inp

<B> DMD Object(s) </B>:\n
The trained DMD object(s) must be located in the directory specified in \b librom.inp
as \a dmd_dirname (#DMDROMObject::m_dirname). For this example, they were generated
using \ref sod_shock_tube_librom_dmd_train.

\b solver.inp:
\include 1D/Euler1D/SodShockTube_libROM_DMD_Predict/solver.inp

\b boundary.inp
\include 1D/Euler1D/SodShockTube_libROM_DMD_Predict/boundary.inp

\b physics.inp
\include 1D/Euler1D/SodShockTube_libROM_DMD_Predict/physics.inp

To generate \b initial.inp, compile and run the 
following code in the run directory:
\include 1D/Euler1D/SodShockTube_libROM_DMD_Predict/aux/init.c

Output:
-------
After running the code, there should be the following output
files:

+ 21 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00020.bin; 
these are the solutions as <B> predicted by the ROM </B>. 

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=0.2\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b Examples/Python/plotSolution_1DBinary.py)
can be used to generate plots from these binary files.
Alternatively, #HyPar::op_file_format can be set to \a text, and GNUPlot 
or something similar can be used to plot the resulting text files.

The animation shows the evolution of the solution (density).
@image html Solution_1DSodShockTubeROM_DMD_Predict.gif

Expected screen output:
\include 1D/Euler1D/SodShockTube_libROM_DMD_Predict/out.log

\page euler2d_vortex_librom_dmd_predict 2D Euler Equations - Isentropic Vortex Convection

See \ref euler2d_vortex to familiarize yourself with this case.
This example uses a DMD object that has already been trained 
(see \ref euler2d_vortex_librom_dmd_train).

Location: \b hypar/Examples/2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Predict
          (This directory contains all the input files needed
          to run this case.)

Reduced Order Modeling: This example predicts the solution from a trained 
DMD object. The code does \b not solve the PDE by discretizing in space and
integrating in time.

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reference: C.-W. Shu, "Essentially Non-oscillatory and Weighted Essentially 
           Non-oscillatory Schemes for Hyperbolic Conservation Laws", 
           ICASE Report 97-65, 1997

Domain: \f$0 \le x,y \le 10\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions.

Initial solution: The freestream flow is given by
\f{equation}{
  \rho_\infty = 1,\ u_\infty = 0.1,\ v_\infty = 0,\ p_\infty = 1
\f}
and a vortex is introduced, specified as
\f{align}{
\rho &= \left[ 1 - \frac{\left(\gamma-1\right)b^2}{8\gamma\pi^2} e^{1-r^2} \right]^{\frac{1}{\gamma-1}},\ p = \rho^\gamma, \\
u &= u_\infty - \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(y-y_c\right),\ v = v_\infty + \frac{b}{2\pi} e^{\frac{1}{2}\left(1-r^2\right)} \left(x-x_c\right),
\f}
where \f$b=0.5\f$ is the vortex strength and \f$r = \left[(x-x_c)^2 + (y-y_c)^2 \right]^{1/2}\f$ is the distance from the vortex center \f$\left(x_c,y_c\right) = \left(5,5\right)\f$.

Reduced Order Modeling:
 + Mode: predict (libROMInterface::m_mode, #_ROM_MODE_PREDICT_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)

Note:
-----

In this mode, HyPar will run just like an usual PDE simulation, except that it will swap out
the numerical spatial discretization and time integration with the ROM-based prediction. The
input files and output files will be the same as a regular simulation with the following 
comments:
+ Numerical method inputs are ignored (eg. those that specify spatial discretization and
  time integration methods).
+ In solver.inp, the values for dt, n_iter, and file_op_iter is used only to compute the 
  simulation times at which to compute and write the solution. The time step size, dt, 
  need not respect any CFL criterion.
+ #HyPar::ConservationCheck is set to "no" since it is not possible to compute conservation
  loss for a general domain (because boundary fluxes are not being computed).

Input files required:
---------------------

\b librom.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Predict/librom.inp

<B> DMD Object(s) </B>:\n
The trained DMD object(s) must be located in the directory specified in \b librom.inp
as \a dmd_dirname (#DMDROMObject::m_dirname). For this example, they were generated
using \ref euler2d_vortex_librom_dmd_train.

\b solver.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Predict/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Predict/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Predict/physics.inp

To generate \b initial.inp (initial solution) and \b exact.inp
(exact solution), compile and run the following code in the run 
directory.
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Predict/aux/exact.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 41 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00040.bin; 
these are the solutions as <B> predicted by the ROM </B>. 

The first of each of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=40\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b Examples/Python/plotSolution_2DBinary.py)
can be used to generate plots from these binary files.
Alternatively, #HyPar::op_file_format can be set to \a text or \a tecplot2d, and 
Tecplot/VisIt or something similar can be used to plot the resulting files.

The animation shows the evolution of the solution (density).
@image html Solution_2DNavStokVortex_libROM_DMD_Predict.gif

Since the exact solution is available at the final time (\a exact.inp), the numerical 
errors are calculated and reported on screen (see below)
as well as \b errors.dat:
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Predict/errors.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity errors (#HyPar::error),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 2D/NavierStokes2D/InviscidVortexConvection_libROM_DMD_Predict/out.log

\page navstok2d_ldsc_librom_dmd_predict 2D Navier-Stokes Equations -  Lid-Driven Square Cavity (Time-Windowed DMD)

See \ref navstok2d_ldsc to familiarize yourself with this case.
This example uses a DMD object that has already been trained 
(see \ref navstok2d_ldsc_librom_dmd_train).

Location: \b hypar/Examples/2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Predict
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Navier-Stokes Equations (navierstokes2d.h)

Reduced Order Modeling: This example predicts the solution from a trained 
DMD object. The code does \b not solve the PDE by discretizing in space and
integrating in time.

Reference: 
+ Erturk, E., Corke, T.C., and Gokcol, C., Numerical Solutions of
  2-D Steady Incompressible Driven Cavity Flow at High Reynolds Numbers,
  International Journal for Numerical Methods in Fluids, 48, 2005,
  http://dx.doi.org/10.1002/fld.953.
+ Ghia, U., Ghia, K.N., Shin, C.T., High-Re Solutions for Incompressible
  Flow using the Navier-Stokes Equations and a Multigrid Method, Journal
  of Computational Physics, 48, 1982, http://dx.doi.org/10.1016/0021-9991(82)90058-4.

Note that this is an incompressible problem being solved here using the compressible
Navier-Stokes equations in terms of non-dimensional flow variables. The density and 
pressure are taken such that the speed of sound is 1.0, and the flow velocities 
specified in the initial and boundary conditions correspond to a characteristic 
Mach number of 0.1 (thus, they are 0.1 times the values in the above reference).

Domain: \f$0 \le x, y \le 1\f$

Boundary conditions:
+ No-slip wall BC on \f$x=0,1, 0 \le y \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ Moving no-slip wall BC on \f$y=1, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with specified 
  wall velocity of 0.1 in the x-direction).

Initial solution: \f$\rho=1, p=1/\gamma\f$. The velocities are specified according to 
the references above, but scaled by a factor of 0.1 to ensure that the characteristic
Mach number is 0.1.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 100, 1000, 3200\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the cavity length and width.
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.1\f$ (characteristic Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Reduced Order Modeling:
 + Mode: predict (libROMInterface::m_mode, #_ROM_MODE_PREDICT_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)

Note:
-----

In this mode, HyPar will run just like an usual PDE simulation, except that it will swap out
the numerical spatial discretization and time integration with the ROM-based prediction. The
input files and output files will be the same as a regular simulation with the following 
comments:
+ Numerical method inputs are ignored (eg. those that specify spatial discretization and
  time integration methods).
+ In solver.inp, the values for dt, n_iter, and file_op_iter is used only to compute the 
  simulation times at which to compute and write the solution. The time step size, dt, 
  need not respect any CFL criterion.
+ #HyPar::ConservationCheck is set to "no" since it is not possible to compute conservation
  loss for a general domain (because boundary fluxes are not being computed).


Input files required:
---------------------

\b librom.inp
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Predict/librom.inp

<B> DMD Object(s) </B>:\n
The trained DMD object(s) must be located in the directory specified in \b librom.inp
as \a dmd_dirname (#DMDROMObject::m_dirname). For this example, they were generated
using \ref navstok2d_ldsc_librom_dmd_train.

\b solver.inp
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Predict/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Predict/boundary.inp

\b physics.inp (\b Note: this file specifies \f$Re = 3200\f$,
change \a Re here for other Reynolds numbers.)
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Predict/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Predict/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processor along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 1 output file \b op.bin; 
this is the solutions as <B> predicted by the ROM </B>. 

It is a binary file.
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b Examples/Python/plotStreamlines_2DBinary.py)
can be used to generate streamline plots from the binary file.
Alternatively, #HyPar::op_file_format can be set to \a tecplot2d, and Tecplot/VisIt
or something similar can be used to plot the resulting files.

The following plot shows the final solution (velocity streamlines colored by
the velocity magnitude):
@image html Solution_2DNavStokLDSC_Re3200_libROM_DMD_Predict.png

Expected screen output:
\include 2D/NavierStokes2D/LidDrivenCavity_libROM_DMD_Predict/out.log

\page vlasov_1d1v_twostreaminstability_librom_dmd_predict 2D (1D-1V) Vlasov Equation - Two-Stream Instability (Time-Windowed DMD)

See \ref vlasov_1d1v_twostreaminstability to familiarize yourself with this case.
This example uses a DMD object that has already been trained 
(see \ref vlasov_1d1v_twostreaminstability_librom_dmd_train).

Location: \b hypar/Examples/2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Predict
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D (1D-1V) Vlasov Equation (vlasov.h)

Reduced Order Modeling: This example predicts the solution from a trained 
DMD object. The code does \b not solve the PDE by discretizing in space and
integrating in time.

Domain: 
  + \f$0 \le x < 2\pi\f$, \a "periodic" (#_PERIODIC_)
  + \f$-7 \le y < 7\f$, \a "dirichlet" (#_DIRICHLET_) (\f$f = 0\f$)

Initial solution: 
\f$f\left(x,v\right) = \frac{4}{\pi T}\left(1+\frac{1}{10}\cos\left(2k\pi\frac{x}{L}\right)\right)\left[\exp\left(-\frac{\left(v-2\right)^2}{2T}\right) + \exp\left(-\frac{\left(v+2\right)^2}{2T}\right)\right]\f$, \f$k=1,T=1,L=2\pi\f$.

Self-consistent electric field is computed by solving the Poisson equation
in a periodic domain using Fourier transforms. This examples *requires* HyPar
to be compiled with FFTW (http://www.fftw.org/).

Reduced Order Modeling:
 + Mode: predict (libROMInterface::m_mode, #_ROM_MODE_PREDICT_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)

Note:
-----

In this mode, HyPar will run just like an usual PDE simulation, except that it will swap out
the numerical spatial discretization and time integration with the ROM-based prediction. The
input files and output files will be the same as a regular simulation with the following 
comments:
+ Numerical method inputs are ignored (eg. those that specify spatial discretization and
  time integration methods).
+ In solver.inp, the values for dt, n_iter, and file_op_iter is used only to compute the 
  simulation times at which to compute and write the solution. The time step size, dt, 
  need not respect any CFL criterion.
+ #HyPar::ConservationCheck is set to "no" since it is not possible to compute conservation
  loss for a general domain (because boundary fluxes are not being computed).

Input files required:
---------------------

\b librom.inp
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Predict/librom.inp

<B> DMD Object(s) </B>:\n
The trained DMD object(s) must be located in the directory specified in \b librom.inp
as \a dmd_dirname (#DMDROMObject::m_dirname). For this example, they were generated
using \ref vlasov_1d1v_twostreaminstability_librom_dmd_train.

\b solver.inp
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Predict/solver.inp

\b boundary.inp
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Predict/boundary.inp

\b physics.inp
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Predict/physics.inp

To generate \b initial.inp (initial solution), compile and run the 
following code in the run directory. 
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Predict/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processors along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 26 output files \b op_00000.bin, \b op_00001.bin, ... \b op_00025.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

The first of these file sets is the solution at \f$t=0\f$ and 
the final one is the solution at \f$t=5\f$. Since #HyPar::op_overwrite is
set to \a no in \b solver.inp, a separate file is written
for solutions at each output time. All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b Examples/Python/plotSolution_2DBinary.py)
can be used to generate plots from the binary files that compare the HyPar and DMD
solutions. Alternatively, #HyPar::op_file_format can be set to \a tecplot2d or \a text, 
and Tecplot/VisIt or something similar can be used to plot the resulting files.

The following plot shows the evolution of the distribution function:
@image html Solution_1D1VVlasov_TwoStreamInstability_libROM_DMD_Predict.gif

Expected screen output:
\include 2D/Vlasov1D1V/TwoStreamInstability_libROM_DMD_Predict/out.log

\page ns3d_cylinder_steady_incompressible_viscous_librom_dmd_predict 3D Navier-Stokes - Steady, incompressible, viscous flow around a cylinder (Time-Windowed DMD)

See \ref ns3d_cylinder_steady_incompressible_viscous to familiarize yourself with this case.
This example uses a DMD object that has already been trained 
(see \ref ns3d_cylinder_steady_incompressible_viscous_librom_dmd_train).

Location: \b hypar/Examples/3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Predict

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Reduced Order Modeling: This example predicts the solution from a trained 
DMD object. The code does \b not solve the PDE by discretizing in space and
integrating in time.

Domain: The domain consists of a fine uniform grid around the cylinder defined by [-2,6] X [-2,2],
        and a stretched grid beyond this zone.
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b hypar/Examples/STLGeometries/cylinder.stl)

See \ref ns3d_cylinder_steady_incompressible_viscous for pictures showing the domain and geometry.

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_
  + xmax: Subsonic outflow #_SUBSONIC_OUTFLOW_
  + ymin and ymax: Subsonic "ambivalent" #_SUBSONIC_AMBIVALENT_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
  + Taneda, S., "Experimental Investigation of the Wakes behind Cylinders and Plates at Low 
    Reynolds Numbers," Journal of the Physical Society of Japan, Vol. 11, 302–307, 1956. 
  + Dennis, S. C. R., Chang, G.-Z., "Numerical solutions for steady flow past a circular
    cylinder at Reynolds numbers up to 100", Journal of Fluid Mechanics, 42 (3), 1970,
    pp. 471-489.

Initial solution: \f$\rho=1, u=0.1, v=w=0, p=1/\gamma\f$ everywhere in the domain.

Other parameters:
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Freestream Mach number \f$M_{\infty} = 0.1\f$ (#NavierStokes3D::Minf)
  + Prandlt number \f$Pr = 0.72\f$ (#NavierStokes3D::Pr)
  + Reynolds number \f$Re = \frac {\rho u L } {\mu} = 10\f$ (#NavierStokes3D::Re) (\b Note: 
    since the diameter of the cylinder is 2.0, the cylinder-diameter-based Reynolds number is 
    \f$Re_D = 2Re = 20\f$.

Reduced Order Modeling:
 + Mode: predict (libROMInterface::m_mode, #_ROM_MODE_PREDICT_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)

Note:
-----

In this mode, HyPar will run just like an usual PDE simulation, except that it will swap out
the numerical spatial discretization and time integration with the ROM-based prediction. The
input files and output files will be the same as a regular simulation with the following 
comments:
+ Numerical method inputs are ignored (eg. those that specify spatial discretization and
  time integration methods).
+ In solver.inp, the values for dt, n_iter, and file_op_iter is used only to compute the 
  simulation times at which to compute and write the solution. The time step size, dt, 
  need not respect any CFL criterion.
+ #HyPar::ConservationCheck is set to "no" since it is not possible to compute conservation
  loss for a general domain (because boundary fluxes are not being computed).

Input files required:
---------------------

<B> DMD Object(s) </B>:\n
The trained DMD object(s) must be located in the directory specified in \b librom.inp
as \a dmd_dirname (#DMDROMObject::m_dirname). For this example, they were generated
using \ref ns3d_cylinder_steady_incompressible_viscous_librom_dmd_train.

\b librom.inp
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Predict/librom.inp

\b solver.inp
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Predict/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Predict/boundary.inp

\b physics.inp : The following file specifies a Reynolds number
of 10 (corresponding to \f$Re_D\f$ of 20). To try other Reynolds 
numbers, change it here.
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Predict/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Predict/aux/init.c

Output:
-------

Note that \b iproc is set to 

      8 4 1

in \b solver.inp (i.e., 8 processors along \a x, 4
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 32 MPI ranks (or change \b iproc).

Please see the original example in the "Immersed Boundaries Examples" section
for a full description of the output files and how to view them. The following only
contains libROM-specific comments.

After running the code, there should be the following output
files:

+ 1 output file \b op.bin; 
this is the \b predicted solutions from the DMD object(s)\b.

The files is in binary format.
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files.
Alternatively, #HyPar::op_file_format can be set to \a tecplot3d, and Tecplot/VisIt
or something similar can be used to plot the resulting text files.

The following figure shows the flow at \f$Re_D=20\f$. The pressure is plotted in the overall
domain, and the wake is shown by plotting the x-velocity, where it is negative.
@image html Solution_3DNavStokCylinder_ReD020_libROM_DMD_Predict.png

Expected screen output:
\include 3D/NavierStokes3D/2D_Cylinder/Steady_Viscous_Incompressible_libROM_DMD_Predict/out.log

\page ns3d_shock_cylinder_interaction_librom_dmd_predict 3D Navier-Stokes - Inviscid Shock-Cylinder Interaction (Time-Windowed DMD)

See \ref ns3d_shock_cylinder_interaction to familiarize yourself with this case.
This example uses a DMD object that has already been trained 
(see \ref ns3d_shock_cylinder_interaction_librom_dmd_train).

Location: \b hypar/Examples/3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Predict

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h - by default
                     #NavierStokes3D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 3D Euler
                     equations are solved.)

Reduced Order Modeling: This example predicts the solution from a trained 
DMD object. The code does \b not solve the PDE by discretizing in space and
integrating in time.

Domain: \f$-2.5 \le x \le 7.5\f$, \f$-5 \le y \le 5\f$
\b Note: This is a 2D flow simulated using a 3-dimensional setup by taking the length of the
         domain along \a z to be very small and with only 3 grid points (the domain size along \a z
         \b must \b be smaller than the cylinder length).

Geometry: A cylinder of radius 1.0 centered at (0,0)
          (\b hypar/Examples/STLGeometries/cylinder.stl)

Boundary conditions:
  + xmin: Subsonic inflow #_SUBSONIC_INFLOW_ (with post-shock flow conditions)
  + xmax: Supersonic inflow #_SUPERSONIC_OUTFLOW_ (with pre-shock flow conditions)
  + ymin and ymax: Slip walls #_SLIP_WALL_
  + zmin and zmax: Periodic #_PERIODIC_ (to simulate a 2D flow in the x-y plane)
  + The immersed body wall is specified as adiabatic (#_IB_ADIABATIC_);
    (this is the default).

Reference:
+ O. Boiron, G. Chiavassa, R. Donat, A high-resolution penalization method for large Mach number 
  flows in the presence of obstacles, Computers & Fluids, 38 (2009), pp. 703-714, Section 5.1

Initial solution: see reference

Other parameters:
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)

Reduced Order Modeling:
 + Mode: predict (libROMInterface::m_mode, #_ROM_MODE_PREDICT_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)

Note:
-----

In this mode, HyPar will run just like an usual PDE simulation, except that it will swap out
the numerical spatial discretization and time integration with the ROM-based prediction. The
input files and output files will be the same as a regular simulation with the following 
comments:
+ Numerical method inputs are ignored (eg. those that specify spatial discretization and
  time integration methods).
+ In solver.inp, the values for dt, n_iter, and file_op_iter is used only to compute the 
  simulation times at which to compute and write the solution. The time step size, dt, 
  need not respect any CFL criterion.
+ #HyPar::ConservationCheck is set to "no" since it is not possible to compute conservation
  loss for a general domain (because boundary fluxes are not being computed).

Input files required:
---------------------

<B> DMD Object(s) </B>:\n
The trained DMD object(s) must be located in the directory specified in \b librom.inp
as \a dmd_dirname (#DMDROMObject::m_dirname). For this example, they were generated
using \ref ns3d_shock_cylinder_interaction_librom_dmd_train.

\b librom.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Predict/librom.inp

\b solver.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Predict/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Predict/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Predict/physics.inp

\b cylinder.stl : the filename "cylinder.stl" \b must match
the input for \a immersed_body in \a solver.inp.\n
Located at \b hypar/Examples/STLGeometries/cylinder.stl

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Predict/aux/init.c

Output:
-------

Note that \b iproc is set to 

      4 4 1

in \b solver.inp (i.e., 4 processors along \a x, 4
processors along \a y, and 1 processor along \a z). Thus, 
this example should be run with 16 MPI ranks (or change \b iproc).

Please see the original example in the "Immersed Boundaries Examples" section
for a full description of the output files and how to view them. The following only
contains libROM-specific comments.

After running the code, there should be the following output
files:

+ 11 output file \b op_00000.bin, ..., op_00010.bin; 
these are the \b predicted solutions from the DMD object(s)\b.

All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files.
Alternatively, #HyPar::op_file_format can be set to \a tecplot3d, and Tecplot/VisIt
or something similar can be used to plot the resulting text files.

The following plot shows the density contours for the final solution (t=2).
@image html Solution_3DNavStokShockCyl_Density_libROM_DMD_Predict.png

Expected screen output:
\include 3D/NavierStokes3D/2D_Shock_Cylinder_Interaction_libROM_DMD_Predict/out.log

\page ns2d_ldsc_petsc_implicit 2D Navier-Stokes Equations -  Lid-Driven Square Cavity

Location: \b hypar/Examples/2D/NavierStokes2D/LidDrivenCavity_PETSc_Implicit
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Navier-Stokes Equations (navierstokes2d.h)

Reference: 
+ Erturk, E., Corke, T.C., and Gokcol, C., "Numerical Solutions of
  2-D Steady Incompressible Driven Cavity Flow at High Reynolds Numbers",
  International Journal for Numerical Methods in Fluids, 48, 2005,
  http://dx.doi.org/10.1002/fld.953.
+ Ghia, U., Ghia, K.N., Shin, C.T., "High-Re Solutions for Incompressible
  Flow using the Navier-Stokes Equations and a Multigrid Method", Journal
  of Computational Physics, 48, 1982, http://dx.doi.org/10.1016/0021-9991(82)90058-4.

Note that this is an incompressible problem being solved here using the compressible
Navier-Stokes equations in terms of non-dimensional flow variables. The density and 
pressure are taken such that the speed of sound is 1.0, and the flow velocities 
specified in the initial and boundary conditions correspond to a characteristic 
Mach number of 0.1 (thus, they are 0.1 times the values in the above reference).

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, but with the entire RHS on the implicit side
(by setting the \b -ts_arkimex_fully_implicit flag in PETSc inputs).
So, the time integration is implicit (a 6-stage, 4-th order DIRK).

Domain: \f$0 \le x, y \le 1\f$

Boundary conditions:
+ No-slip wall BC on \f$x=0,1, 0 \le y \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ No-slip wall BC on \f$y=0, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with 0 wall velocity).
+ Moving no-slip wall BC on \f$y=1, 0 \le x \le 1\f$ (#_NOSLIP_WALL_ with specified 
  wall velocity of 0.1 in the x-direction).

Initial solution: \f$\rho=1, p=1/\gamma\f$. The velocities are specified according to 
the references above, but scaled by a factor of 0.1 to ensure that the characteristic
Mach number is 0.1.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 3200\f$ (Reynolds number) (#NavierStokes2D::Re), 
    where \f$L=1\f$ is the cavity length and width.
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes2D::Pr)
  + \f$M_\infty = 0.1\f$ (characteristic Mach number) (#NavierStokes2D::Minf)

\b Note: Pressure is taken as \f$1/\gamma\f$ in the above so that the freestream 
speed of sound is 1.

Numerical method:
 + Spatial discretization (hyperbolic): 5th order upwind (Interp1PrimFifthOrderUpwind())
 + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                        non-conservative 2-stage (ParabolicFunctionNC2Stage())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (fully-implicit mode) (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: implicit part of <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX4.html)

Input files required:
---------------------

<B>.petscrc</B>
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_Implicit/petscrc

\b solver.inp
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_Implicit/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_Implicit/boundary.inp

\b physics.inp (\b Note: this file specifies \f$Re = 3200\f$,
change \a Re here for other Reynolds numbers.)
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_Implicit/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_Implicit/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4

in \b solver.inp (i.e., 4 processors along \a x, and 4
processor along \a y). Thus, this example should be run
with 16 MPI ranks (or change \b iproc).

After running the code, there should be one output file
\b op.bin, since #HyPar::op_overwrite is set to \a yes in \b solver.inp.
Since #HyPar::op_file_format is set to \a binary in \b solver.inp.

The provided Python script (\b Examples/Python/plotStreamlines_2DBinary.py)
can be used to generate streamline plots from the binary file.
Alternatively, #HyPar::op_file_format can be set to \a tecplot2d, and Tecplot/VisIt
or something similar can be used to plot the resulting files.

The following plot shows the final solution (velocity streamlines colored by
the velocity magnitude):
@image html Solution_2DNavStokLDSC_Re3200_PETSc_Implicit.png

The file <B>function_counts.dat</B> reports the computational expense
(in terms of the number of function counts):
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_Implicit/function_counts.dat
The numbers are, respectively,
+ Time iterations
+ Number of times the hyperbolic term was evaluated
+ Number of times the parabolic term was evaluated
+ Number of times the source term was evaluated
+ Number of calls to the explicit right-hand-side function (PetscRHSFunctionIMEX() or PetscRHSFunctionExpl())
+ Number of calls to the implicit right-hand-side function (PetscIFunctionIMEX() or PetscRHSFunctionImpl())
+ Number of calls to the Jacobian (PetscIJacobianIMEX() or PetscIJacobian())
+ Number of calls to the matrix-free Jacobian function (PetscJacobianFunctionIMEX_Linear(), PetscJacobianFunctionIMEX_JFNK(), PetscJacobianFunction_JFNK(), or PetscJacobianFunction_Linear()).

Expected screen output (for Reynolds number 3200):
\include 2D/NavierStokes2D/LidDrivenCavity_PETSc_Implicit/out.log

\page euler2d_igwave_petsc_librom_dmd_train 2D Euler Equations (with gravitational force) - Inertia-Gravity Waves (Time-Windowed DMD) with PETSc IMEX time integration

See \ref euler2d_igwave_petsc to familiarize yourself with this case.

Location: \b hypar/Examples/2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 2D Euler Equations (navierstokes2d.h - By default,
                     #NavierStokes2D::Re is set to \b -1 which makes the
                     code skip the parabolic terms, i.e., the 2D Euler
                     equations are solved.)

Reduced Order Modeling: This example trains a DMD object and then
predicts the solution using the DMD at the same times that the
actual HyPar solution is written at.

Reference:
  + W. C. Skamarock and J. B. Klemp, "Efficiency and accuracy of 
    the Klemp-Wilhelmson timesplitting technique", Monthly Weather 
    Review, 122 (1994), pp. 2623–2630.
  + Giraldo, F.X., Restelli, M., "A study of spectral element and
    discontinuous Galerkin methods for the Navier–Stokes equations
    in nonhydrostatic mesoscale atmospheric modeling: Equation sets
    and test cases", J. Comput. Phys., 227, 2008, 3849--3877, 
    (Section 3.1).

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See the following reference:
+ Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
  Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
  Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
  http://dx.doi.org/10.1137/15M1044369.

Domain: \f$0 \le x \le 300,000\,m, 0 \le y \le 10,000\,m\f$, \a "periodic" (#_PERIODIC_)
        boundary conditions along \f$x\f$, "slip-wall" (#_SLIP_WALL_) boundary conditions along \f$y\f$.

Initial solution: See references above.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes2D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes2D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes2D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes2D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes2D::p0)
  + Hydrostatic balance type 3 (#NavierStokes2D::HB)
  + Brunt-Vaisala frequency 0.01 (#NavierStokes2D::N_bv)

Numerical method:
 + Spatial discretization (hyperbolic): 5th order compact upwind (Interp1PrimFifthOrderCompactUpwind())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>Kennedy-Carpenter ARK4</B> (TSARKIMEX4 - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX4.html)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 32 (DMDROMObject::m_rdim)
 + Sampling frequency: 1 (libROMInterface::m_sampling_freq)
 + Number of samples per time window: 50 (DMDROMObject::m_num_window_samples)


Input files required:
---------------------

<B>librom.inp</B>
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/librom.inp

<B>.petscrc</B>
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/petscrc

\b solver.inp
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/physics.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/aux/init.c

Output:
-------
Note that \b iproc is set to 

      24 1

in \b solver.inp (i.e., 24 processors along \a x, and 1
processor along \a y). Thus, this example should be run
with 24 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 16 output file \b op_00000.bin, ..., \b op_00015.bin; 
this is the \b HyPar solutions\b. 

+ 16 output file \b op_rom_00000.bin, ..., \b op_rom_00015.bin; 
this is the \b predicted solutions from the DMD object(s)\b.

All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

#HyPar::op_file_format is set to \a binary in \b solver.inp, and
thus, all the files are written out in the binary format, see 
WriteBinary(). The binary file contains the conserved variables
\f$\left(\rho, \rho u, \rho v, e\right)\f$. The following code
converts these variables to the primitive variables of interest
to atmospheric flows \f$\left(\rho, u, v, p, \theta\right)\f$.
It also writes out the hydrostatically balanced quantities 
\f$\left(\rho_0,\pi_0, \theta_0\right)\f$ for this case that
can be used to compute and plot the temperature and density
perturbations. These variables are then written to either
a tecplot2d or text file.
(compile and run it in the run directory):
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/aux/PostProcess.c

Alternatively, the provided Python script (\b plotSolution.py)
can be used to generate plots from the binary files of the atmospheric
flow variables and compare the HyPar and DMD
solutions. It will plot the potential temperature perturbation but can 
be modified easily to plot other quantities.

The following plot shows the potential temperature perturbation
contours at the final time t=3000 (FOM is the HyPar solution, ROM is 
the DMD prediction, and diff is the difference between the two):
@image html Solution_2DNavStokIGWavePETSc_libROMDMD_Train.png

\b Wall \b clock \b times:
- PDE solution: 521 seconds
- DMD training time: 15.1 seconds
- DMD prediction/query time: 8.7 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 2D/NavierStokes2D/InertiaGravityWave_PETSc_IMEX_libROM_DMD_Train/out.log

\page ns3d_bubble_petsc_librom_dmd_train 3D Navier-Stokes Equations (with gravitational force) - Rising Thermal Bubble (Time-Windowed DMD) with PETSc IMEX time integration

See \ref ns3d_bubble_petsc to familiarize yourself with this case.

Location: \b hypar/Examples/3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: \f$0 \le x,y,z < 1000\,{\rm m}\f$, \a "slip-wall" (#_SLIP_WALL_) boundaries 
        everywhere, with zero wall velocity.

Reference:
  + Kelly, J. F., Giraldo, F. X., "Continuous and discontinuous Galerkin methods for a scalable
  three-dimensional nonhydrostatic atmospheric model: Limited-area mode", J. Comput. Phys., 231,
  2012, pp. 7988-8008 (see section 5.1.2).
  + Giraldo, F. X., Kelly, J. F., Constantinescu, E. M., "Implicit-Explicit Formulations of a
  Three-Dimensional Nonhydrostatic Unified Model of the Atmosphere (NUMA)", SIAM J. Sci. Comput., 
  35 (5), 2013, pp. B1162-B1194 (see section 4.1).

The problem is solved here using <B>implicit-explicit (IMEX)</B> time
integration, where the hyperbolic flux is partitioned into its entropy
and acoustic components with the former integrated explicitly and the
latter integrated implicitly. See the following reference:
+ Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of 
  Atmospheric Flows with Characteristic-Based Flux Partitioning", SIAM 
  Journal on Scientific Computing, 38 (3), 2016, A1848-A1875, 
  http://dx.doi.org/10.1137/15M1044369.

Initial solution: A warm bubble in cool ambient atmosphere. Note that in this example, the
gravitational forces and rising of the bubble is along the \a y-axis.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes3D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes3D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes3D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes3D::p0)
  + Hydrostatic balance type 2 (#NavierStokes3D::HB)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: PETSc (SolvePETSc()) 
   - Method Class: <B>Additive Runge-Kutta method</B> (TSARKIMEX - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html)
   - Specific method: <B>ARK 2e</B> (TSARKIMEX2E - https://petsc.org/release/docs/manualpages/TS/TSARKIMEX2E.html)

Reduced Order Modeling:
 + Mode: train (libROMInterface::m_mode, #_ROM_MODE_TRAIN_)
 + Component Mode: monolithic (libROMInterface::m_comp_mode, #_ROM_COMP_MODE_MONOLITHIC_)
 + Type: Dynamic Mode Decomposition (DMD) with time windowing (libROMInterface::m_rom_type)
 + Latent subspace dimension: 16 (DMDROMObject::m_rdim)
 + Sampling frequency: 1 (libROMInterface::m_sampling_freq)
 + Number of samples per time window: 50 (DMDROMObject::m_num_window_samples)

Input files required:
---------------------

<B>.petscrc</B>
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/petscrc

\b librom.inp
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/librom.inp

\b solver.inp
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/physics.inp

\b weno.inp (optional)
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4 4

in \b solver.inp (i.e., 4 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 64 MPI ranks (or change \b iproc).

After running the code, there should be the following output
files:

+ 11 output file \b op_00000.bin, ..., \b op_00010.bin; 
this is the \b HyPar solutions\b. 

+ 11 output file \b op_rom_00000.bin, ..., \b op_rom_00010.bin; 
this is the \b predicted solutions from the DMD object(s)\b.

All the files are binary
(#HyPar::op_file_format is set to \a binary in \b solver.inp).

The binary solution file contains the conserved variables
(\f$\rho,\rho u,\rho v,\rho w,e\f$). 

The Python script plotSolution.py calculates the primitive
and reference variables 
\f$\rho,u,v,w,P,\theta,\rho_0,P_0,\pi,\theta_0\f$ (where the
subscript \f$0\f$ indicates the hydrostatic mean value) for both
the op_* and op_rom_* files. It plots 3D figures for \f$\theta\f$ 
(potential temperature) for the PDE solution, the DMD predictions,
and the difference between the two.

The following figure shows the potential temperature at the final time. Colored
contours are plotted in the X-Y plane, and contour lines are plotted in the Z-Y
plane (FOM is the HyPar solution, ROM is the DMD prediction, and diff is the 
difference between the two):
@image html Solution_3DNavStok_Bubble3D_PETSc_libROMDMD_Train.png

\b Wall \b clock \b times:
- PDE solution: 906 seconds
- DMD training time: 22.5 seconds
- DMD prediction/query time: 4.5 seconds

The L1, L2, and Linf norms of the diff between the HyPar and ROM solution 
at the final time are calculated and reported on screen (see below)
as well as \b pde_rom_diff.dat:
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/pde_rom_diff.dat
The numbers are: number of grid points in each dimension (#HyPar::dim_global), 
number of processors in each dimension (#MPIVariables::iproc),
time step size (#HyPar::dt),
L1, L2, and L-infinity norms of the diff (#HyPar::rom_diff_norms),
solver wall time (seconds) (i.e., not accounting for initialization,
and cleaning up),
and total wall time.

Expected screen output:
\include 3D/NavierStokes3D/RisingThermalBubble_PETSc_IMEX_libROM_DMD_Train/out.log

\page gpu_examples CUDA-enabled Examples

Implementation of HyPar on GPU is a work-in-progress. Currently, only the following 
parts of the code has been implemented on NVidia GPUs with CUDA:

+ Spatial discretization schemes:
  + 5th order WENO scheme (Interp1PrimFifthOrderWENO())
  
+ Time integration schemes:
  + Runge-Kutta time integration (TimeRK())

+ Boundary conditions:
  + Periodic (#_PERIODIC_)
  + "Slip-wall" (#_SLIP_WALL_)

+ Physical models:
  + 3D Navier-Stokes (#NavierStokes3D)

Implementation note: 
+ The CUDA implementation uses a different variable ordering than the
  CPU code. While the CPU code uses the "Array of Structures" variable ordering, where
  the innermost loop is over the vector components at a grid point, the CUDA implementation
  uses the "Structure of Arrays" ordering, where the innermost loops are over the grid
  points.
+ This implementation is designed for 1 GPU per MPI rank. Thus, on a node, the
  simulation should be run with as many MPI ranks as there are GPUs.

\subpage ns3d_isoturb_gpu \n
\subpage ns3d_bubble_gpu

\page ns3d_isoturb_gpu 3D Navier-Stokes Equations - Isotropic Turbulence Decay

See \ref ns3d_isoturb to familiarize yourself with this case.

Location: \b hypar/Examples/3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay_CUDA
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: \f$0 \le x,y,z < 2\pi\f$, "periodic" (#_PERIODIC_) boundaries 
        everywhere.

Initial solution: Isotropic turbulent flow - The initial solution is 
specified in the Fourier space (with an energy
distribution similar to that of turbulent flow), and then transformed
to the physical space through an inverse transform.

Other parameters:
  + \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + \f$Re = \frac {\rho u L } {\mu} = 333.33\f$ (#NavierStokes3D::Re)
  + \f$Pr = 0.72\f$ (Prandtl number) (#NavierStokes3D::Pr)
  + \f$M_\infty = 0.3\f$ (turbulence fluctuation Mach number) (#NavierStokes3D::Minf)

Numerical Method:
  + Spatial discretization (hyperbolic): 5th order CRWENO (Interp1PrimFifthOrderCRWENO())
  + Spatial discretization (parabolic) : 4th order (FirstDerivativeFourthOrderCentral()) 
                                         non-conservative 2-stage (ParabolicFunctionNC2Stage())
  + Time integration: RK4 (TimeRK(), #_RK_44_)

Hardware Details:
  + GPU: NVIDIA V100 (Volta)
  + CPU: IBM Power9 
  + Configuration: 16 nodes, 4 GPUs and 4 MPI ranks per node (64 GPUs and 64 MPI ranks).

Input files required:
---------------------

\b solver.inp
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay_CUDA/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay_CUDA/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay_CUDA/physics.inp

\b weno.inp (optional)
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay_CUDA/weno.inp

To generate \b initial.inp (initial solution), compile 
and run the following code in the run directory.
\b Note: this code requires the \b FFTW library installed (http://www.fftw.org/).
To compile:

    gcc -I/path/to/fftw3.h -L/path/to/libfftw3.a -lfftw3 init.c

(see the FFTW website on ways to install it).
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay_CUDA/aux/init.c

Output:
-------
Note that \b iproc is set to 

      4 4 4

in \b solver.inp (i.e., 4 processors along \a x, 4
processors along \a y, and 4 processor along \a z). Thus, 
this example should be run with 64 MPI ranks (or change \b iproc)
and 64 GPUs.

After running the code, there should be 11 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00010.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=5\f$.

See \ref ns3d_isoturb for detailed postprocessing instructions. The 
following figure shows the initial and final (t=5) energy spectra:
@image html Solution_3DNavStok_IsoTurb_CUDA_Spectrum.png

The following figures show the density and internal energy contours
on the X-Y plane mid-slice at the final time:
@image html Solution_3DNavStok_IsoTurb_CUDA_density.png
@image html Solution_3DNavStok_IsoTurb_CUDA_energy.png

Expected screen output:
\include 3D/NavierStokes3D/DNS_IsotropicTurbulenceDecay_CUDA/out.log

\page ns3d_bubble_gpu 3D Navier-Stokes Equations - Rising Thermal Bubble

See \ref ns3d_bubble to familiarize yourself with this case.

Location: \b hypar/Examples/3D/NavierStokes3D/RisingThermalBubble_Config1_CUDA
          (This directory contains all the input files needed
          to run this case.)

Governing equations: 3D Navier-Stokes Equations (navierstokes3d.h)

Domain: \f$0 \le x,y,z < 1000\,{\rm m}\f$, \a "slip-wall" (#_SLIP_WALL_) boundaries 
        everywhere, with zero wall velocity.

Reference:
  + Kelly, J. F., Giraldo, F. X., "Continuous and discontinuous Galerkin methods for a scalable
  three-dimensional nonhydrostatic atmospheric model: Limited-area mode", J. Comput. Phys., 231,
  2012, pp. 7988-8008 (see section 5.1.2).
  + Giraldo, F. X., Kelly, J. F., Constantinescu, E. M., "Implicit-Explicit Formulations of a
  Three-Dimensional Nonhydrostatic Unified Model of the Atmosphere (NUMA)", SIAM J. Sci. Comput., 
  35 (5), 2013, pp. B1162-B1194 (see section 4.1).

Initial solution: A warm bubble in cool ambient atmosphere. Note that in this example, the
gravitational forces and rising of the bubble is along the \a y-axis.

Other parameters (all dimensional quantities are in SI units):
  + Specific heat ratio \f$\gamma = 1.4\f$ (#NavierStokes3D::gamma)
  + Universal gas constant \f$R = 287.058\f$ (#NavierStokes3D::R)
  + Gravitational force per unit mass \f$g = 9.8\f$ along \a y-axis (#NavierStokes3D::grav_y)
  + Reference density (at zero altitude) \f$\rho_{ref} = 1.1612055171196529\f$ (#NavierStokes3D::rho0)
  + Reference pressure (at zero altitude) \f$P_{ref} = 100000\f$ (#NavierStokes3D::p0)
  + Hydrostatic balance type 2 (#NavierStokes3D::HB)

Numerical Method:
 + Spatial discretization (hyperbolic): 5th order WENO (Interp1PrimFifthOrderWENO())
 + Time integration: SSP RK3 (TimeRK(), #_RK_SSP3_)

Hardware Details:
  + GPU: NVIDIA V100 (Volta)
  + CPU: IBM Power9 
  + Configuration: 1 node, 1 GPU and 1 MPI rank

Input files required:
---------------------

\b solver.inp
\include 3D/NavierStokes3D/RisingThermalBubble_Config1_CUDA/solver.inp

\b boundary.inp
\include 3D/NavierStokes3D/RisingThermalBubble_Config1_CUDA/boundary.inp

\b physics.inp
\include 3D/NavierStokes3D/RisingThermalBubble_Config1_CUDA/physics.inp

\b weno.inp (optional)
\include 3D/NavierStokes3D/RisingThermalBubble_Config1_CUDA/weno.inp

To generate \b initial_par.inp.* (initial solution), compile 
and run the following code in the run directory.
\include 3D/NavierStokes3D/RisingThermalBubble_Config1_CUDA/aux/init_parallel.c

Output:
-------
Note that \b iproc is set to 

      1 1 1

in \b solver.inp This example should be run with 1 MPI rank
and 1 GPU.

After running the code, there should be 17 output
files \b op_00000.bin, \b op_00001.bin, ... \b op_00016.bin; 
the first one is the solution at \f$t=0\f$ and the final one
is the solution at \f$t=200\,{\rm s}\f$. 

See \ref ns3d_bubble for detailed postprocessing instructions.  
The following figure is generated using aux/plotSolution.py and shows the 
potential temperature at the final time. Colored
contours are plotted in the X-Y plane, and contour lines are plotted in the Z-Y
plane:
@image html Solution_3DNavStok_Bubble3D_CUDA.png

Expected screen output:
\include 3D/NavierStokes3D/RisingThermalBubble_Config1_CUDA/out.log

