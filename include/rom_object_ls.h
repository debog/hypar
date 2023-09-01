#ifdef with_librom

/*! @file rom_object_ls.h
    @brief Linear Subpsace ROM object
    @author Ping-Hsuan Tsai
*/

#ifndef _ROM_OBJECT_LS_H_
#define _ROM_OBJECT_LS_H_

/*! LS-type ROM */
#define _ROM_TYPE_LS_ "LS"

#include <climits>
#include <vector>
#include <utility>
#include <linalg/Matrix.h>
#include "linalg/Vector.h"
#include "linalg/BasisGenerator.h"
#include "linalg/Options.h"
#include <rom_object.h>
#include <simulation_object.h>
#include <fstream>

#ifndef serial
#include <mpi.h>
#endif

/*! Interval type */
typedef std::pair<double,double> Interval;

class LSROMObject : public ROMObject
{
  public:
    /*! Constructor */
    LSROMObject(  const int, 
                  const double, 
                  const int, 
                  const int,
                  const int,
                  const int, 
                  const int,
                  const int,
                  const int a_var_idx = -1);

    /*! Destructor */
    virtual ~LSROMObject()
    {
//    for (int i = 0; i < m_projected_init.size(); i++) delete m_projected_init[i];
//    m_projected_init.clear();
      for (int i = 0; i < m_snapshots.size(); i++) delete m_snapshots[i];
      m_snapshots.clear();
      for (int i = 0; i < m_basis.size(); i++) delete m_basis[i];
      m_basis.clear();
      for (int i = 0; i < m_romhyperb.size(); i++) delete m_romhyperb[i];
      m_romhyperb.clear();
      for (int i = 0; i < m_romcoef.size(); i++) delete m_romcoef[i];
      m_romcoef.clear();

      for (int i = 0; i < m_projected_init_phi.size(); i++) delete m_projected_init_phi[i];
      m_projected_init_phi.clear();
      for (int i = 0; i < m_snapshots_phi.size(); i++) delete m_snapshots_phi[i];
      m_snapshots_phi.clear();
      for (int i = 0; i < m_basis_phi.size(); i++) delete m_basis_phi[i];
      m_basis_phi.clear();
      for (int i = 0; i < m_romrhs_phi.size(); i++) delete m_romrhs_phi[i];
      m_romrhs_phi.clear();
      for (int i = 0; i < m_romlaplace_phi.size(); i++) delete m_romlaplace_phi[i];
      m_romlaplace_phi.clear();

      for (int i = 0; i < m_options_e.size(); i++) delete m_options_e[i];
      m_options_e.clear();
      for (int i = 0; i < m_generator_e.size(); i++) delete m_generator_e[i];
      m_generator_e.clear();
      for (int i = 0; i < m_projected_init_e.size(); i++) delete m_projected_init_e[i];
      m_projected_init_e.clear();
      for (int i = 0; i < m_snapshots_e.size(); i++) delete m_snapshots_e[i];
      m_snapshots_e.clear();
      for (int i = 0; i < m_basis_e.size(); i++) delete m_basis_e[i];
      m_basis_e.clear();

      for (int i = 0; i < m_romhyperb_x.size(); i++) delete m_romhyperb_x[i];
      m_romhyperb_x.clear();
      for (int i = 0; i < m_romhyperb_v.size(); i++) {
        for (int j = 0; j < m_romhyperb_v[i].size(); j++) {
          delete m_romhyperb_v[i][j];
        }
        m_romhyperb_v[i].clear();
      }
      m_romhyperb_v.clear();

      /* Delete allocation used in time integration */
      free(A);
      free(b);
      free(c);
      for (int i = 0; i < m_U.size(); i++) delete m_U[i];
      m_U.clear();
      for (int i = 0; i < m_Udot.size(); i++) delete m_Udot[i];
      m_Udot.clear();

      m_intervals.clear();
      m_ls_is_trained.clear();
			printf("Calling LSROMObject desctructor\n");
    }

    /*! Project initial solution for prediction */
    virtual void projectInitialSolution(  CAROM::Vector&, void* );

    /*! take a sample (solution snapshot) */
    virtual void takeSample( const CAROM::Vector&, const double, void* );

    /*! train the LS object */
    virtual void train(void* );

    /*! compute prediction at given time */
    virtual 
    const CAROM::Vector* predict(const double, void* );

    /*! save LS object to file */
    virtual void save(const std::string& a_fname_root /*!< Filename root*/) const;

    /*! load LS object from file */
    virtual void load(const std::string& a_fname_root /*!< Filename root*/);

    /*! Copy a vector to HyPar */
    void copyToHyPar( const CAROM::Vector&, void*, int ) const;

    /*! Write out libROM field */
    void OutputlibROMfield(double*, SimulationObject&, char*);

    /*! Check the LS object */
    virtual void check(void* );

    /*! Reconstruct libROM field from coefficient */
    CAROM::Vector* ReconlibROMfield(const CAROM::Vector*, const CAROM::Matrix*, const int );

    /*! Project a libROM field onto reduced space*/
    CAROM::Vector* ProjectToRB(CAROM::Vector*, const CAROM::Matrix*, const int );

    /*! Initialize the time integration for LSROM */
    int TimeInitialize (int);
    /* RK Time integration for LSROM */
    int TimeExplicitRKInitialize();

    int TimeRK(const double, void*, int);

    /* Construct reduced hyperbolic operator */
    void ConstructROMHy(void*, const CAROM::Matrix*,int);

    /* Dump ROM basis */
    void OutputROMBasis(void*, const CAROM::Matrix*,int);

    int CalSnapROMDiff(void*, void*, double*, double*, char*);

    void CheckSolProjError(void*,int);

    void CheckHyProjError(void*,int);

    /* Construct reduced rhs operator for potential */
    void ConstructPotentialROMRhs(void*, const CAROM::Matrix*, const CAROM::Matrix*, int);

    /* Construct reduced laplace operator for potential */
    void ConstructPotentialROMLaplace(void*, const CAROM::Matrix*, int);

    void CheckPotentialProjError(void*,int);

    void CheckLaplaceProjError(void*,int);

    void CheckRhsProjError(void*,int);

    int CalSnapROMDiff_phi(void*, void*, double*, double*, char*);

    /*! Project initial solution for phi */
    void projectInitialSolution_phi(  CAROM::Vector& );

    /*! Evaluate Rhs */
    void EvaluatePotentialRhs(void*, CAROM::Vector*, double*);

    void CheckEProjError(void*,int);

    /* Dump ROM basis phi */
    void OutputROMBasisPhi(void*, const CAROM::Matrix*,int);

    void ConstructEBasis(void*,int);

    /* Construct reduced hyperbolic operator in x direction */
    void ConstructROMHy_x(void*, const CAROM::Matrix*,int);

    /* Construct reduced hyperbolic tensor in v direction */
    void ConstructROMHy_v(void*, const CAROM::Matrix*, const CAROM::Matrix*,int);

    /*! Write Snapshot matrix */
    virtual void writeSnapshot(void* );

    /*! Merge stage */
    virtual void merge(void* );

    /*! Online stage */
    virtual void online(void* );

    /* Read Time Window Paramters */
    void ReadTimeWindows(void*);

    /* Read Time interval for each time window */
    void ReadTimeIntervals(void*);

    /* Clean up allocation related to LSROM */
    virtual void cleanup(void*);

    /* Clean up allocation related to time integration */
    int TimeCleanup();

  protected:

    std::vector<CAROM::Options*> m_options; /*!< Vector of Options objects */
    std::vector<CAROM::BasisGenerator*> m_generator; /*!< Vector of BasisGenerator objects */
    std::vector<CAROM::Vector*> m_projected_init; /*!< Vector of Vector objects */
    std::vector<const CAROM::Matrix*> m_snapshots; /*!< Snapshot Matrix */
    std::vector<const CAROM::Matrix*> m_basis; /*!< Snapshot Matrix */
    std::vector<CAROM::Matrix*> m_romhyperb; /*!< Reduced Hyperbolic Term Matrix */
    const CAROM::Vector* m_S; /*!< Vector of Singular value */
    std::vector<CAROM::Vector*> m_romcoef; /*!< Vector of rom coefficients */
    CAROM::Vector* m_romhycoef; /*!< Vector of rom coefficients */
    CAROM::Vector* m_recon; /*!< Vector of rom coefficients */

    std::vector<bool> m_ls_is_trained; /*!< Flag to indicate if LS is trained */
    std::vector<Interval> m_intervals; /*!< Time intervals for each LS object */

    int m_rank;   /*!< MPI rank */
    int m_nproc;  /*!< Number of MPI ranks */

    int     m_vec_size; /*!< Local size of solution vector */
    double  m_dt;       /*!< Time step size */
    double  m_t_final;  /*!< Final time */
    int     m_rdim;     /*!< Latent space dimension */
    int     m_rdim_phi;     /*!< Phi Latent space dimension */
    std::vector<int> m_rdims;     /*!< Latent space dimension */
    std::vector<int> m_rdims_phi;     /*!< Phi Latent space dimension */

    int m_num_window_samples; /*!< Number of samples per LS for time-windowing */

    int m_tic; /*!< private ticker to count number of samples taken */
    int m_curr_win; /*!< private index for current window */

    int m_sim_idx; /*!< simulation index of this object for ensemble simulations */
    int m_var_idx; /*!< component index of this object if component-wise ROMs are being used */

    bool m_write_snapshot_mat;  /*!< Write snapshot matrix to file or not */
    bool m_direct_comp_hyperbolic;  /*!< Compute Hyperbolic term directly or not */
    bool m_solve_phi;  /*!< Solve potential or not */
    bool m_c_err_snap;  /*!< Compute error of each snapshot*/

    double m_f_energy_criteria;
    double m_phi_energy_criteria;

    std::string m_dirname; /*!< Subdirectory where LS objects are written to or read from */

    /* parameters for generator and option */
    int max_num_snapshots = 100000;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "basis";
    const std::string basisFileName;
    const std::string basisName_phi = "basis_phi";
    const std::string basisFileName_phi;
    int numRowRB, numColumnRB;

    /* parameters for RK, similar as in timeintegration_struct.h */
    int nstages; /*!< number of stages */
    double *A=0; /*!< Stage computation coefficients (Butcher tableau form),
                   saved as a 1D-array in row-major form */

    double *b=0; /*!< Step completion coefficients (Butcher tableau form) */
    double *c=0; /*!< Stage time coefficients (Butcher tableau form) */
    /*! Arrays to store stage values for a multi-stage time-integration method */
    std::vector<CAROM::Vector*> m_U;

    /*! Arrays to store stage right-hand-sides for a multi-stage time-integration method */
    std::vector<CAROM::Vector*> m_Udot;

    std::vector<int> m_snap; /*!< private ticker to count number of samples taken */

    /* Below for precomputation idea */
    std::vector<CAROM::Options*> m_options_phi; /*!< Vector of Options objects */
    std::vector<CAROM::BasisGenerator*> m_generator_phi; /*!< Vector of BasisGenerator objects */
    std::vector<CAROM::Vector*> m_projected_init_phi; /*!< Vector of Vector objects */
    std::vector<CAROM::Matrix*> m_snapshots_phi; /*!< Snapshot Matrix */
    std::vector<const CAROM::Matrix*> m_basis_phi; /*!< Snapshot Matrix */
    std::vector<CAROM::Matrix*> m_romrhs_phi; /*!< Reduced Hyperbolic Term Matrix */
    std::vector<CAROM::Matrix*> m_romlaplace_phi; /*!< Reduced Hyperbolic Term Matrix */
    const CAROM::Vector* m_S_phi; /*!< Vector of Singular value */

    /* Below for verify projection error in E*/
    std::vector<CAROM::Options*> m_options_e; /*!< Vector of Options objects */
    std::vector<CAROM::BasisGenerator*> m_generator_e; /*!< Vector of BasisGenerator objects */
    std::vector<CAROM::Vector*> m_projected_init_e; /*!< Vector of Vector objects */
    std::vector<CAROM::Matrix*> m_snapshots_e; /*!< Snapshot Matrix */
    std::vector<CAROM::Matrix*> m_basis_e; /*!< Basis Matrix */

    std::vector<CAROM::Matrix*> m_romhyperb_x; /*!< Reduced Hyperbolic Term Matrix */
    std::vector<std::vector<CAROM::Matrix*>> m_romhyperb_v; /*!< Reduced Hyperbolic Term Matrix */
    int m_sampling_freq;       /*!< Time step size */

    std::vector<CAROM::Matrix*> m_fullscale;
    CAROM::Matrix* m_matrix;

    int m_parametric_id; /*!< frequency (iterations) of writing ROM solution to file */

    /* rightnow m_nsets only for LSROM */
    int m_nsets; /*!< Number of parametric snapshot sets */

    std::ofstream outfile_twp;
    std::string outputPath = ".";
    const char *twpfile = "twp.csv";
    const char *twintervalfile = "twinterval.csv";
    int m_numwindows, m_numwindows_phi;

  private:
};

#endif

#endif
