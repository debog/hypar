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
                  const int a_var_idx = -1);

    /*! Destructor */
    virtual ~LSROMObject()
    {
//    for (int i = 0; i < m_ls.size(); i++) delete m_ls[i];
//    m_ls.clear();
      m_intervals.clear();
      m_ls_is_trained.clear();
    }

    /*! Project initial solution for prediction */
    virtual void projectInitialSolution(  CAROM::Vector& );

    /*! take a sample (solution snapshot) */
    virtual void takeSample( const CAROM::Vector&, const double, void* );

    /*! train the LS object */
    virtual void train(void* );

    /*! compute prediction at given time */
    virtual 
    const CAROM::Vector* predict(const double );

    /*! save LS object to file */
    virtual void save(const std::string& a_fname_root /*!< Filename root*/) const;

    /*! load LS object from file */
    virtual void load(const std::string& a_fname_root /*!< Filename root*/);

    /*! Copy a vector to HyPar */
    void copyToHyPar( const CAROM::Vector&, void*, int ) const;

    /*! Write out libROM field */
    void OutputlibROMfield(double*, SimulationObject&, char*);

    /*! Initialize the time integration for LSROM */
    int TimeInitialize ();
    /* RK Time integration for LSROM */
    int TimeExplicitRKInitialize();
    int TimeRK(const double );

  protected:

    std::vector<CAROM::Options*> m_options; /*!< Vector of Options objects */
    std::vector<CAROM::BasisGenerator*> m_generator; /*!< Vector of BasisGenerator objects */
    std::vector<CAROM::Vector*> m_projected_init; /*!< Vector of Vector objects */
    const CAROM::Vector* m_S; /*!< Vector of Singular value */
    CAROM::Vector* m_romcoef; /*!< Vector of rom coefficients */
    CAROM::Matrix* m_snapshots; /*!< Snapshot Matrix */
    CAROM::Matrix* m_romhyperb; /*!< Reduced Hyperbolic Term Matrix */

    std::vector<bool> m_ls_is_trained; /*!< Flag to indicate if LS is trained */
    std::vector<Interval> m_intervals; /*!< Time intervals for each LS object */

    int m_rank;   /*!< MPI rank */
    int m_nproc;  /*!< Number of MPI ranks */

    int     m_vec_size; /*!< Local size of solution vector */
    double  m_dt;       /*!< Time step size */
    double  m_t_final;  /*!< Final time */
    int     m_rdim;     /*!< Latent space dimension */

    int m_num_window_samples; /*!< Number of samples per LS for time-windowing */

    int m_tic; /*!< private ticker to count number of samples taken */
    int m_curr_win; /*!< private index for current window */

    int m_sim_idx; /*!< simulation index of this object for ensemble simulations */
    int m_var_idx; /*!< component index of this object if component-wise ROMs are being used */

    bool m_write_snapshot_mat;  /*!< Write snapshot matrix to file or not */

    std::string m_dirname; /*!< Subdirectory where LS objects are written to or read from */

    /* parameters for generator and option */
    int max_num_snapshots = 100000;
    bool update_right_SV = false;
    bool isIncremental = false;
    const std::string basisName = "basis";
    const std::string basisFileName;
    int numRowRB, numColumnRB;

    /* parameters for RK, similar as in timeintegration_struct.h */
    int nstages; /*!< number of stages */
    double *A, /*!< Stage computation coefficients (Butcher tableau form),
                   saved as a 1D-array in row-major form */
           *b, /*!< Step completion coefficients (Butcher tableau form) */
           *c; /*!< Stage time coefficients (Butcher tableau form) */
    /*! Arrays to store stage values for a multi-stage time-integration method */
    std::vector<CAROM::Vector*> m_U;

    /*! Arrays to store stage right-hand-sides for a multi-stage time-integration method */
    std::vector<CAROM::Vector*> m_Udot;

  private:
};

#endif

#endif
