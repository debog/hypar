#ifdef with_librom

/*! @file rom_object_dmd.h
    @brief Dynamic Mode Decomposition ROM object
    @author Debojyoti Ghosh
*/

#ifndef _ROM_OBJECT_DMD_H_
#define _ROM_OBJECT_DMD_H_

/*! DMD-type ROM */
#define _ROM_TYPE_DMD_ "DMD"

#include <climits>
#include <vector>
#include <utility>
#include <linalg/Matrix.h>
#include <algo/DMD.h>
#include <rom_object.h>

#ifndef serial
#include <mpi.h>
#endif

/*! Interval type */
typedef std::pair<double,double> Interval;

/*! \class DMDROMObject
 *  \brief ROM object of type DMD (see libROM)
 *  DMD-type ROM (see libROM for details)
*/
/*! \brief ROM object of type DMD (see libROM)
 *
 *  DMD-type ROM (see libROM for details). DMD stands
 *  for "Dynamic Mode Decomposition".
*/
class DMDROMObject : public ROMObject
{
  public:

    /*! Constructor */
    DMDROMObject( const int,
                  const double,
                  const int,
                  const int,
                  const int,
                  const int,
                  const int a_var_idx = -1);

    /*! Destructor */
    virtual ~DMDROMObject()
    {
      for (int i = 0; i < m_dmd.size(); i++) delete m_dmd[i];
      m_dmd.clear();
      m_intervals.clear();
      m_dmd_is_trained.clear();
    }

    /*! Project initial solution for prediction */
    virtual void projectInitialSolution(  CAROM::Vector& a_U /*!< solution vector */ )
    {
      if (m_dmd.size() == 0) {
        if (!m_rank) {
          printf("ERROR in DMDROMObject::projectInitialSolution() - m_dmd is a vector of size 0.\n");
        }
        return;
      }

      m_dmd[0]->projectInitialCondition( &a_U );
      for (int i = 1; i < m_dmd.size(); i++) {
        m_dmd[i]->projectInitialCondition( m_dmd[i-1]->predict(m_intervals[i].first) );
      }
      return;
    }

    /*! take a sample (solution snapshot) */
    virtual void takeSample( const CAROM::Vector&, const double );

    /*! train the DMD object */
    virtual void train();

    /*! compute prediction at given time */
    virtual
    const CAROM::Vector* const predict(const double a_t /*!< time at which to predict solution */ ) const
    {
      for (int i = 0; i < m_dmd.size(); i++) {
        if (   (a_t >= m_intervals[i].first)
            && (  (a_t < m_intervals[i].second) || (m_intervals[i].second < 0)  ) ){
          return m_dmd[i]->predict(a_t);
        }
      }
      printf("ERROR in DMDROMObject::predict(): m_dmd is of size zero or interval not found!");
      return nullptr;
    }

    /*! save DMD object to file */
    virtual void save(const std::string& a_fname_root /*!< Filename root*/) const;

    /*! load DMD object from file */
    virtual void load(const std::string& a_fname_root /*!< Filename root*/);

  protected:

    std::vector<CAROM::DMD*> m_dmd; /*!< Vector of DMD objects */
    std::vector<bool> m_dmd_is_trained; /*!< Flag to indicate if DMD is trained */
    std::vector<Interval> m_intervals; /*!< Time intervals for each DMD object */

    int m_rank;   /*!< MPI rank */
    int m_nproc;  /*!< Number of MPI ranks */

    int     m_vec_size; /*!< Local size of solution vector */
    double  m_dt;       /*!< Time step size */
    double  m_t_final;  /*!< Final time */
    int     m_rdim;     /*!< Latent space dimension */

    int m_num_window_samples; /*!< Number of samples per DMD for time-windowing */

    int m_tic; /*!< private ticker to count number of samples taken */
    int m_curr_win; /*!< private index for current window */

    int m_sim_idx; /*!< simulation index of this object for ensemble simulations */
    int m_var_idx; /*!< component index of this object if component-wise ROMs are being used */

    bool m_write_snapshot_mat;  /*!< Write snapshot matrix to file or not */

    std::string m_dirname; /*!< Subdirectory where DMD objects are written to or read from */

  private:
};

#endif

#endif
