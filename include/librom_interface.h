#ifdef with_librom

/*! @file librom_interface.h
    @brief libROM interface class
    @author Debojyoti Ghosh
*/

#ifndef _LIBROM_INTERFACE_H_
#define _LIBROM_INTERFACE_H_

#include <climits>
#include <cfloat>
#include <string>
#include <vector>
#include <utility>
#include <sys/time.h>

#include <linalg/Vector.h>
#include <linalg/Matrix.h>
#include <algo/DMD.h>

/*! Filename for libROM-related inputs */
#define _LIBROM_INP_FNAME_ "librom.inp"

/*! DMD-type ROM */
#define _ROM_TYPE_DMD_ "DMD"

/*! Interval type */
typedef std::pair<double,double> Interval;

/*! \class ROMObject
 *  \brief Base class defining a ROM object.
 *  This is the base class that defines a 
 *  ROM object.
*/
/*! \brief Base class defining a ROM object.
 *
 *  This is the base class that defines a 
 *  ROM object.
*/
class ROMObject
{
  public:

    /*! take a sample (solution snapshot) */
    virtual void takeSample( double*, const double ) = 0;
    /*! train the ROM object */
    virtual void train() = 0;
    /*! compute prediction at given time */
    virtual const CAROM::Vector* const predict( const double ) = 0;
    /*! save ROM object to file */
    virtual void save( const std::string& ) = 0;

  protected:

  private:

};

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
    DMDROMObject( const int, const double, const int, const int, const int);

    /*! Destructor */
    virtual ~DMDROMObject()
    {
      m_dmd.clear();
      m_intervals.clear();
    }

    /*! take a sample (solution snapshot) */
    virtual void takeSample(  double* a_U, /*!< solution vector */
                              const double a_time /*!< sample time */ )
    {
      if (m_tic == 0) {

        m_dmd.push_back( new CAROM::DMD(m_vec_size, m_dt) );
        m_intervals.push_back( Interval(a_time, DBL_MAX) );
        if (!m_rank) {
          printf("DMDROMObject::takeSample() - creating new DMD object, t=%f (total: %d).\n",
                 m_intervals[m_curr_win].first, m_dmd.size());
        }
        m_dmd[m_curr_win]->takeSample( a_U, a_time );

      } else {

        m_dmd[m_curr_win]->takeSample( a_U, a_time );
        if (m_tic%m_num_window_samples == 0) {
          m_intervals[m_curr_win].second = a_time;
          m_curr_win++;
          m_dmd.push_back( new CAROM::DMD(m_vec_size, m_dt) );
          m_intervals.push_back( Interval(a_time, DBL_MAX) );
          m_dmd[m_curr_win]->takeSample( a_U, a_time );
          if (!m_rank) {
            printf("DMDROMObject::takeSample() - creating new DMD object, t=%f (total: %d).\n",
                   m_intervals[m_curr_win].first, m_dmd.size());
          }
        }

      }

      m_tic++;
      return;
    }

    /*! train the DMD object */
    virtual void train();

    /*! compute prediction at given time */
    virtual 
    const CAROM::Vector* const predict(const double a_t /*!< time at which to predict solution */ )
    {
      for (int i = 0; i < m_dmd.size(); i++) {
        if ((a_t >= m_intervals[i].first) && (a_t < m_intervals[i].second)) {
          return m_dmd[i]->predict(a_t);
        }
      }
      printf("ERROR in DMDROMObject::predict(): m_dmd is of size zero or interval not found!");
      return nullptr;
    }

    /*! save DMD object to file */
    virtual void save(const std::string& a_fname_root /*!< Filename root*/);

  protected:

    std::vector<CAROM::DMD*> m_dmd; /*!< Vector of DMD objects */
    std::vector<Interval> m_intervals; /*!< Time intervals for each DMD object */

    int m_rank;   /*!< MPI rank */
    int m_nproc;  /*!< Number of MPI ranks */

    int         m_vec_size; /*!< Local size of solution vector */
    double      m_dt;       /*!< Time step size */
    int         m_rdim;     /*!< Latent space dimension */

    int m_num_window_samples; /*!< Number of samples per DMD for time-windowing */

    int m_tic; /*! private ticker to count number of samples taken */
    int m_curr_win; /*! private index for current window */

    std::string m_dirname; /*!< Subdirectory where DMD objects are written to */

  private:
};

/*! \class libROMInterface
    \brief Class implementing interface with libROM
    This class implements the interface between HyPar
    and libROM. For details on libROM, see:
    https://github.com/LLNL/libROM.
*/

/*! \brief Class implementing interface with libROM
 
    This class implements the interface between HyPar
    and libROM. For details on libROM, see:
    https://github.com/LLNL/libROM.
*/
class libROMInterface
{
  public:

    /*! Default constructor */
    libROMInterface()
    {
      m_is_defined = false;
      m_rank = -1;
      m_nproc = -1;
      m_nsims = -1;
      m_vec_size.clear();
      m_rdim = -1;
      m_sampling_freq = 1;
      m_mode = "train";
      m_rom.clear();
      m_U.clear();
      m_train_wctime = 0;
      m_predict_wctime = 0;
      m_save_ROM = true;
    }

    /*! Constructor */
    libROMInterface(  void*   a_s, /*!< Array of simulation objects of type #SimulationObject */
                      int     a_nsims, /*!< number of simulation objects */
                      int     a_rank,  /*!< MPI rank of this process */
                      int     a_nproc, /*!< Number of MPI processes */
                      double  a_dt     /*!< Time step size */ )
    {
      define( a_s, a_nsims, a_rank, a_nproc, a_dt );
    }

    /*! Define the interface */
    void define( void*, int, int, int, double);

    /*! Destructor */
    ~libROMInterface()
    {
      if (m_mode == "train") {
        for (int ns = 0; ns < m_nsims; ns++) {
          delete m_rom[ns];
          delete m_U[ns];
        }
      }
    }

    /*! return the sampling frequency */
    inline int samplingFrequency() const { return m_sampling_freq; }
    /*! return the mode */
    inline const std::string& mode() const { return m_mode; }
    /*! return the ROM type */
    inline const std::string& type() const { return m_rom_type; }
    /*! return the training wall clock time */
    inline double trainWallclockTime() const { return m_train_wctime; }
    /*! return the prediction wall clock time */
    inline double predictWallclockTime() const { return m_predict_wctime; }

    /*! Take a sample for training */
    inline void takeSample( void*   a_s, /*!< Array of simulation objects of type #SimulationObject */
                            double  a_t /*!< Current simulation time */ )
    {
      copyFromHyPar( m_U, a_s );
      for (int ns = 0; ns < m_nsims; ns++) {
        m_rom[ns]->takeSample( m_U[ns], a_t );
      }
    }

    /*! Train the ROM object */
    inline void train()
    {
#ifndef serial
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      gettimeofday(&m_train_start, NULL);
      for (int ns = 0; ns < m_nsims; ns++) {
        if (m_nsims > 1) {
          if (!m_rank) {
            printf("Training ROM object for domain %d.\n", ns);
          }
        }
        m_rom[ns]->train();
      }
#ifndef serial
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      gettimeofday(&m_train_end, NULL);

      long long walltime;
      walltime = (  (m_train_end.tv_sec*1000000 + m_train_end.tv_usec)
                  - (m_train_start.tv_sec*1000000 + m_train_start.tv_usec) );
      m_train_wctime = (double) walltime / 1000000.0;
    }

    /*! Predict the solution at a given time */
    inline void predict(  void*  a_s, /*!< Array of simulation objects of type #SimulationObject */
                          double a_t  /*!< time at which to predict solution */)
    {
      m_predict_wctime = 0.0;
      for (int ns = 0; ns < m_nsims; ns++) {
#ifndef serial
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        gettimeofday(&m_predict_start, NULL);
        const CAROM::Vector* const u_predicted = m_rom[ns]->predict(a_t);
#ifndef serial
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        gettimeofday(&m_predict_end, NULL);
        copyToHyPar( u_predicted->getData(), a_s, ns );

        long long walltime;
        walltime = (  (m_predict_end.tv_sec*1000000 + m_predict_end.tv_usec)
                    - (m_predict_start.tv_sec*1000000 + m_predict_start.tv_usec) );
        m_predict_wctime += (double) walltime / 1000000.0;
      }
    }

    /*! Save ROM object to file */
    inline void saveROM(const std::string& a_fname_root=""/*!< filename root */)
    {
      if (m_save_ROM) {      
        for (int ns = 0; ns < m_nsims; ns++) {
          std::string fname_root = a_fname_root;
          if (m_nsims > 1) {
            char idx_string[_MAX_STRING_SIZE_];
            sprintf(idx_string, "sim%03d", ns);
            fname_root += std::string(idx_string);
          }
          if (!m_rank) {
            printf("libROMInterface::saveROM() - saving ROM objects.\n");
          }
          m_rom[ns]->save(fname_root);
        }
      }
      return;
    }

    /*! Copy HyPar solution to the work vectors m_U */
    void copyFromHyPar( std::vector<double*>&, void* );

    /*! Copy the work vectors m_U to HyPar */
    void copyToHyPar( std::vector<double*>&, void* );

    /*! Copy a vector to HyPar */
    void copyToHyPar( double*, void*, int );

  protected:

    bool m_is_defined;  /*!< Boolean to show if this object is defined */

    int m_rank;   /*!< MPI rank */
    int m_nproc;  /*!< Number of MPI ranks */
    int m_nsims;  /*!< Number of simulations */

    std::vector<int> m_vec_size; /*!< Solution vector size for each simulation*/

    int m_rdim; /*!< Reduced model dimensionality */
    int m_sampling_freq; /*!< Frequency at which to take samples */

    std::string m_mode; /*!< Mode: none, train */
    std::string m_rom_type; /*!< Type of ROM (eg: DMD) */

    std::vector<ROMObject*> m_rom; /*!< ROM objects */
    std::vector<double*> m_U; /*!< Work vectors */

    struct timeval m_train_start; /*<! Training start time */
    struct timeval m_train_end; /*<! Training end time */
    struct timeval m_predict_start; /*<! Prediction start time */
    struct timeval m_predict_end; /*<! Prediction end time */
    double m_train_wctime; /*!< Wallclock time for training */
    double m_predict_wctime; /*!< Wallclock time for prediction */

    bool m_save_ROM; /*!< Save ROM objects to file (default: yes) */

  private:

};

#endif

#endif
