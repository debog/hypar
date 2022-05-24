#ifdef with_librom

/*! @file librom_interface.h
    @brief libROM interface class
    @author Debojyoti Ghosh
*/

#ifndef _LIBROM_INTERFACE_H_
#define _LIBROM_INTERFACE_H_

/*! Don't use any ROM */
#define _ROM_MODE_NONE_ "none"
/*! Training mode */
#define _ROM_MODE_TRAIN_ "train"
/*! Prediction mode */
#define _ROM_MODE_PREDICT_ "predict"

#include <string>
#include <vector>
#include <sys/time.h>

#include <linalg/Vector.h>
#include <rom_object.h>

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
    inline void takeSample( void* a_s,/*!< Array of simulation objects of type #SimulationObject */
                            const double a_t /*!< Current simulation time */ )
    {
      copyFromHyPar( m_U, a_s );
      gettimeofday(&m_train_start, NULL);
      for (int ns = 0; ns < m_nsims; ns++) {
        m_rom[ns]->takeSample( *(m_U[ns]), a_t );
      }
      gettimeofday(&m_train_end, NULL);

      long long walltime;
      walltime = (  (m_train_end.tv_sec*1000000 + m_train_end.tv_usec)
                  - (m_train_start.tv_sec*1000000 + m_train_start.tv_usec) );
      m_train_wctime += (double) walltime / 1000000.0;

#ifndef serial
      MPI_Allreduce(  MPI_IN_PLACE,
                      &m_train_wctime,
                      1,
                      MPI_DOUBLE,
                      MPI_MAX,
                      MPI_COMM_WORLD );
#endif
    }

    /*! Project initial solution for prediction */
    inline void projectInitialSolution( void* a_s  /*!< Array of simulation objects of 
                                                        type #SimulationObject */ )
    {
      copyFromHyPar( m_U, a_s );
      for (int ns = 0; ns < m_nsims; ns++) {
        m_rom[ns]->projectInitialSolution( *(m_U[ns]) );
      }
    }

    /*! Train the ROM object */
    inline void train()
    {
      gettimeofday(&m_train_start, NULL);
      for (int ns = 0; ns < m_nsims; ns++) {
        m_rom[ns]->train();
      }
      gettimeofday(&m_train_end, NULL);

      long long walltime;
      walltime = (  (m_train_end.tv_sec*1000000 + m_train_end.tv_usec)
                  - (m_train_start.tv_sec*1000000 + m_train_start.tv_usec) );
      m_train_wctime += (double) walltime / 1000000.0;
#ifndef serial
      MPI_Allreduce(  MPI_IN_PLACE,
                      &m_train_wctime,
                      1,
                      MPI_DOUBLE,
                      MPI_MAX,
                      MPI_COMM_WORLD );
#endif
    }

    /*! Predict the solution at a given time */
    inline void predict(  void*  a_s, /*!< Array of simulation objects of type #SimulationObject */
                          const double a_t  /*!< time at which to predict solution */) const
    {
      m_predict_wctime = 0.0;
      for (int ns = 0; ns < m_nsims; ns++) {
        gettimeofday(&m_predict_start, NULL);
        const CAROM::Vector* const u_predicted = m_rom[ns]->predict(a_t);
        gettimeofday(&m_predict_end, NULL);
        copyToHyPar( *u_predicted, a_s, ns );

        long long walltime;
        walltime = (  (m_predict_end.tv_sec*1000000 + m_predict_end.tv_usec)
                    - (m_predict_start.tv_sec*1000000 + m_predict_start.tv_usec) );
        m_predict_wctime += (double) walltime / 1000000.0;
      }
#ifndef serial
      MPI_Allreduce(  MPI_IN_PLACE,
                      &m_predict_wctime,
                      1,
                      MPI_DOUBLE,
                      MPI_MAX,
                      MPI_COMM_WORLD );
#endif
    }

    /*! Save ROM object to file */
    inline void saveROM(const std::string& a_fname_root=""/*!< filename root */) const
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

    /*! load ROM object from file */
    inline void loadROM(const std::string& a_fname_root=""/*!< filename root */)
    {
      for (int ns = 0; ns < m_nsims; ns++) {
        std::string fname_root = a_fname_root;
        if (m_nsims > 1) {
          char idx_string[_MAX_STRING_SIZE_];
          sprintf(idx_string, "sim%03d", ns);
          fname_root += std::string(idx_string);
        }
        if (!m_rank) {
          printf("libROMInterface::loadROM() - loading ROM objects.\n");
        }
        m_rom[ns]->load(fname_root);
      }
      return;
    }

    /*! Copy HyPar solution to the work vectors m_U */
    void copyFromHyPar( std::vector<CAROM::Vector*>&, void* );

    /*! Copy the work vectors m_U to HyPar */
    void copyToHyPar( const std::vector<CAROM::Vector*>&, void* ) const;

    /*! Copy a vector to HyPar */
    void copyToHyPar( const CAROM::Vector&, void*, int ) const;

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
    std::vector<CAROM::Vector*> m_U; /*!< Work vectors */

    mutable struct timeval m_train_start; /*<! Training start time */
    mutable struct timeval m_train_end; /*<! Training end time */
    mutable struct timeval m_predict_start; /*<! Prediction start time */
    mutable struct timeval m_predict_end; /*<! Prediction end time */
    mutable double m_train_wctime; /*!< Wallclock time for training */
    mutable double m_predict_wctime; /*!< Wallclock time for prediction */

    bool m_save_ROM; /*!< Save ROM objects to file (default: yes) */

  private:

};

#endif

#endif
