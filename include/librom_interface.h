#ifdef with_librom

/*! @file librom_interface.h
    @brief libROM interface class
    @author Debojyoti Ghosh
*/

#ifndef _LIBROM_INTERFACE_H_
#define _LIBROM_INTERFACE_H_

#include <string>
#include <vector>
#include <sys/time.h>
#include <linalg/Vector.h>
#include <algo/DMD.h>

/*! Filename for libROM-related inputs */
#define _LIBROM_INP_FNAME_ "librom.inp"

/*! DMD-type ROM */
#define _ROM_TYPE_DMD_ "DMD"

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
    DMDROMObject( const int     a_vec_size, /*! vector size */
                  const double  a_dt        /*! time step size */,
                  const int     a_rdim      /*! latent space dimension */ )
    {
      m_dmd = new CAROM::DMD( a_vec_size, a_dt );
      m_rdim = a_rdim;
    }

    /*! Destructor */
    virtual ~DMDROMObject()
    {
      delete m_dmd;
    }

    /*! take a sample (solution snapshot) */
    virtual void takeSample(  double* a_U, /*! solution vector */
                              const double a_time /*! sample time */ )
    {
      m_dmd->takeSample( a_U, a_time );
    }

    /*! train the DMD object */
    virtual void train()
    {
      m_dmd->train(m_rdim);
    }

    /*! compute prediction at given time */
    virtual 
    const CAROM::Vector* const predict(const double a_t /*!< time at which to predict solution */ )
    {
      return m_dmd->predict(a_t);
    }

  protected:

    CAROM::DMD *m_dmd;  /*! DMD object */
    int m_rdim;         /*! Latent space dimension */

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
      m_vec_size = -1;
      m_vec_offsets.clear();
      m_rdim = -1;
      m_sampling_freq = 1;
      m_mode = "train";
      m_rom = nullptr;
      m_U = nullptr;
      m_train_wctime = 0;
      m_predict_wctime = 0;
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
        delete m_rom;
        delete m_U;
      }
    }

    /*! return vector size */
    inline int vectorSize() const { return m_vec_size; }
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
      m_rom->takeSample( m_U, a_t );
    }

    /*! Train the ROM object */
    inline void train()
    {
#ifndef serial
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      gettimeofday(&m_train_start, NULL);
      m_rom->train();
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
#ifndef serial
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      gettimeofday(&m_predict_start, NULL);
      const CAROM::Vector* const u_predicted = m_rom->predict(a_t);
#ifndef serial
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      gettimeofday(&m_predict_end, NULL);
      copyToHyPar( u_predicted->getData(), a_s );

      long long walltime;
      walltime = (  (m_predict_end.tv_sec*1000000 + m_predict_end.tv_usec)
                  - (m_predict_start.tv_sec*1000000 + m_predict_start.tv_usec) );
      m_predict_wctime = (double) walltime / 1000000.0;
    }

    /*! Copy HyPar solution to the work vector m_U */
    void copyFromHyPar( double*, void* );

    /*! Copy the work vector m_U to HyPar */
    void copyToHyPar( double*, void* );

  protected:

    bool m_is_defined;  /*!< Boolean to show if this object is defined */

    int m_rank;   /*!< MPI rank */
    int m_nproc;  /*!< Number of MPI ranks */
    int m_nsims;  /*!< Number of simulations */

    int m_vec_size; /*!< Solution vector size */
    std::vector<int> m_vec_offsets; /*!< Offset for each simulation */

    int m_rdim; /*!< Reduced model dimensionality */
    int m_sampling_freq; /*!< Frequency at which to take samples */

    std::string m_mode; /*!< Mode: none, train */
    std::string m_rom_type; /*!< Type of ROM (eg: DMD) */

    ROMObject* m_rom; /*!< ROM object */
    double* m_U; /*!< Work vector */

    struct timeval m_train_start; /*<! Training start time */
    struct timeval m_train_end; /*<! Training end time */
    struct timeval m_predict_start; /*<! Prediction start time */
    struct timeval m_predict_end; /*<! Prediction end time */
    double m_train_wctime; /*!< Wallclock time for training */
    double m_predict_wctime; /*!< Wallclock time for prediction */

  private:

};

#endif

#endif
