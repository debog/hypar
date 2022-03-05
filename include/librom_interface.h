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

/*! \class libROMInterface
    \brief Class implementing interface with libROM
    This class implements the interface between HyPar
    and libROM. For details on libROM, see:
    https://github.com/LLNL/libROM.
*/

/*! \brief Class implementing interface with libROM
 *
    This class implements the interface between HyPar
    and libROM. For details on libROM, see:
    https://github.com/LLNL/libROM.
*/
class libROMInterface
{
  public:

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

    inline int vectorSize() const { return m_vec_size; }
    inline int samplingFrequency() const { return m_sampling_freq; }
    inline const std::string& mode() const { return m_mode; }
    inline double trainWallclockTime() const { return m_train_wctime; }
    inline double predictWallclockTime() const { return m_predict_wctime; }

    inline void takeSample( void*   a_s, /*!< Array of simulation objects of type #SimulationObject */
                            double  a_t /*!< Current simulation time */ )
    {
      copyFromHyPar( m_U, a_s );
      m_rom->takeSample( m_U, a_t );
    }

    inline void train()
    {
#ifndef serial
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      gettimeofday(&m_train_start, NULL);
      m_rom->train(m_rdim);
#ifndef serial
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      gettimeofday(&m_train_end, NULL);

      long long walltime;
      walltime = (  (m_train_end.tv_sec*1000000 + m_train_end.tv_usec)
                  - (m_train_start.tv_sec*1000000 + m_train_start.tv_usec) );
      m_train_wctime = (double) walltime / 1000000.0;
    }

    inline void predict(  void*  a_s, /*!< Array of simulation objects of type #SimulationObject */
                          double a_t  /*!< time at which to predict solution */)
    {
#ifndef serial
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      gettimeofday(&m_predict_start, NULL);
      CAROM::Vector* u_predicted = m_rom->predict(a_t);
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

    CAROM::DMD *m_rom; /*!< ROM object */

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
