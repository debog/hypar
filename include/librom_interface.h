#ifdef with_librom

/*! @file librom_interface.h
    @brief libROM interface class
    @author Debojyoti Ghosh
*/

#ifndef _LIBROM_INTERFACE_H_
#define _LIBROM_INTERFACE_H_

#include <vector>
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
      m_rom = nullptr;
      m_U = nullptr;
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
      delete m_rom;
      delete m_U;
    }

    inline int vectorSize() const { return m_vec_size; }
    inline int samplingFrequency() const { return m_sampling_freq; }

    inline void takeSample( void*   a_s, /*!< Array of simulation objects of type #SimulationObject */
                            double  a_t /*!< Current simulation time */ )
    {
      copyFromHyPar( m_U, a_s );
      m_rom->takeSample( m_U, a_t );
    }

    inline void train()
    {
      m_rom->train(m_rdim);
    }

    inline void predict(  void*  a_s, /*!< Array of simulation objects of type #SimulationObject */
                          double a_t  /*!< time at which to predict solution */)
    {
      CAROM::Vector* u_predicted = m_rom->predict(a_t);
      copyToHyPar( u_predicted->getData(), a_s );
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

    CAROM::DMD *m_rom; /*!< ROM object */

    double* m_U; /*!< Work vector */

  private:

};

#endif

#endif
