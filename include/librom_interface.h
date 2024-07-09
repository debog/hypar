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
/*! OFFLINE mode */
#define _ROM_MODE_OFFLINE_ "offline"
/*! Merge mode */
#define _ROM_MODE_MERGE_ "merge"
/*! ONLINE mode */
#define _ROM_MODE_ONLINE_ "online"
/*! Prediction mode */
#define _ROM_MODE_PREDICT_ "predict"
/*! Initial guess mode - use ROM prediction as initial guess
 *  in implicit time integration */
#define _ROM_MODE_INITIAL_GUESS_ "initial_guess"

/*! Monolithic ROM for all components */
#define _ROM_COMP_MODE_MONOLITHIC_ "monolithic"
/*! Separate ROM objects for each component */
#define _ROM_COMP_MODE_COMPONENTWISE_ "component-wise"

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
      m_mode = _ROM_MODE_TRAIN_;
      m_comp_mode = "monolithic";
      m_ncomps.clear();
      m_rom.clear();
      m_U.clear();
      m_train_wctime = 0;
      m_predict_wctime = 0;
      m_save_ROM = true;
      m_parametric_id = 0;
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
      for (int i = 0; i < m_rom.size(); i++) delete m_rom[i];
      m_rom.clear();
      for (int i = 0; i < m_U.size(); i++) delete m_U[i];
      m_U.clear();
      for (int i = 0; i < m_U_stages.size(); i++) delete m_U_stages[i];
      m_U_stages.clear();
    }

    /*! return the sampling frequency */
    inline int samplingFrequency() const { return m_sampling_freq; }
    /*! return the mode */
    inline const std::string& mode() const { return m_mode; }
    /*! return the component mode */
    inline const std::string& componentMode() const { return m_comp_mode; }
    /*! return the ROM type */
    inline const std::string& type() const { return m_rom_type; }
    /*! return the training wall clock time */
    inline double trainWallclockTime() const { return m_train_wctime; }
    /*! return the prediction wall clock time */
    inline double predictWallclockTime() const { return m_predict_wctime; }
    /*! return the file op iter */
    inline int file_op_iter() const { return m_file_op_iter; }

    /*! Take a sample for training */
    void takeSample( void* a_s, const double a_t, const int a_nstep, void* ts);

    /*! Project initial solution for prediction */
    void projectInitialSolution( void* a_s );

    /*! Train the ROM object */
    void train(void* a_s);

    /*! Predict the solution at a given time */
    void predict(void*  a_s, const double a_t) const;

    /*! Save ROM object to file */
    void saveROM(const std::string& a_fname_root=""/*!< filename root */) const;

    /*! load ROM object from file */
    void loadROM(const std::string& a_fname_root=""/*!< filename root */);

    /*! Copy HyPar solution to the work vectors m_U */
    void copyFromHyPar( std::vector<CAROM::Vector*>&, void* );

    /*! Copy HyPar multi-stages to the work vectors m_U_stages */
    void copyFromHyPar_stages( std::vector<CAROM::Vector*>&, void*, void* );

    /*! Copy a vector to HyPar */
    void copyToHyPar( const CAROM::Vector&, void*, int ) const;

    /*! Copy a vector to HyPar */
    void copyToHyPar( const CAROM::Vector&, void*, int, int ) const;

    /*! Write Snapshot matrix */
    void writeSnapshot(void* a_s);

    /*! Merge stage */
    void merge(void* a_s);

    /*! Online stage */
    void online(void* a_s);

    /*! Write final solution */
    void writeFinal(void* a_s);

    /*! Read final solution */
    void readFinal(void* a_s);

    /* Clean up allocation related to ROMObject */
    void cleanup(void* a_s);

  protected:

    bool m_is_defined;  /*!< Boolean to show if this object is defined */

    int m_rank;   /*!< MPI rank */
    int m_nproc;  /*!< Number of MPI ranks */
    int m_nsims;  /*!< Number of simulations */

    std::vector<int> m_vec_size; /*!< Solution vector size for each simulation*/

    int m_rdim; /*!< Reduced model dimensionality */
    int m_sampling_freq; /*!< Frequency at which to take samples */
    int m_file_op_iter; /*!< frequency (iterations) of writing ROM solution to file */
    int m_parametric_id; /*!< frequency (iterations) of writing ROM solution to file */

    std::string m_mode; /*!< Mode: none, train, predict */
    std::string m_comp_mode; /*!< Component mode: monolithic, component-wise */
    std::string m_rom_type; /*!< Type of ROM (eg: DMD) */

    std::vector<int> m_ncomps; /*!< Number of vector components for each simulation */

    std::vector<ROMObject*> m_rom; /*!< ROM objects */
    std::vector<CAROM::Vector*> m_U; /*!< Work vectors */
    std::vector<CAROM::Vector*> m_U_stages; /*!< Work vectors */

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
