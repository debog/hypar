#ifdef with_librom

/*! @file rom_object.h
    @brief Base class for a libROM object
    @author Debojyoti Ghosh
*/

#ifndef _ROM_OBJECT_H_
#define _ROM_OBJECT_H_

/*! No ROM type */
#define _ROM_TYPE_NONE_ "none"

/*! Filename for libROM-related inputs */
#define _LIBROM_INP_FNAME_ "librom.inp"

#include <string>
#include <linalg/Vector.h>

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

    /*! virtual destructor */
    virtual ~ROMObject() { }

    /*! Project initial solution for prediction */
    virtual void projectInitialSolution( CAROM::Vector&, void* ) = 0;
    /*! take a sample (solution snapshot) */
    virtual void takeSample( const CAROM::Vector&, const std::vector<CAROM::Vector*>&, const double, const int, void* ) = 0;
    /*! train the ROM object */
    virtual void train(void*) = 0;
    /*! compute prediction at given time */
    virtual const CAROM::Vector* predict( const double, void* ) = 0;
    /*! save ROM object to file */
    virtual void save( const std::string& ) const = 0;
    /*! load ROM object from file */
    virtual void load( const std::string& ) = 0;
    /*! Write Snapshot matrix */
    virtual void writeSnapshot(void*) = 0;
    /*! Merge stage */
    virtual void merge(void*) = 0;
    /*! Online stage */
    virtual void online(void*) = 0;
    /* Clean up allocation related to ROMObject */
    virtual void cleanup(void*) = 0;

  protected:

  private:

};

#endif

#endif
