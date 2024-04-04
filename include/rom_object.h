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
    virtual void projectInitialSolution( CAROM::Vector& ) = 0;
    /*! take a sample (solution snapshot) */
    virtual void takeSample( const CAROM::Vector&, const double ) = 0;
    /*! train the ROM object */
    virtual void train() = 0;
    /*! compute prediction at given time */
    virtual const CAROM::Vector* const predict( const double ) const = 0;
    /*! save ROM object to file */
    virtual void save( const std::string& ) const = 0;
    /*! load ROM object from file */
    virtual void load( const std::string& ) = 0;

  protected:

  private:

};

#endif

#endif
