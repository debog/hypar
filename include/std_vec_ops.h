/*! @file std_vec_ops.h
    @brief Contains some vector ops
    @author Debojyoti Ghosh
 */

#ifndef _STD_VEC_OPS_H_
#define _STD_VEC_OPS_H_

#include <math.h>
#include <vector>
#include <basic.h>

/*! Vector operations defined for the std::vector class */
namespace StdVecOps {

  /*! Sum of all elements */
  inline
  long sum( const std::vector<int>& a_iv /*!< Integer vector */ )
  {
    long retval(0);
    for (int i=0; i<a_iv.size(); i++) retval += a_iv[i];
    return retval;
  }

  /*! Product of all elements */
  inline
  long product( const std::vector<int>& a_iv /*!< Integer vector */ )
  {
    long retval(1);
    for (int i=0; i<a_iv.size(); i++) retval *= a_iv[i];
    return retval;
  }

  /*! Copy from a C-array */
  inline
  void copyFrom(  std::vector<int>& a_iv,       /*!< C++ vector of ints */
                  const int* const  a_iv_carr,  /*!< C array of ints */
                  int               a_n         /*!< size */
               )
  {
    a_iv.resize(a_n);
    for (int i=0; i<a_n; i++) a_iv[i] = a_iv_carr[i];
    return;
  }

  /*! Add constant to all elements of a vector */
  inline
  void add( std::vector<int>& a_iv, /*!< C++ vector of ints */
            const int         a_a   /*!< constant to add */ )
  {
    for (int i=0; i<a_iv.size(); i++) a_iv[i] += a_a;
    return;
  }

  /*! Copy from a C-array */
  inline
  void copyFrom(  std::vector<double>& a_iv,       /*!< C++ vector of doubles */
                  const double* const  a_iv_carr,  /*!< C array of doubles */
                  int                  a_n         /*!< size */
               )
  {
    a_iv.resize(a_n);
    for (int i=0; i<a_n; i++) a_iv[i] = a_iv_carr[i];
    return;
  }

  /*! Create a normal vector of doubles (of unit magnitude) from
   *  a C-array of integers */
  inline
  void createNormalVector(std::vector<double>&  a_normal_vec, /*!< Normal vector */
                          const int*            a_vec,        /*!< C-array of integers*/
                          const int             a_size        /*!< size of C-array */)
  {
    a_normal_vec = std::vector<double>(a_size, 0.0);
    double magn = 0.0;
    for (int i=0; i<a_size; i++) {
      magn += (double)(a_vec[i]*a_vec[i]);
    }
    magn = sqrt(magn);
    if (magn > _MACHINE_ZERO_) {
      for (int i=0; i<a_size; i++) {
        a_normal_vec[i] = a_vec[i] / magn;
      }
    }
    return;
  }

  /*! Create a normal vector of doubles (of unit magnitude) from
   *  a vector of integers */
  inline
  void createNormalVector(std::vector<double>&    a_normal_vec, /*!< Normal vector */
                          const std::vector<int>& a_vec         /*!< Integer vector */)
  {
    createNormalVector( a_normal_vec,
                        a_vec.data(),
                        a_vec.size() );
    return;
  }

  /*! Compute norm between two vectors */
  inline
  double compute2Norm(  const std::vector<double>& a_a, /*!< input vector */
                        const std::vector<double>& a_b  /*!< input vector */ )
  {
    double retval = 0.0;

    for (int i = 0; i < std::min(a_a.size(),a_b.size()); i++) {
      retval += ( (a_a[i]-a_b[i]) * (a_a[i]-a_b[i]) );
    }
    retval = sqrt(retval);

    return retval;
  }

}

#endif
