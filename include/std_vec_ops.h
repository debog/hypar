/*! @file std_vec_ops.h
    @brief Contains some vector ops
    @author Debojyoti Ghosh
 */

#ifndef _STD_VEC_OPS_H_
#define _STD_VEC_OPS_H_

#include <vector>
#include <basic.h>
#include <mathfunctions.h>

/*! Vector operations defined for the std::vector class */
namespace StdVecOps {

  /*! Create a normal vector of doubles (of unit magnitude) from
   *  a C-array of integers */
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
  void createNormalVector(std::vector<double>&    a_normal_vec, /*!< Normal vector */
                          const std::vector<int>& a_vec         /*!< Integer vector */)
  {
    createNormalVector( a_normal_vec,
                        a_vec.data(),
                        a_vec.size() );
    return;
  }
                          
  /*! Compute norm between two vectors */
  double compute2Norm(  const std::vector<double>& a_a, /*!< input vector */
                        const std::vector<double>& a_b  /*!< input vector */ )
  {
    double retval = 0.0;

    for (int i = 0; i < min(a_a.size(),a_b.size()); i++) {
      retval += ( (a_a[i]-a_b[i]) * (a_a[i]-a_b[i]) );
    }
    retval = sqrt(retval);

    return retval;
  }
                          
}

#endif
