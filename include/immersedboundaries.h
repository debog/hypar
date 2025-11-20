/*! @file immersedboundaries.h
    @brief Structures and function definitions for immersed boundaries
    @author Debojyoti Ghosh
*/

#ifndef _IB_H_
#define _IB_H_

/*! Immersed boundaries are implemented only for 3-dimensional simulations. */
#define _IB_NDIMS_ 3
/*! Number of grid points surrounding a random point in space, i.e., number of
    corners of a cube. */
#define _IB_NNODES_ 8

/*! "Pseudo-2D" simulation in the x-y plane */
#define _IB_XY_ "2d (xy)"
/*! "Pseudo-2D" simulation in the x-z plane */
#define _IB_XZ_ "2d (xz)"
/*! "Pseudo-2D" simulation in the y-z plane */
#define _IB_YZ_ "2d (yz)"
/*! 3D simulation */
#define _IB_3D_ "3d"

#include <basic.h>

/*! \def Facet3D
    \brief Structure defining a facet.

    A "facet" is the basic unit of a tessellated surface in 3D: It is a triangle
    in an unstructured triangulated surface, defined by its vertices and surface
    normal. \sa https://en.wikipedia.org/wiki/STL_(file_format)
*/
/*!\brief Structure defining a facet.

    A "facet" is the basic unit of a tessellated surface in 3D: It is a triangle
    in an unstructured triangulated surface, defined by its vertices and surface
    normal. \sa https://en.wikipedia.org/wiki/STL_(file_format)
*/
typedef struct _facet_3d_{
  double m_x1, /*!< x-coordinate of vertex 1 */
         m_x2, /*!< x-coordinate of vertex 2 */
         m_x3, /*!< x-coordinate of vertex 3 */
         m_y1, /*!< y-coordinate of vertex 1 */
         m_y2, /*!< y-coordinate of vertex 2 */
         m_y3, /*!< y-coordinate of vertex 3 */
         m_z1, /*!< z-coordinate of vertex 1 */
         m_z2, /*!< z-coordinate of vertex 2 */
         m_z3, /*!< z-coordinate of vertex 3 */
         m_nx, /*!< x-component of surface normal */
         m_ny, /*!< y-component of surface normal */
         m_nz; /*!< z-component of surface normal */
} Facet3D;

/*! \def FacetMap
    \brief Structure defining a facet map.

    A facet map contains information for facets that lie
    within the local computational domain of this MPI
    rank.
*/
/*! \brief Structure defining a facet map.

    A facet map contains information for facets that lie
    within the local computational domain of this MPI
    rank.
*/
typedef struct _facet_map_{
  Facet3D   *m_facet; /*!< pointer to the facet */
  int       m_index;  /*!< index of this facet in the array #Body3D::m_surface */
  int       m_interp_nodes    [_IB_NNODES_];/*!< indices of grid points surrounding the facet centroid */
  double    m_interp_coeffs   [_IB_NNODES_];/*!< interpolation coefficients corresponding to #FacetMap::m_interp_nodes */
  int       m_interp_nodes_ns [_IB_NNODES_];/*!< indices of grid points surrounding the "near-surface" point near the centroid */
  double    m_interp_coeffs_ns[_IB_NNODES_];/*!< interpolation coefficients corresponding to #FacetMap::m_interp_nodes_ns */

  double    m_xc, /*!< x-coordinate of centroid of #FacetMap::m_facet */
            m_yc, /*!< y-coordinate of centroid of #FacetMap::m_facet */
            m_zc, /*!< z-coordinate of centroid of #FacetMap::m_facet */
            m_xns,/*!< x-coordinate of "near surface" point of #FacetMap::m_facet */
            m_yns,/*!< y-coordinate of "near surface" point of #FacetMap::m_facet */
            m_zns;/*!< z-coordinate of "near surface" point of #FacetMap::m_facet */

  double    m_dx, /*!< #FacetMap::m_xns - #FacetMap::m_xc */
            m_dy, /*!< #FacetMap::m_yns - #FacetMap::m_yc */
            m_dz; /*!< #FacetMap::m_zns - #FacetMap::m_zc */
} FacetMap;

/*! \def Body3D
    \brief Structure defining a body.

    A 3D body whose surface is represented as a collection
    of faces of type #Facet3D.
    \sa https://en.wikipedia.org/wiki/STL_(file_format)
*/
/*! \brief Structure defining a body.

    A 3D body whose surface is represented as a collection
    of faces of type #Facet3D.
    \sa https://en.wikipedia.org/wiki/STL_(file_format)
*/
typedef struct _body_3d_{
  int     m_nfacets;    /*!< number of surface facets   */
  Facet3D *m_surface;   /*!< array of surface facets    */
  /* coordinates of bounding box */
  double m_xmin, /*!< x-coordinate of lower end of bounding box */
         m_xmax, /*!< x-coordinate of higher end of bounding box */
         m_ymin, /*!< y-coordinate of lower end of bounding box */
         m_ymax, /*!< y-coordinate of higher end of bounding box */
         m_zmin, /*!< z-coordinate of lower end of bounding box */
         m_zmax; /*!< z-coordinate of higher end of bounding box */
} Body3D;

/*! \def IBNode
    \brief Structure defining an immersed boundary node

    An immersed boundary node is a grid point inside the immersed
    body but within stencil-width-distance of a point outside
    the body.
*/
/*! \brief Structure defining an immersed boundary node

    An immersed boundary node is a grid point inside the immersed
    body but within stencil-width-distance of a point outside
    the body.
*/
typedef struct _boundary_node_{
  int     m_i,      /*!< grid index along x */
          m_j,      /*!< grid index along y */
          m_k,      /*!< grid index along z */
          m_p;      /*!< array index of this point */
  double  m_x,      /*!< x-coordinate of the boundary node */
          m_y,      /*!< y-coordinate of the boundary node */
          m_z;      /*!< z-coordinate of the boundary node */
  Facet3D *m_face;  /*!< the nearest facet of #ImmersedBoundary::m_body */

  int     m_interp_nodes[_IB_NNODES_];   /*!< indices of interior nodes from which to extrapolate */
  double  m_interp_coeffs[_IB_NNODES_],  /*!< interpolation coefficients corresponding to #IBNode::m_interp_nodes */
          m_interp_node_distance,        /*!< Distance from the immersed body surface to the interior node from which to extrapolate */
          m_surface_distance;            /*!< Distance from this node to the immersed body surface */

} IBNode;

/*! \def ImmersedBoundary
    \brief Structure containing variables for immersed boundary implementation.

    This structure contains all the variables needed to implement immersed
    boundaries.
*/
/*! \brief Structure containing variables for immersed boundary implementation.

    This structure contains all the variables needed to implement immersed
    boundaries.
*/
typedef struct immersed_boundary{
  Body3D    *m_body;      /*!< immersed body */
  IBNode    *m_boundary;  /*!< immersed boundary nodes */
  FacetMap  *m_fmap;      /*!< list of "local" facets */

  double  m_tolerance; /*!< zero tolerance */
  double  m_delta;     /*!< small number */
  int     m_itr_max;   /*!< maximum intersections in ray-tracing method */
  int     m_n_boundary_nodes; /*!< number of immersed boundary nodes */
  int     m_nfacets_local;    /*!< number of "local" facets */

  char    m_mode[_MAX_STRING_SIZE_]; /*!< identifies if the simulation is 2D along a plane
                                        or truly 3D. \sa IBIdentifyMode() */
} ImmersedBoundary;


int IBReadBodySTL (Body3D**,char*,void*,int*);
int IBWriteBodySTL(Body3D*,char*,void*,int,int*);

int IBCleanup           (void*);
int IBComputeBoundingBox(Body3D*);
int IBCreateFacetMapping(void*,void*,double*,int*,int);
int IBIdentifyBody      (void*,int*,int*,int,void*,double*,double*);
int IBIdentifyBoundary  (void*,void*,int*,int,double*);
int IBIdentifyMode      (double*,int*,void*);
int IBNearestFacetNormal(void*,void*,double*,double,int*,int);
int IBInterpCoeffs      (void*,void*,double*,int*,int,double*);

int IBAssembleGlobalFacetData(void*,void*,const double* const, double** const,int);

int IBComputeNormalGradient(void*,void*,const double* const, int, double** const);
int IBComputeFacetVar(void*,void*,const double* const, int, double** const);

#endif
