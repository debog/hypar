/*! @file IBWriteBodySTL.c
    @author Debojyoti Ghosh
    @brief Writes a 3D body surface to STL file
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mpivars.h>
#include <immersedboundaries.h>


/*! Function to write a 3D surface from a STL file. See
    https://en.wikipedia.org/wiki/STL_(file_format)
    for details of the STL file format.
    \b Notes:
    + This function only writes ASCII STL files.

    It is assumed that all MPI ranks have the body data.
    The MPI rank specified as the input \a rank will actually
    write the file.
*/
int IBWriteBodySTL(
                    Body3D *body,     /*!< 3D body to write to file */
                    char   *filename, /*!< Filename */
                    void   *m,        /*!< MPI object of type #MPIVariables */
                    int    rank,      /*!< MPI rank that does the writing */
                    int    *stat      /*!< Status (0: success; non-0: failure) */
                  )
{
  MPIVariables  *mpi = (MPIVariables*) m;

  if (mpi->m_rank == rank) {
    FILE *out;
    out = fopen(filename,"w");
    if (!out) {
      fprintf(stderr,"Error in IBWriteBodySTL(): Unable to open file ");
      fprintf(stderr,"%s for writing on rank %d.\n",filename,mpi->m_rank);
      *stat = 1;
    } else {
      fprintf(out,"solid TEST\n");
      int n;
      for (n = 0; n < body->m_nfacets; n++) {
        fprintf(out,"  facet normal %1.16e %1.16e %1.16e\n",
                body->m_surface[n].m_nx,body->m_surface[n].m_ny,body->m_surface[n].m_nz);
        fprintf(out,"    outer loop\n");
        fprintf(out,"      vertex %1.16e %1.16e %1.16e\n",
                body->m_surface[n].m_x1,body->m_surface[n].m_y1,body->m_surface[n].m_z1);
        fprintf(out,"      vertex %1.16e %1.16e %1.16e\n",
                body->m_surface[n].m_x2,body->m_surface[n].m_y2,body->m_surface[n].m_z2);
        fprintf(out,"      vertex %1.16e %1.16e %1.16e\n",
                body->m_surface[n].m_x3,body->m_surface[n].m_y3,body->m_surface[n].m_z3);
        fprintf(out,"    endloop\n");
        fprintf(out,"  endfacet\n");
      }
      fprintf(out,"endsolid TEST\n");
      *stat = 0;
    }
  }
  return(0);
}

