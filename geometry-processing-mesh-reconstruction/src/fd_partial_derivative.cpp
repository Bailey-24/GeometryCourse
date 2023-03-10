#include "fd_partial_derivative.h"
#include <assert.h>

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
    int nx_D = nx, ny_D = ny, nz_D = nz;

    switch (dir)
    {
    case 0:
        nx_D--;
        break;
    case 1:
        ny_D--;
        break;
    case 2:
        nz_D--;
        break;

    default:
        assert(0);
    }
    // resize D so that each row corresponds to a node in the staggered grid
    // and each column corresponds to a node in the regular grid
    D.resize(nx_D * ny_D * nz_D, nx * ny * nz);

    // loop through each of the nodes in the staggered grid 
    // and assign -1 to the regular grid node "left"-neighbour
    // and 1 to the regular grid node "right"-neighbour
    for (int i=0; i<nx_D; i++)
    {
        for (int j=0; j<ny_D; j++)
        {
            for (int k=0; k<nz_D; k++)
            {
                int staggered_ijk = i + j * nx_D + k * ny_D * nx_D;
                int l_prev = i + j * nx + k * ny * nx;
                int l_curr;

                switch (dir)
                {
                case 0:
                    l_curr = (i + 1) + j * nx + k * ny * nx;
                    break;
                case 1:
                    l_curr = i + (j + 1) * nx + k * ny * nx;
                    break;
                case 2:
                    l_curr = i + j * nx + (k + 1) * ny * nx;
                    break;
                }
                D.insert(staggered_ijk, l_prev) = -1;
                D.insert(staggered_ijk, l_curr) = 1;
            }
        }
    }
}
