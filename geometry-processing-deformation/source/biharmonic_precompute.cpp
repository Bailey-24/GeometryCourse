#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

// Precompute data needed to efficiently solve for a biharmonic deformation
// given a mesh with vertices `V` and faces `F` and a list of selected vertices
// as indices `b` into `V`. The output should be a prefacorized system using
// the `data` struct employed by `igl::min_quad_with_fixed`.
//
// Inputs:
//   V  #V by dim vertex positions
//   F  #F by simplex-size list of element indices
//   b  #b indices into V of handle vertices
// Outputs:
//   data  pre-factorized system matrix etc. (see `igl::min_quad_with_fixed`)

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
    // mass and cot matrix
    Eigen::SparseMatrix<double> laplacian;
    igl::cotmatrix(V, F, laplacian);
    
    // Q = L^T * M^-1 * L
    Eigen::SparseMatrix<double> mass, mass_inv, Q;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);
    igl::invert_diag(mass, mass_inv);
    Q = laplacian.transpose() * mass_inv * laplacian;

    // precompute data
    igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), false, data);
}
