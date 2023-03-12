#include "../include/arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

// Given precomputed data (`data` and `K`), handle _positions_ `bc` and current
// positions of all vertices `U`, conduct a _single_ iteration of the
// local-global solver for minimizing the as-rigid-as-possible energy. Output
// the _positions_ of all vertices of the mesh (by overwriting `U`).
//
// Inputs:
//   data  pre-factorized system matrix etc. (see `arap_precompute`
//   K  pre-constructed bi-linear term of energy combining rotations and
//     positions
//   U  #V by dim list of current (not necessarily at rest) positions (see output)
//   bc  #b by dim list of positions for each handle vertex
// Outputs:
//   U  #V by dim list of new positions (see input)

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
    // local
    // compute the matrix of stacked weighted covaiane matrices
    Eigen::MatrixXd C, R, Beq;
    C = K.transpose() * U;

    // compute the closest rotation matrix Rk for each Ck 
    // and stack them
    R.resize(data.n * 3, 3);
    for (int i = 0; i < data.n; i++)
    {
        Eigen::Matrix3d Rk, Ck;
        Ck = C.block(i*3, 0, 3, 3);
        igl::polar_svd3x3(Ck, Rk);
        R.block(i*3, 0, 3, 3) = Rk;
    }

    // global
    // solve the minimization problem with bc constraints using precomputed data
    igl::min_quad_with_fixed_solve(data, K*R, bc, Beq, U);
    
}
