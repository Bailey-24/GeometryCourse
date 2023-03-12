#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>

// Given precomputation data and a list of handle _displacements_ determine
// _displacements_ for all vertices in the mesh.
//
// Inputs:
//   data  pre-factorized system matrix etc. (see `biharmonic_precompute` Q_uu
//   bc  #b by 3 list of displacements for each handle vertex
// Outputs:
//   D  #V by 3 list of displacements

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
    // solove the minization problem with bc constraints
    // using the precompute data
    Eigen::VectorXd X, B;
    B = Eigen::VectorXd::Zero(data.n);
    D.resize(data.n, 3);

    igl::min_quad_with_fixed_solve(data, B, bc, X, D);
}
