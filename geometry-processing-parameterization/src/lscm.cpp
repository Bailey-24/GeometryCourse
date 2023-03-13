#include "lscm.h"
#include "vector_area_matrix.h"

#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>
#include <Eigen/SVD>
#include <iostream>
// Given a 3D mesh (`V`,`F`) with boundary compute a 2D parameterization that
// minimizes the "least squares conformal" energy:
//
// \\[
// ∫_Ω ‖ ∇v - (∇u)^⊥ ‖² dA,
// \\]
//
// where u and v are the unknown (output) coordinates in the parametric domain
// `U`.
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of triangle indices into V
// Outputs:
//   U  #U by 2 list of mesh UV parameterization coordinates
//
void lscm(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
    // Solve optimization as a generalized Eigen value problem
    //      min_U U'QU  subject to  U'BU = 1

    // Angle distortion setup, construct matrix Q
    Eigen::SparseMatrix<double> lapacian, diag_lapacian, vector_area;
    igl::cotmatrix(V, F, lapacian);
    igl::repdiag(lapacian, 2, diag_lapacian);
    vector_area_matrix(F, vector_area);
    Eigen::SparseMatrix<double> Q = diag_lapacian - vector_area;

    // Free boundary constraint setup, construct matrix B
    Eigen::SparseMatrix<double> mass, B;
    igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, mass);
    igl::repdiag(mass, 2, B);

    // solve the generalized Eigen value problem
    Eigen::MatrixXd eigen_vectors;
    Eigen::VectorXd eigen_values;
    igl::eigs(Q, B, 3, igl::EigsType::EIGS_TYPE_SM, eigen_vectors, eigen_values);

    // extract the Fiedler vector: 3rd -- first 2 trivial
    // Somehow, first 2 eigen value ~e^{-13}
    //      3, 4, ... looks more reasonable
    int V_size = V.rows();
    U.resize(V_size, 2);
    U.setZero();

    U.col(0) = eigen_vectors.col(2).head(V_size);
    U.col(1) = eigen_vectors.col(2).tail(V_size);

    // Find Canonical rotation V, using PCA with SVD on U^T U
    //      Right singular vectors V are eigenvectors of U^T U
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(U.transpose() * U, Eigen::ComputeThinU | Eigen::ComputeThinV);
    U = U * svd.matrixU();
}
