#include "../include/arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/arap_linear_block.h>
#include <igl/cotmatrix.h>

// Precompute data needed to efficiently conduct local-global iterations for an
// arap deformation. This includes the `data` struct employed by
// `igl::min_quad_with_fixed` to solve the global step" and constructing the
// bi-linear form `K` that mixes rotation degrees of freedom with unknown
// positions for preparing the covariance matrices of the local step and the
// linear term of the global step.
//
// Inputs:
//   V  #V by dim vertex positions at rest.
//   F  #F by simplex-size list of element indices
//   b  #b indices into V of handle vertices
// Outputs:
//   data  pre-factorized system matrix etc. (see `igl::min_quad_with_fixed`)
//   K  #R*dim by #V

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
    // get cotangent laplacian
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    // do precomputation. Factorize
    igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);

    // the cotangents of angles in cotangent laplacian
    // directly accessible by face idx
    Eigen::MatrixXd CE;
    igl::cotmatrix_entries(V, F, CE);

    // Initialize triplets used to fill the K matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(F.rows() * 3 * 6);

    // first loop through all faces
    for (int f = 0; f<F.rows(); f++)
    {
        // get the vertices of the current face
        // loop over each edge on this face
        for (int v = 0; v < 3; v++)
        {
            // for triangle meshes
            int i = F(f, (v + 1) % 3);
            int j = F(f, (v + 2) % 3);

            // compute the product of these into K
            Eigen::Vector3d eij = CE(f, v) * (V.row(i) - V.row(j));
            eij = eij / 3.0;

            // the alpha loop
            for (int k = 0; k < 3; k++)
            {
                // the beta loop
                for (int beta = 0; beta < 3; beta++)
                {
                    // multiple values are inserted into the same k
                    tripletList.push_back(T(i, 3*F(f, k) + beta, eij(beta)));
                    tripletList.push_back(T(j, 3*F(f, k) + beta, -eij(beta)));

                }
            }
        }
    }
    K.resize(V.rows(), 3*V.rows());
    K.setFromTriplets(tripletList.begin(), tripletList.end());

}
