#include "vector_area_matrix.h"
#include <igl/boundary_loop.h>
#include <iostream>
// Constructs the symmetric area matrix A, s.t.  [V.col(0)' V.col(1)'] * A *
// [V.col(0); V.col(1)] is the **vector area** of the mesh (V,F).
//
// Inputs:
//   F  #F by 3 list of mesh faces (must be triangles)
// Outputs:
//   A  #Vx2 by #Vx2 sparse area matrix
//
void vector_area_matrix(
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double>& A)
{

    // Convert asymetric Ã -> symmetric A
    // A = 1/2 (Ã + Ã^T)
    //  where Ã is asymetric area matrix
    //  s.t. for each (i,j) \in bnd(S),
    //      Ã[i, j+n] = 1
    //      Ã[i+n, j] = 1

    int V_size = F.maxCoeff()+1;
    A.resize(V_size*2,V_size*2);

    // Determine the longest boundary loop
    std::vector<std::vector<int>> bnds;
    igl::boundary_loop(F, bnds);

    // Each edge (i, j) in the loop contributes 1/2 (x_i y_j - x_j y_i) to the vector area
    std::vector<Eigen::Triplet<double>> entries;
    int v1, v2;
    for (std::vector<int> bnd : bnds) {
        for (int i = 0; i < bnd.size(); i++) {
        v1 = bnd[i];
        v2 = bnd[(i + 1) % bnd.size()];

        // ensure symmetry
        // A
        entries.push_back(Eigen::Triplet<double>(v1, v2 + V_size, 1));
        entries.push_back(Eigen::Triplet<double>(v2, v1 + V_size, -1));

        // A^T
        entries.push_back(Eigen::Triplet<double>(v2 + V_size, v1, 1));
        entries.push_back(Eigen::Triplet<double>(v1 + V_size, v2, -1));
        }
    }

    A.setFromTriplets(entries.begin(), entries.end());
    A *= 0.5;
}

