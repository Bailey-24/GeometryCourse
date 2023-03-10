// #include "fd_interpolate.h"

// void fd_interpolate(
//   const int nx,
//   const int ny,
//   const int nz,
//   const double h,
//   const Eigen::RowVector3d & corner,
//   const Eigen::MatrixXd & P,
//   Eigen::SparseMatrix<double> & W)
// {
//     typedef Eigen::Triplet<double> T;
// 	std::vector<T> weight_list;
// 	weight_list.reserve(P.rows()*8*2);
//     for(int row_index = 0; row_index < P.rows(); row_index++)
//     {
//         // normalize, relative position to corner
//         double x0 = (P(row_index, 0) - corner(0))/h;
//         double y0 = (P(row_index, 1) - corner(1))/h;
//         double z0 = (P(row_index, 2) - corner(2))/h;

//         // next on grid position
//         int x1 = x0 + 1;
//         int y1 = y0 + 1;
//         int z1 = z0 + 1;

//         // relative position to next grid position
//         double ratio_x = (P(row_index, 0) - (x0*h + corner[0]))/h;
//         double ratio_y = (P(row_index, 1) - (y0*h + corner[1]))/h;
//         double ratio_z = (P(row_index, 2) - (z0*h + corner[2]))/h;

//         // calculate weight for each of eight surrounding point
//         weight_list.push_back(Eigen::Triplet<double>(row_index, (x0 + y0*nx + z0*ny*nx), (1 - ratio_x) * (1 - ratio_y) * (1 - ratio_z)));
//         weight_list.push_back(Eigen::Triplet<double>(row_index, (x1 + y0*nx + z0*ny*nx), (    ratio_x) * (1 - ratio_y) * (1 - ratio_z)));
//         weight_list.push_back(Eigen::Triplet<double>(row_index, (x0 + y1*nx + z0*ny*nx), (1 - ratio_x) * (    ratio_y) * (1 - ratio_z)));
//         weight_list.push_back(Eigen::Triplet<double>(row_index, (x1 + y1*nx + z0*ny*nx), (    ratio_x) * (    ratio_y) * (1 - ratio_z)));
//         weight_list.push_back(Eigen::Triplet<double>(row_index, (x0 + y0*nx + z1*ny*nx), (1 - ratio_x) * (1 - ratio_y) * (    ratio_z)));
//         weight_list.push_back(Eigen::Triplet<double>(row_index, (x1 + y0*nx + z1*ny*nx), (    ratio_x) * (1 - ratio_y) * (    ratio_z)));
//         weight_list.push_back(Eigen::Triplet<double>(row_index, (x0 + y1*nx + z1*ny*nx), (1 - ratio_x) * (    ratio_y) * (    ratio_z)));
//         weight_list.push_back(Eigen::Triplet<double>(row_index, (x1 + y1*nx + z1*ny*nx), (    ratio_x) * (    ratio_y) * (    ratio_z)));
//     }
//     W.resize(P.rows(), nx*ny*nz);
//     W.setFromweight_list(weight_list.begin(), weight_list.end());

// }


#include "fd_interpolate.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
	using namespace std;
	typedef Eigen::Triplet<double> T;
	vector<T> weight_list;
	weight_list.reserve(P.rows()*8*2);
	for (int row = 0; row < P.rows(); row++) {
		double dx = P(row, 0) - corner(0);
		double dy = P(row, 1) - corner(1);
		double dz = P(row, 2) - corner(2);
		int i = int(floor(dx / h));
		double x_value = fmod(dx, h) / h;
		int j = int(floor(dy / h));
		double y_value = fmod(dy, h) / h;
		int k = int(floor(dz / h));
		double z_value = fmod(dz, h) / h;

		// compute the value for the 8 neighboring nodes
		weight_list.push_back(T(row, i + nx * (j + k * ny), 
			(1.0 - x_value) * (1.0 - y_value) * (1.0 - z_value)));
		weight_list.push_back(T(row, (i + 1) + nx * (j + k * ny),
			x_value * (1.0 - y_value) * (1.0 - z_value)));
		weight_list.push_back(T(row, i + nx * ((j + 1) + k * ny),
			(1.0 - x_value) * y_value * (1.0 - z_value)));
		weight_list.push_back(T(row, (i + 1) + nx * ((j + 1) + k * ny),
			x_value * y_value * (1.0 - z_value)));
		weight_list.push_back(T(row, i + nx * (j + (k + 1) * ny),
			(1.0 - x_value) * (1.0 - y_value) * z_value));
		weight_list.push_back(T(row, (i + 1) + nx * (j + (k + 1) * ny),
			x_value * (1.0 - y_value) * z_value));
		weight_list.push_back(T(row, i + nx * ((j + 1) + (k + 1) * ny),
			(1.0 - x_value) * y_value * z_value));
		weight_list.push_back(T(row, (i + 1) + nx * ((j + 1) + (k + 1) * ny),
			x_value * y_value * z_value));
	}
	W.resize(P.rows(), nx*ny*nz);
	W.setFromTriplets(weight_list.begin(), weight_list.end());
}