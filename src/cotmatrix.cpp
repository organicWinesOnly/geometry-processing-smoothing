#include "cotmatrix.h"
#include "igl/doublearea.h"
#include <cmath>
#include <vector>


using namespace std;
// cotangent
// Calculate the cotangent of an angle using the length of the triangle only.
// Inputs: 
//     a, b, c lengths of the triangle with a being the length opposite of the
//             angle of interest.
//     d_area  double area of triangle
double cotangent(
    const double & a,
    const double & b,
    const double & c, 
    const double & d_area)
{
  double numerator, denom;
  numerator = pow(b, 2) + pow(c, 2) - pow(a, 2);
  denom = 2.0 * d_area;
  return numerator / denom;
}


typedef Eigen::Triplet<double> T;
void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Compute the areas of the triangles
  Eigen::VectorXd dArea;
  igl::doublearea(l, dArea);

  int number_vertices = F.maxCoeff() + 1;
  std::vector<T> triplet_list;
  triplet_list.reserve(4 * 3 * number_vertices);
  for (int q = 0; q < F.rows(); q++)
  {
    double a, b, c;
    a = l(q, 0);
    b = l(q, 1);
    c = l(q, 2);
    vector<double> sides = {a, b, c};

    for ( int r = 0; r < 3; r++)
    {
      int i = F(q, r);
      int j = F(q, (r+1) % 3);
      double weight =  0.5 * cotangent(sides[(r+2) % 3], sides[(r+1) % 3],
                                       sides[r], dArea(q));
      // two insertions because L is symmetric
      triplet_list.push_back(T(i, j, weight));
      triplet_list.push_back(T(j, i, weight));
      triplet_list.push_back(T(i, i, -1.0 * weight));
      triplet_list.push_back(T(j, j, -1.0 * weight));
    }
  }
  L.resize(number_vertices, number_vertices);
  L.setFromTriplets(triplet_list.begin(), triplet_list.end());
}

