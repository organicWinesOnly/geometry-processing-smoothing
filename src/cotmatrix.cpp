#include "cotmatrix.h"
#include "igl/doublearea.h"
#include <cmath>
#include <map>
#include <utility>
#include <algorithm>


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
  denom = 2 * d_area;
  return numerator / denom;
}


typedef Eigen::Triplet<double> T;
typedef pair<int, int> P;
void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Compute the areas of the triangles
  Eigen::VectorXd dArea;
  igl::doublearea(l, dArea);

  map<P, double> laplacian_values;
  double sum_of_weights = 0;
  // for each face, compute the 3 cotangent values
  // for each race create a pair of idicies and put them in a map object; 
  for (int q = 0; q < F.rows(); q++)
  {
    double a, b, c;
    a = l(q, 0);
    b = l(q, 1);
    c = l(q, 2);
    vector<double> sides = {a, b, c};

    for (int r = 0; r < 3; r++)
    {
      int i = min(F(q, r), F(q, (r+1) % 3));
      int j = max(F(q, r), F(q, (r+1) % 3));
      P index(i, j);

      double weight =  0.5 * cotangent(sides[r], sides[(r+1) % 3],
	                               sides[(r+2) % 3], dArea(q));

      sum_of_weights += weight;
      map<P, double>::iterator it = laplacian_values.find(index);
      if (it == laplacian_values.end())
      {
	laplacian_values[index] = weight;
      } else
      {
	it->second += weight;
      }
    }
  }

  int number_vertices = F.maxCoeff() + 1;
  std::vector<T> triplet_list;
  triplet_list.reserve(2 * laplacian_values.size() + number_vertices);

  for (map<P, double>::iterator it = laplacian_values.begin();
       it != laplacian_values.end(); it++)
  {
    int i = it->first.first;
    int j = it->first.second;
    double value = it->second;
    // two insertions because L is symmetric
    triplet_list.push_back(T(i, j, value));
    triplet_list.push_back(T(j, i, value));
  }

  for (int i = 0; i < number_vertices; i++)
  {
    triplet_list.push_back(T(i, i, -1 * sum_of_weights));
  }
  L.resize(number_vertices, number_vertices);
  L.setFromTriplets(triplet_list.begin(), triplet_list.end());
}

