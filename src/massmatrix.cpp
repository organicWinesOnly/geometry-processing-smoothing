#include "massmatrix.h"
#include "igl/doublearea.h"
#include <map>

typedef Eigen::Triplet<double> T;
typedef std::map<int, double> Map;
void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  Eigen::VectorXd Area;
  igl::doublearea(l, Area);
  Area = Area / 2.0;

  int number_of_vertices = F.maxCoeff() + 1;
  M.resize(number_of_vertices);

  for (int i = 0; i < F.rows(); i++)
  {
    for (int j = 0; j < 3; j++)
    {
      M.diagonal()[F(i, j)] += 1.0 / 3.0 * Area(i);
    }
  }
}

