#include "smooth.h"
#include <iostream>
#include <Eigen/SparseCholesky>
#include "igl/edge_lengths.h"
#include "cotmatrix.h"
#include "massmatrix.h"

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  using Sparse = Eigen::SparseMatrix<double>;
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);

  Sparse L;
  cotmatrix(l, F, L);
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
  massmatrix(l, F, M);

  Sparse A;
  Sparse M_sparse(L.rows(), L.cols());

  A = Sparse(M) - lambda * L;
  Eigen::SimplicialLDLT<Sparse> solver(A);
  U = solver.solve(M * G);
}
