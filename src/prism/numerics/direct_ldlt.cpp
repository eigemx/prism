#include "direct_ldlt.h"

namespace prism::solver {

void DirectLDLT::prepare(const SparseMatrix& A) {
    if (!_analyzed || A.rows() != _n || A.nonZeros() != _nnz) {
        _ldlt.analyzePattern(A);
        _analyzed = true;
        _n = A.rows();
        _nnz = A.nonZeros();
    }
    _ldlt.factorize(A);
}

auto DirectLDLT::step(const VectorXd& x, const VectorXd& b) -> VectorXd {
    (void) x;
    return _ldlt.solve(b);
}

} // namespace prism::solver
