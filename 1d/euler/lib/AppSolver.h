#ifndef _SOLVER_1D_H_
#define _SOLVER_1D_H_

class Params1D;
class DblArray;

namespace AppSolver
{
  double Solver(const Params1D& params1D,
                const DblArray& QL,
                const DblArray& QR,
                DblArray& s,
                DblArray& wave,
                DblArray& amdq,
                DblArray& apdq);
}
#endif