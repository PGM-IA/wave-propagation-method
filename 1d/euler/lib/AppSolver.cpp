#include "AppSolver.h"
#include "DogArrays.h"
#include "Params1D.h"
#include "dog_math.h"
#include <cmath>
#include <stdio.h>

// gravitational constant
const double grav = 1.0;

// -------------------------------------------------------------------------- //
double AppSolver::Solver(const Params1D& params1D,
                         const DblArray& QL,
                         const DblArray& QR,
                         DblArray& s,
                         DblArray& wave,
                         DblArray& amdq,
                         DblArray& apdq)
{
  // Roe averages
  const double hl = QL.get(1); 
  const double ul = QL.get(2)/hl; 
   
  const double hr = QR.get(1);
  const double ur = QR.get(2)/hr;

  const double hav = 0.5*(hl + hr);
  const double uav = (sqrt(hl)*ul + sqrt(hr)*ur)/(sqrt(hr)+sqrt(hl));

  // speeds
  s.fetch(1) = uav - sqrt(grav*hav);
  s.fetch(2) = uav + sqrt(grav*hav);
  double smax = dog_math::Max(fabs(s.get(1)),fabs(s.get(2)));

  // right eigenvectors
  DblArray rmat(2,2);

  rmat.fetch(1,1) = 1.0;
  rmat.fetch(2,1) = s.get(1);

  rmat.fetch(1,2) = 1.0;
  rmat.fetch(2,2) = s.get(2);

  // left eigenvectors
  DblArray lmat(2,2);

  lmat.fetch(1,1) = s.get(2)/(s.get(2)-s.get(1));
  lmat.fetch(1,2) = -1.0/(s.get(2)-s.get(1));

  lmat.fetch(2,1) = -s.get(1)/(s.get(2)-s.get(1));
  lmat.fetch(2,2) = 1.0/(s.get(2)-s.get(1));

  // wave strengths
  DblArray alpha(2); alpha.setall(0.0);
  for (int m1=1; m1<=2; m1++)
    for (int m2=1; m2<=2; m2++)
      { alpha.fetch(m1) += lmat.get(m1,m2)*(QR.get(m2)-QL.get(m2)); }

  // waves
  wave.setall(0.0);
  for (int m1=1; m1<=2; m1++)
    for (int m2=1; m2<=2; m2++)
      { wave.fetch(m1,m2) = alpha.get(m2)*rmat.get(m1,m2); } 

  // fluctuations
  amdq.setall(0.0);
  apdq.setall(0.0);
  for (int m2=1; m2<=2; m2++)
    {
      if (s.get(m2)<=0.0)
        { 
          for (int m1=1; m1<=2; m1++)
            { amdq.fetch(m1) += s.get(m2)*wave.get(m1,m2); }
        }
      else
        { 
          for (int m1=1; m1<=2; m1++)
            { apdq.fetch(m1) += s.get(m2)*wave.get(m1,m2); }
        }
    }

  // entropy fix
  if (params1D.get_use_entropy_fix())
    {
      DblArray QM(2);
      QM.fetch(1) = QL.get(1) + wave.get(1,1);
      QM.fetch(2) = QL.get(2) + wave.get(2,1);

      // first wave
      {
        double eig1L = QL.get(2)/QL.get(1) - sqrt(grav*QL.get(1));
        double eig1R = QM.get(2)/QM.get(1) - sqrt(grav*QM.get(1));
        double beta1 = (eig1R - s.get(1))/(eig1R - eig1L);
        if (eig1L<0.0 && eig1R>0.0)
          {
            for (int m=1; m<=2; m++)
              { 
                apdq.fetch(m) += ( eig1R*(1.0-beta1) - dog_math::Max(0.0, s.get(1)) ) * wave.get(m,1);
                amdq.fetch(m) += ( eig1L*(beta1) - dog_math::Min(0.0, s.get(1)) ) * wave.get(m,1);
              }
          }
      }

      // second wave
      {
        double eig2L = QM.get(2)/QM.get(1) + sqrt(grav*QM.get(1));;
        double eig2R = QR.get(2)/QR.get(1) + sqrt(grav*QR.get(1));;
        double beta2 = (eig2R - s.get(2))/(eig2R - eig2L);
        if (eig2L<0.0 && eig2R>0.0)
          {
            for (int m=1; m<=2; m++)
              { 
                apdq.fetch(m) += ( eig2R*(1.0-beta2) - dog_math::Max(0.0, s.get(2)) ) * wave.get(m,2);
                amdq.fetch(m) += ( eig2L*(beta2) - dog_math::Min(0.0, s.get(2)) ) * wave.get(m,2);
              }
          }
      }
    }

  // return local max speed
  return smax;
}
// -------------------------------------------------------------------------- //