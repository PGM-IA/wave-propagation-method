#include "AppSolver.h"
#include "DogArrays.h"
#include "Params1D.h"
#include "dog_math.h"
#include <cmath>
#include <stdio.h>

// gravitational constant
const double gamma_gas = 1.4;

// -------------------------------------------------------------------------- //
double AppSolver::Solver(const Params1D& params1D,
                         const DblArray& QL,
                         const DblArray& QR,
                         DblArray& s,
                         DblArray& wave,
                         DblArray& amdq,
                         DblArray& apdq)
{
  // Assume Euler 3 conserved variables
  // Q = [rho, rho*u, E]
  // Gas constant is defined above
  
  // Left & right conserved primitive variables
  const double rhoL = QL.get(1);
  const double mL   = QL.get(2); 
  const double EL   = QL.get(3); 
   
  const double rhoR = QR.get(1);
  const double mR   = QR.get(2);
  const double ER   = QR.get(3);

  const double uL = mL/rhoL;
  const double uR = mR/rhoR;

  const double pL = (gamma_gas - 1.0) * (EL - 0.5 * rhoL * uL * uL);
  const double pR = (gamma_gas - 1.0) * (ER - 0.5 * rhoR * uR * uR);

  const double HL = (EL + pL) / rhoL;  // total enthalpy left
  const double HR = (ER + pR) / rhoR;  // total enthalpy right

  // Roe-averaged states
  const double sqL   = sqrt(rhoL);
  const double sqR   = sqrt(rhoR);
  const double denom = sqL + sqR;

  const double uHat = (sqL * uL + sqR * uR) / denom;
  const double HHat = (sqL * HL + sqR * HR) / denom;

  const double cHat = sqrt((gamma_gas - 1.0) * (HHat - 0.5 * uHat * uHat));


  // wave speeds
  s.fetch(1) = uHat - cHat;
  s.fetch(2) = uHat;
  s.fetch(3) = uHat + cHat;

  double smax = std::fabs(s.get(1));
  smax = dog_math::Max(smax, std::fabs(s.get(2)));
  smax = dog_math::Max(smax, std::fabs(s.get(3)));

  // Right eigenvectors
  // r1 = [1, u-c, H-u*c]^T
  // r2 = [1,   u, 0.5*u^2]^T
  // r3 = [1, u+c, H+u*c]^T
  double R[3][3];

  R[0][0] = 1.0;
  R[1][0] = uHat - cHat;
  R[2][0] = HHat - uHat * cHat;

  R[0][1] = 1.0;
  R[1][1] = uHat;
  R[2][1] = 0.5 * uHat * uHat;

  R[0][2] = 1.0;
  R[1][2] = uHat + cHat;
  R[2][2] = HHat + uHat * cHat;

  // Decompose jump deltaQ = Q_R - Q_L into waves
  double dQ[3];
  dQ[0] = rhoR - rhoL;
  dQ[1] = mR - mL;
  dQ[2] = ER - EL;

  // Solve R * alpha = dQ via small 3x3 Gaussian elimination
  double A[3][4];
  for (int i = 0; i < 3; ++i)
  {
      for (int j = 0; j < 3; ++j)
          A[i][j] = R[i][j];
      A[i][3] = dQ[i];
  }

  for (int i = 0; i < 3; ++i)
  {
      double piv = A[i][i];
      if (std::fabs(piv) < 1e-14)
          piv = (piv >= 0.0 ? 1e-14 : -1e-14);

      for (int j = i; j < 4; ++j)
          A[i][j] /= piv;

      for (int k = 0; k < 3; ++k)
      {
          if (k == i) continue;
          double fac = A[k][i];
          for (int j = i; j < 4; ++j)
              A[k][j] -= fac * A[i][j];
      }
  }

  double alpha[3];
  for (int i = 0; i < 3; ++i)
      alpha[i] = A[i][3];

  // Build waves and fluctuations
  wave.setall(0.0);
  amdq.setall(0.0);
  apdq.setall(0.0);

  for (int p = 0; p < 3; ++p)       // p = wave index: 0,1,2
  {
      const double sp = s.get(p + 1);   // corresponding speed

      for (int m = 0; m < 3; ++m)    // m = equation index: 0,1,2
      {
          const double Wmp = alpha[p] * R[m][p];  // component of wave p in eqn m

          // In this code base, wave(m+1, p+1) = component m of wave p
          wave.fetch(m + 1, p + 1) = Wmp;

          if (sp < 0.0)
          {
              amdq.fetch(m + 1) += sp * Wmp;
          }
          else if (sp > 0.0)
          {
              apdq.fetch(m + 1) += sp * Wmp;
          }
      }
  }

  return smax;
}