#include "waveprop1D.h"
#include "DogArrays.h"
#include "constants.h"
#include <cmath>

// -------------------------------------------------------------------------- //
void waveprop1D::set_initial(const DblArray& x,
                             DblArray& qsoln)
{
  const int     Nx = x.get_ind_length(1);
  const double  hl = 3.0, ul = 0.0;
  const double  hr = 1.0, ur = 0.0;

  for (int i=1; i<=Nx; i++)
    { 
      if (x.get(i)<0.0)
        {
          qsoln.fetch(i,1) = hl;
          qsoln.fetch(i,2) = hl*ul;
        }
      else
        {
          qsoln.fetch(i,1) = hr;
          qsoln.fetch(i,2) = hr*ur;
        }
    }
}
// -------------------------------------------------------------------------- //