#include "waveprop1D.h"
#include "DogArrays.h"
#include "constants.h"
#include <cmath>

// -------------------------------------------------------------------------- //
void waveprop1D::set_initial(const DblArray& x,
    DblArray& qsoln)
{
    const int Nx = x.get_ind_length(1);
    const double gamma_gas = 1.4;

    for (int i = 1; i <= Nx; i++)
    {
        const double xi = x.get(i);

        double rho, u, p;
        if (xi < 0.0)
        {
            // Left state
            rho = 1.0;
            u = 0.0;
            p = 1.0;
        }
        else
        {
            // Right state
            rho = 0.125;
            u = 0.0;
            p = 0.1;
        }

        const double kinetic = 0.5 * rho * u * u;
        const double E = p / (gamma_gas - 1.0) + kinetic;

        // Conserved variables: [rho, rho*u, E]
        qsoln.fetch(i, 1) = rho;
        qsoln.fetch(i, 2) = rho * u;
        qsoln.fetch(i, 3) = E;
    }
}
//--------------------------------------------------------------------- //