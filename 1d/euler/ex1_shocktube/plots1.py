#----------------------------------------------------------
def plots1(output_directory,        # directory where output data lives
           show_all_frames_at_once, # True ==> "show all frames at once"
           nframe,                  # current frame
           ini_params,              # collection of all parameters read from "parameters.ini"
           time,                    # current time
           xc,                      # spatial mesh (cell centers)
           xl,                      # spatial mesh (cell corners)
           qsoln):                  # q solution

    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    matplotlib.rcParams.update({'font.size': 16, 'font.family': 'sans-serif'})
    matplotlib.rcParams.update({'text.usetex': 'true'})

    # ---------------- basic info ----------------
    time_string = '%4.2f' % time
    nframe_str  = f'{nframe:03}'

    Nx    = ini_params['Nx']
    xlow  = ini_params['xlow']
    xhigh = ini_params['xhigh']
    Neqn  = ini_params['Neqn']
    dx    = (xhigh - xlow) / Nx

    assert Neqn == 3, "Euler solver should have Neqn = 3 (rho, rho*u, E)"

    gamma_gas = 1.4

    # qsoln shape: (Nx, Neqn)
    rho = qsoln[:,0]
    m   = qsoln[:,1]
    E   = qsoln[:,2]

    u = m / rho
    p = (gamma_gas - 1.0) * (E - 0.5 * rho * u * u)

    # ---------------- density plot ----------------
    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0] - 0.5*dx, xc[Nx-1] + 0.5*dx])
    plt.plot(xc, rho, 'bo', linewidth=2.0, label='WPM')
    plt.grid(True)
    title = r"Density: $\rho(t,x)$ at $t = " + time_string + r"$"
    plt.title(title)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\rho$')
    plt.gca().set_xlim([xlow, xhigh])
    plt.legend()

    if not show_all_frames_at_once:
        plt.draw()
    else:
        plt.savefig('density_' + nframe_str + '.pdf',
                    format='pdf', bbox_inches='tight')

    # ---------------- velocity plot ----------------
    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0] - 0.5*dx, xc[Nx-1] + 0.5*dx])
    plt.plot(xc, u, 'bo', linewidth=2.0, label='WPM')
    plt.grid(True)
    title = r"Velocity: $u(t,x)$ at $t = " + time_string + r"$"
    plt.title(title)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u$')
    plt.gca().set_xlim([xlow, xhigh])
    plt.legend()

    if not show_all_frames_at_once:
        plt.draw()
    else:
        plt.savefig('velocity_' + nframe_str + '.pdf',
                    format='pdf', bbox_inches='tight')

    # ---------------- pressure plot ----------------
    plt.figure(3)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0] - 0.5*dx, xc[Nx-1] + 0.5*dx])
    plt.plot(xc, p, 'bo', linewidth=2.0, label='WPM')
    plt.grid(True)
    title = r"Pressure: $p(t,x)$ at $t = " + time_string + r"$"
    plt.title(title)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$p$')
    plt.gca().set_xlim([xlow, xhigh])
    plt.legend()

    if not show_all_frames_at_once:
        plt.draw()
    else:
        plt.savefig('pressure_' + nframe_str + '.pdf',
                    format='pdf', bbox_inches='tight')
#----------------------------------------------------------