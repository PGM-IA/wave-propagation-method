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
    matplotlib.rcParams.update({'font.size': 16, 'font.family': 'sans-serif'});
    matplotlib.rcParams.update({'text.usetex': 'true'});
    import pylab
    import sys;
    sys.path.append('../lib/')
    import shllw_exact_rs as exact_rs

    # time string
    time_string = '%4.2f' % time;

    # frame string
    nframe_str = f'{nframe:03}';

    Nx    = ini_params['Nx' ];
    xlow  = ini_params['xlow'];
    xhigh = ini_params['xhigh'];
    Neqn  = ini_params['Neqn'];
    dx  = (xhigh-xlow)/Nx;

    # grab left and right states and solve exact Riemann problem
    hL = qsoln[0,0];
    hR = qsoln[Nx-1,0];
    uL = qsoln[0,1]/hL;
    uR = qsoln[Nx-1,1]/hR;
    grav = 1.0;
    h_ex, u_ex = exact_rs.exact_riemann_solution(0.0,hL,uL,hR,uR,xc,time,grav);

    # plot 1
    plt.figure(1);
    plt.clf();
    plt.gca().set_aspect('auto');
    plt.gca().set_xlim([xc[0]-0.5*dx,xc[Nx-1]+0.5*dx]);
    plt.plot(xc,qsoln[:,0],'bo',linewidth=2.0,label='WPM');
    plt.plot(xc,h_ex,'r-',linewidth=2.0,label='exact');
    plt.grid(True);
    title = r"Height: $h(t,x)$ at $t = " + time_string + r"$";
    plt.title(title);
    plt.xlabel(r'$x$');
    plt.ylabel(r'$h$');
    plt.gca().set_xlim([xlow,xhigh])
    plt.legend();

    if (show_all_frames_at_once==False):
        plt.draw();
    else:
        plt.savefig('height_'+nframe_str\
                   +'.pdf', format='pdf', bbox_inches='tight');

    # plot 2
    plt.figure(2);
    plt.clf();
    plt.gca().set_aspect('auto');
    plt.gca().set_xlim([xc[0]-0.5*dx,xc[Nx-1]+0.5*dx]);
    plt.plot(xc,qsoln[:,1]/qsoln[:,0],'bo',linewidth=2.0,label='WPM');
    plt.plot(xc,u_ex,'r-',linewidth=2.0,label='exact')
    plt.grid(True);
    title = r"Velocity: $u(t,x)$ at $t = " + time_string + r"$";
    plt.title(title);
    plt.xlabel(r'$x$');
    plt.ylabel(r'$u$');
    plt.gca().set_xlim([xlow,xhigh])
    plt.legend();

    if (show_all_frames_at_once==False):
        plt.draw();
    else:
        plt.savefig('velocity_'+nframe_str\
                   +'.pdf', format='pdf', bbox_inches='tight');
#----------------------------------------------------------