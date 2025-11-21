def exact_riemann_solution(shock_location,h_l,u_l,h_r,u_r,x=None,t=None,grav=1.0):
    """
    Modified from the code written by David Ketcheson (http://www.davidketcheson.info/) as downloaded from
    https://gist.github.com/ketch/08ce0845da0c8f3fa9ff.

    This code is made available under the Creative Commons Attribution 4.0 International
    License (https://creativecommons.org/licenses/by/4.0/).

    Return the exact solution to the Riemann problem with initial states (hl,ul) and (hr,ur).
       The solution is computed at time t and points x (where x may be a 1D numpy array).
    """

    import numpy as np;
    from scipy.optimize import fsolve;

    # Compute left and right state gravity wave speeds
    c_l = np.sqrt(grav*h_l)
    c_r = np.sqrt(grav*h_r)
    R1 = u_l+2.0*c_l;
    R2 = u_r-2.0*c_r;

    # Check for cavitation
    if (R1-R2)**2 < 0.0:
        print('Cavitation detected!  Exiting.');
        return None;

    # Define the integral curves and hugoniot loci
    integral_curve_1 = lambda h : u_l + 2.0*(c_l-np.sqrt(grav*h));
    hugoniot_locus_1 = lambda h : u_l - (h-h_l)*np.sqrt(0.5*grav*(1.0/h+1.0/h_l));
    integral_curve_2 = lambda h : u_r - 2.0*(c_r-np.sqrt(grav*h));
    hugoniot_locus_2 = lambda h : u_r + (h-h_r)*np.sqrt(0.5*grav*(1.0/h+1.0/h_r));

    # Check whether the 1-wave is a shock or rarefaction
    def phi_l(h):
        if h<0.0: return np.nan;
        elif h>=h_l: return hugoniot_locus_1(h)
        else: return integral_curve_1(h)

    # Check whether the 1-wave is a shock or rarefaction
    def phi_r(h):
        if h<0.0: return np.nan;
        elif h>=h_r: return hugoniot_locus_2(h)
        else: return integral_curve_2(h)

    phi = lambda h : phi_l(h)-phi_r(h)

    # Compute middle state p, u by finding curve intersection
    h,info, ier, msg = fsolve(phi, 0.5*(h_l+h_r),full_output=True,xtol=1.e-14)
    # For strong rarefactions, sometimes fsolve needs help
    if ier!=1:
        h,info, ier, msg = fsolve(phi, 0.5*(h_l+h_r),full_output=True,factor=0.1,xtol=1.e-10)
        # This should not happen:
        if ier!=1:
            print('Warning: fsolve did not converge.');
            print(msg);

    u = phi_l(h);

    # compute the wave speeds
    ws = np.zeros(4)

    # Find shock and rarefaction speeds
    if h>h_l:
        ws[0] = (h_l*u_l - h*u)/(h_l - h)
        ws[1] = ws[0]
    else:
        c_l_star = np.sqrt(grav*h)
        ws[0] = u_l - c_l
        ws[1] = u - c_l_star

    if h>h_r:
        ws[3] = (h_r*u_r - h*u)/(h_r - h)
        ws[2] = ws[3]
    else:
        c_r_star = np.sqrt(grav*h)
        ws[2] = u + c_r_star
        ws[3] = u_r + c_r


    # Compute return values

    # Choose a time based on the wave speeds
    if x is None: x = np.linspace(shock_location-1.,shock_loation+1.,1000)
    if t is None: t = 0.8*max(np.abs(x-shock_location))/max(np.abs(ws))

    xs = shock_location + ws*t # Wave locations

    # Find solution inside rarefaction fans
    if (t>0.0):
        xi = (x-shock_location)/t

        h1 = (1.0/(9.0*grav))*(R1-xi)**2;
        h2 = (1.0/(9.0*grav))*(R2-xi)**2;
        u1 = R1 - 2.0*np.sqrt(grav*h1);
        u2 = R2 + 2.0*np.sqrt(grav*h2);

        h_out = (x<=xs[0])*h_l + (x>xs[0])*(x<=xs[1])*h1 + (x>xs[1])*(x<=xs[2])*h + (x>xs[2])*(x<=xs[3])*h2 + (x>xs[3])*h_r
        u_out = (x<=xs[0])*u_l + (x>xs[0])*(x<=xs[1])*u1 + (x>xs[1])*(x<=xs[2])*u + (x>xs[2])*(x<=xs[3])*u2 + (x>xs[3])*u_r
    else:
        h_out = (x<=shock_location)*h_l + (x>shock_location)*h_r;
        u_out = (x<=shock_location)*u_l + (x>shock_location)*u_r;

    return h_out, u_out

if __name__ == "__main__":

    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    matplotlib.rcParams.update({'font.size': 16, 'font.family': 'sans-serif'});
    matplotlib.rcParams.update({'text.usetex': 'true'});

    # exact Riemann solution
    shock_location = 0.5; grav = 1.0;
    #
    hL = 3.0;  uL = 0.0;
    hR = 1.0;  uR = 0.0;
    #
    time = 0.2;
    xc = np.linspace(0.0,1.0,5001);

    h_ex, u_ex = exact_riemann_solution(shock_location,hL,uL,hR,uR,xc,time,grav);

    # time string
    time_string = '%4.2f' % time;

    # height
    plt.figure(1)
    plt.clf()
    plt.grid()
    xlow  =  0.0;
    xhigh =  1.0;
    xticks = np.linspace(xlow,xhigh,11);
    plt.plot(xc,h_ex,'r-',linewidth=3.0);
    plt.grid(True);
    title = "Height: " + r'$h(t,x)$' + " at $t = " + time_string + "$";
    plt.title(title);
    plt.xlabel(r'$x$')
    plt.ylabel(r'height: $h$')
    plt.xlim(xlow,xhigh)
    plt.xticks(xticks)
    plt.savefig("shllw_dambreak_height.pdf", format='pdf')

    # velocity
    plt.figure(2)
    plt.clf()
    plt.grid()
    plt.plot(xc,u_ex,'r-',linewidth=3.0);
    plt.grid(True);
    title = "Velocity: " + r'$u(t,x)$' + " at $t = " + time_string + "$";
    plt.title(title);
    plt.xlabel(r'$x$')
    plt.ylabel(r'velocity: $u$')
    plt.xlim(xlow,xhigh)
    plt.xticks(xticks)
    plt.savefig("shllw_dambreak_velocity.pdf", format='pdf')

    #plt.show()