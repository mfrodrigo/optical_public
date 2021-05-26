def ssprop(u0, dt, dz, nz, alpha, betap, gamma=0, maxiter=4, tol=1e-5):

    # This function solves the nonlinear Schrodinger equation for
    # pulse propagation in an optical fiber using the split-step
    # Fourier method.
    #
    # The following effects are included in the model: group velocity
    # dispersion (GVD), higher order dispersion, loss, and self-phase
    # modulation (gamma).
    #
    # USAGE
    #
    # u1 = ssprop(u0,dt,dz,nz,alpha,betap,gamma);
    # u1 = ssprop(u0,dt,dz,nz,alpha,betap,gamma,maxiter);
    # u1 = ssprop(u0,dt,dz,nz,alpha,betap,gamma,maxiter,tol);
    #
    # INPUT
    #
    # u0 - starting field amplitude (vector)
    # dt - time step
    # dz - propagation stepsize
    # nz - number of steps to take, ie, ztotal = dz*nz
    # alpha - power loss coefficient, ie, P=P0*exp(-alpha*z)
    # betap - dispersion polynomial coefs, [beta_0 ... beta_m]
    # gamma - nonlinearity coefficient
    # maxiter - max number of iterations (default = 4)
    # tol - convergence tolerance (default = 1e-5)
    #
    # OUTPUT
    #
    # u1 - field at the output
    #
    # NOTES  The dimensions of the input and output quantities can
    # be anything, as long as they are self consistent.  E.g., if
    # |u|^2 has dimensions of Watts and dz has dimensions of
    # meters, then gamma should be specified in W^-1*m^-1.
    # Similarly, if dt is given in picoseconds, and dz is given in
    # meters, then beta(n) should have dimensions of ps^(n-1)/m.
    #
    # See also:  sspropc (compiled MEX routine)
    #
    # AUTHOR: Rodrigo Machado Fonseca (rodrigomf43@gmail.com)

    # Libraries
    import math
    import numpy as np
    from math import e
    import warnings

    nt = u0.shape[0]
    T = nt*dt

    w = list(map(lambda x: 2*math.pi*x/T, np.arange(0, nt/2))) + \
        list(map(lambda x: 2*math.pi*x/T, np.arange(-nt/2, 0)))

    len_w = len(w)
    w = np.array(w).reshape(len_w, 1)

    halfstep = -alpha/2 + np.zeros((len_w, 1))

    for ii in range(0, len(betap)):

        # halfstep = halfstep - ((1j)**(ii-1))*betap[ii, 0]*(w**ii)/math.factorial(ii)
        halfstep = halfstep - 1j * betap[ii, 0] * (w ** ii) / math.factorial(ii)

    halfstep = e**(halfstep*dz/2)

    u1 = u0
    ufft = np.fft.fft(u0, axis=0)

    for iz in range(0, nz):

        uhalf = np.fft.ifft(halfstep*ufft, axis=0)

        for ii in range(0, maxiter):

            # uv = uhalf*e**(1j*gamma*(abs(u1)**2 + abs(u0)**2)*dz/2)
            uv = uhalf * e ** (-1j * gamma * (abs(u1) ** 2 + abs(u0) ** 2) * dz / 2)
            uv = np.fft.fft(uv, axis=0)
            ufft = halfstep*uv
            uv = np.fft.ifft(ufft, axis=0)

            if np.linalg.norm((uv-1))/np.linalg.norm(u1) < tol:

                u1 = uv
                break

            else:

                u1 = uv

        if ii == maxiter-1:

            warnings.warn(
                "Failed to converge to {} in {} iterations".format(tol, maxiter))
        u0 = u1

    return list(u1)
