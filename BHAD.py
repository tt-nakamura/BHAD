import numpy as np
from scipy.optimize import newton
from constants import G,c,Msun,MEsun,a_rad,stefan,kB,mH,me,pc
from radiation import kappa_es,kappa_0,kappa_ff,nu_0,nu_E
from radiation import BlackBody,InvCompton

class StandardModel:
    """ Standard model of accretion disks
    reference: N.I.Shakura and R.A.Sunyaev,
    Astronomy and Astrophysics 24 (1973) 337
    """
    def __init__(self, M, M_dot, alpha,
                 r_out, r_in=3, N=256):
        """
        M = mass of black hole / solar mass
        M_dot = mass accretion rate / critical rate
        alpha = viscosity parameter (0<alpha<=1)
        r_out = radius of disk's outer edge / 2GM/c^2
        r_in = radius of disk's inner edge / 2GM/c^2
        N = number of grid points in radial coordinate
            (excluding inner boundary)
        """
        rg = 2*G*M*Msun/c**2
        self.M = M*Msun
        self.M_dot = M_dot*M*MEsun
        self.alpha = alpha
        self.r_in = r_in*rg
        self.r_out = r_out*rg
        self.r = np.geomspace(self.r_in, self.r_out, N+1)[1:]
        self.y = StandardModel.solve(self, self.r)

    def solve(self, r, Sigma=None): 
        """
        r = distance from disk center (black hole) / cm
          r can be an array (assume increasing order)
        Sigma = initial guess for Surface density at r / g/cm^2;
          if r is array, Sigma is used only at r[-1];
          if Sigma is None, asymptotic solution at large r
          is used as an initial guess
        return y = (Sigma,p,T,H,F,tau) where
          Sigma = surface density at r / g/cm^2
          p = pressure at r / dyne / cm^2
          T = temperature on disk at r / K
          H = disk height at r / cm
          F = flux from one side of disk / erg/cm^2
          tau = optical depth at r
        if r is scalar, y is tuple of length 6
        else if r has shape (n,), y has shape (6,n)
        """
        if hasattr(r, '__len__'):
            Y,y = [],[Sigma]
            for r in r[::-1]:
                y = StandardModel.solve(self, r, y[0])
                Y.append(y)
            return np.array(Y[::-1]).T

        if r <= self.r_in: return (0,)*6
        Omega = np.sqrt(G*self.M/r**3) # Keplerian angular velocity
        f = self.M_dot*(1 - (self.r_in/r)**0.5)*Omega/2/np.pi
        kBmH = 2*kB/mH
        Sigma0 = (256*stefan*f**7/Omega**2# asymptotic solution
                  /(9*self.alpha**8*kappa_0*kBmH**7.5))**0.1

        def equation(Sigma):
            nonlocal p,T,H,F,tau
            if Sigma <= 0: Sigma = Sigma0
            cs2 = f/Sigma/self.alpha # sound spped squared
            H = cs2**0.5/Omega # disk height
            rho = Sigma/2/H # density
            p = cs2*rho # pressure
            if T is None: T = cs2/kBmH # temperature
            T = newton(lambda T: rho*kBmH*T + a_rad/3*T**4 - p, T,
                       lambda T: rho*kBmH + 4*a_rad/3*T**3)
            kappa = kappa_es + kappa_ff(rho,T) # opacity
            tau = kappa*Sigma/2 # optical depth
            Q = 1.5*f*Omega # heat generation
            F = 16/3*stefan*T**4/tau # flux from one side
            return Q - 2*F # equation for radiative equilibrium

        if Sigma is None: Sigma = Sigma0
        p,T,H,F,tau = [None]*5
        Sigma = newton(equation, Sigma)
        return (Sigma,p,T,H,F,tau)

    def radius(self): return self.r
    def surface_density(self): return self.y[0]
    def pressure(self): return self.y[1]
    def temperature(self): return self.y[2]
    def height(self): return self.y[3]
    def flux(self): return self.y[4]
    def optical_depth(self): return self.y[5]
    def infall_velocity(self): return self.M_dot/2/np.pi/self.r/self.y[0]
    def rotation_velocity(self): return (G*self.M/self.r)**0.5
    def photon_pressure(self): return a_rad*self.y[2]**4/3
    def gas_pressure(self): return 2*kB/mH*self.density()*self.y[2]
    def density(self): return self.y[0]/2/self.y[3]
    def effective_temperature(self): return (self.y[4]/stefan)**0.25
    def opacity(self):
        return kappa_es + kappa_ff(self.density(), self.y[2])

    def spectrum(self, nu, D=10*pc, i=0):
        """ spectrum of emitted photons from disk surface
        nu = photon frequency / Hz
        D = distance to obersever / cm
        i = inclination angle wrt line of sight / radian
        return d(intensity)/dln(nu) / erg/s/cm^2
        """
        Teff = self.effective_temperature()
        I = BlackBody(np.expand_dims(nu, -1), Teff)
        I = np.trapz(I*self.r, self.r)
        return 2*np.pi*np.cos(i)/D**2*I


class SLEModel(StandardModel):
    """ Optically thin model of accretion disks
    reference: S.L.Shapiro, A.P.Lightman and D.M.Eardley
    The Astrophysical Journal 204 (1976) 187
    """
    def __init__(self, M, M_dot, alpha,
                 r_out, r_in=3, N=256):
        """ override StandardModel.__init__ """
        super().__init__(M, M_dot, alpha, r_out, r_in, N)
        self.y = np.vstack([self.y, self.temperature()])
        p_rad = self.photon_pressure()
        p_gas = self.gas_pressure()
        hot = (p_rad > 3*p_gas)
        self.is_hot = np.any(hot)
        if self.is_hot:
            self.hot = len(hot) - np.argmax(hot[::-1])
            y = SLEModel.solve(self, self.r[:self.hot])
            self.y[:,:self.hot] = y

    def solve(self, r, Sigma=None):
        """ override StandardModel.solve
        return y = (Sigma,p,Ti,H,F,tau,Te) where
          Ti = ion temperature on disk at r / K
          Te = electron temperature on disk at r / K
        if r is scalar, y is tuple of length 7
        else if r has shape (n,), y has shape (7,n)
        """
        if hasattr(r, '__len__'):
            Y,y = [],[Sigma]
            for r in r[::-1]:
                y = SLEModel.solve(self, r, y[0])
                Y.append(y)
            return np.array(Y[::-1]).T

        if r <= self.r_in: return (0,)*7
        Omega = np.sqrt(G*self.M/r**3) # Keplerian angular velocity
        f = self.M_dot*(1 - (self.r_in/r)**0.5)*Omega/2/np.pi
        mckB = me*c**2/kB
        mcke = me*c**2/kappa_es
        Sigma0 = ((f*mcke/kB)**3/self.alpha # solution for tau<1
                  /2/nu_0**2/(f/self.alpha - mcke/mH)**2)**(1/6)

        def equation(Sigma):
            nonlocal p,Te,Ti,H,F,tau
            if Sigma <= 0: Sigma = Sigma0
            cs2 = f/Sigma/self.alpha # sound spped squared
            H = cs2**0.5/Omega # disk height
            rho = Sigma/2/H # density
            p = cs2*rho # pressure
            tau = kappa_es*Sigma/2 # optical depth
            Te = mckB/4/tau/max(1,tau) # electron temperature
            Ti = cs2*mH/kB - Te # ion temperature
            F = 0.75*nu_E(rho,Te)*Sigma*kB*(Ti-Te)/mH # flux from one side
            Q = 1.5*f*Omega # heat generation
            return Q - 2*F # equation for radiative equilibrium

        if Sigma is None: Sigma = Sigma0
        p,Te,Ti,H,F,tau = [None]*6
        Sigma = newton(equation, Sigma)
        return (Sigma,p,Ti,H,F,tau,Te)

    def ion_temperature(self): return self.y[2]
    def electron_temperature(self): return self.y[6]
    def gas_pressure(self):
        return kB/mH*self.density()*(self.y[2] + self.y[6])

    def spectrum(self, nu, Es, D=10*pc, i=0):
        """ override StandardModel.spectrum
        Es = energy of photons supplied / erg
        """
        I = np.empty((len(nu), len(self.r)))
        nu = np.expand_dims(nu, -1)
        Teff = (self.y[4,self.hot:]/stefan)**0.25
        I[:,self.hot:] = BlackBody(nu, Teff)
        if self.is_hot:
            Te = self.y[6,:self.hot]
            I0 = self.y[4,:self.hot]/np.pi
            I[:,:self.hot] = I0 * InvCompton(nu, Te, Es)

        I = np.trapz(I*self.r, self.r)
        return 2*np.pi*np.cos(i)/D**2*I
