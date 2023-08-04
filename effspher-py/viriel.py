"""
Calcul de la surdensité viriel.

1. Calculer R tel que 2Ec + Epot = 0 ou R = Rmax/2
2. Calculer le temps auquel ce rayon est atteint
3. Calculer à quel delta cela correspond
4. Faire évoluer delta selon delta = delta x (a(teff) / a(tvir)) ** 3
"""

import numpy as np

from cosmofunc import H, G, eq_diff, rho_m, rho, a
from solvers import integrateur


class Surdensite:
    def __init__(self, di, ti, equation, tf, dt=1e-5, max_density=4 * 1e7, masse=1e16):
        self.di = di
        self.ti = ti
        self.equation = equation
        self.tf = tf
        self.dt = dt
        self.max_density = max_density
        solution = integrateur(di, ti, eq_diff, tf, dt=dt, max_density=max_density, array=True)
        self. delta, self.ddelta, self.ttab = solution
        self.masse = masse
        self.R = ((3 * masse) / (4 * np.pi * rho_m(self.ttab) * (1 + self.delta))) ** (1 / 3)
        self.tvir = None
        self.rvir = None

    def vitesse(self, r):
        """ Vitesse d'évolution du rayon de la sphère, v = dR/dt"""
        return H(self.ttab) * r - (r * self.ddelta) / (3 * (1 + self.delta))

    def energie_cin(self, r):
        v = self.vitesse(r)
        return (3 * self.masse * v**2) / 10

    def energie_pot(self, r):
        return - (3 * G * self.masse**2) / (5 * r)

    def find_rta(self, energie=False):
        """ Fonction pour calculer le rayon de volte-face (ou donc rayon max)
            et le temps associé.
            Si energie=True, trouve Rta via Ecin = 0
        """
        if not energie:
            imax = self.R.argmax()
            return self.R[imax], self.ttab[imax]

        # trouver Ecin = 0 <=> v = 0
        imax = np.abs(self.vitesse(self.R)).argmin()
        return self.R[imax], self.ttab[imax]

    def find_rvir(self, energie=True):
        """ Fonction pour calculer le rayon viriel défini par Rta/2
            et le temps associé.
            Si energie=True, Trouve R via 2Ecin + Epot = 0, soit v**2 = GM/r
        """
        if not energie:
            raise NotImplementedError("Actuellement impossible de trouver tvir sans la méthode énergie")
            # rta = find_rta(r, t, energie=False)[0]
            # ivir = np.abs(r - rta/2).argmin()
            # return r[ivir], t[ivir]

        vvir = -np.sqrt(G * self.masse / self.R)  # on veut la racine négative puisque v < 0 lorsque rayon se rétracte
        ivir = np.abs(self.vitesse(self.R) - vvir).argmin()
        return self.R[ivir], self.ttab[ivir]

    def surd_vir(self, _rho=False):
        self.rvir, self.tvir = self.find_rvir()
        if _rho:
            return rho(self.rvir, self.masse) / rho_m(self.tvir) - 1
        sol = integrateur(self.di, self.ti, self.equation, t_max=self.tvir, dt=self.dt, max_density=self.max_density)
        return sol[0]

    def surd_finale(self):
        dvir = self.surd_vir()
        if self.tvir:
            return 1 + dvir * (a(self.ttab[-1]) / a(self.tvir)) ** 3
        else:
            return None
