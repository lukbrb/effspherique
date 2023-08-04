"""
Calcul de la surdensité viriel.

1. Calculer R tel que 2Ec + Epot = 0 ou R = Rmax/2
2. Calculer le temps auquel ce rayon est atteint
3. Calculer à quel delta cela correspond
4. Faire évoluer delta selon delta = delta x (a(teff) / a(tvir)) ** 3
"""

import numpy as np
import matplotlib.pyplot as plt

from cosmofunc import H, G, surd_mini, ti, eq_diff, tf, milliard_annee, rho_m, rho, a
from solvers import integrateur

PLOT = True

tf /= milliard_annee
ti /= milliard_annee
dt = 1e-5
masse = 1e16
delta, ddelta, ttab = integrateur(4 * surd_mini, ti, eq_diff, tf, dt=dt, max_density=4 * 1e7, array=True)

R = ((3 * masse) / (4 * np.pi * rho_m(ttab) * (1 + delta))) ** (1 / 3)


def vitesse(r):
    """ Vitesse d'évolution du rayon de la sphère, v = dR/dt"""
    return H(ttab) * r - (r * ddelta) / (3 * (1 + delta))


def energie_cin(r):
    v = vitesse(r)
    return (3 * masse * v**2) / 10


def energie_pot(r):
    return - (3 * G * masse**2) / (5 * r)


def find_rta(r, t, energie=False):
    """ Fonction pour calculer le rayon de volte-face (ou donc rayon max)
        et le temps associé.
        Si energie=True, trouve Rta via Ecin = 0
    """
    if not energie:
        imax = r.argmax()
        return r[imax], t[imax]

    # trouver Ecin = 0 <=> v = 0
    imax = np.abs(vitesse(r)).argmin()
    return r[imax], t[imax]


def find_rvir(r, t, energie=False):
    """ Fonction pour calculer le rayon viriel défini par Rta/2
        et le temps associé.
        Si energie=True, Trouve R via 2Ecin + Epot = 0, soit v**2 = GM/r
    """
    if not energie:
        raise NotImplementedError("Actuellement impossible de trouver tvir sans la méthode énergie")
        # rta = find_rta(r, t, energie=False)[0]
        # ivir = np.abs(r - rta/2).argmin()
        # return r[ivir], t[ivir]

    vvir = -np.sqrt(G * masse / r)  # on veut la racine négative puisque v < 0 lorsque rayon se rétracte
    ivir = np.abs(vitesse(r) - vvir).argmin()
    return r[ivir], t[ivir]


Rta, tmax = find_rta(R, ttab, energie=True)
Rvir, tvir = find_rvir(R, ttab, energie=True)
delta_vir, _, _ = integrateur(4 * surd_mini, ti, eq_diff, t_max=tvir, dt=dt, max_density=4 * 1e7)

delta_vir2 = rho(Rvir, masse) / rho_m(tvir) - 1
print(delta_vir, delta_vir2, (18 * np.pi**2) / ((a(ttab[-1]) / a(tvir)) ** 3))
delta_fin = 1 + delta_vir * (a(ttab[-1]) / a(tvir)) ** 3
delta_fin2 = 1 + delta_vir2 * (a(ttab[-1]) / a(tvir)) ** 3

print(delta_fin, delta_fin2, 18 * np.pi**2)

if PLOT:
    plt.figure()
    plt.title("Évolution du rayon")
    plt.plot(ttab, R, 'r', label="Rayon")
    plt.plot(tmax, Rta, 'x', label="Rayon max")
    plt.plot(tvir, Rvir, 'x', label="Rayon viriel")
    plt.xlabel('Temps (Gyr)')
    plt.ylabel('R(t)')
    # plt.yscale('log')
    plt.legend()

    drdt = vitesse(R)
    plt.figure()
    plt.title("Profil des vitesses")
    plt.plot(ttab, drdt, '--k', label="Vitesse")
    plt.vlines(tmax, drdt.min(), drdt.max(), colors='red', label="Vitesse au rayon max")
    plt.vlines(tvir, drdt.min(), drdt.max(), colors='orange', label="Vitesse au rayon viriel")
    plt.xlabel('Temps (Gyr)')
    plt.ylabel('R(t)')
    # plt.yscale('log')
    plt.legend()

    plt.figure()
    plt.title("Évolution de la surdensité")
    plt.plot(ttab, delta, '--k')
    plt.xlabel('Temps (Gyr)')
    plt.ylabel('Surdensité $\delta$')
    plt.yscale('log')
    plt.show()
