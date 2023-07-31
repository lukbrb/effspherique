import numpy as np
import matplotlib.pyplot as plt

from cosmofunc import rho_m, G, a, H, milliard_annee, surd_mini, eq_diff, tf, ti
from solvers import rk4

PLOT = True
tf /= milliard_annee
ti /= milliard_annee

init = (4 * surd_mini, 4 * surd_mini * H(ti), ti)
delta, p, ttab = rk4(init, eq_diff, tf, dt=1e-5, max_density=3 * 1e4)

teff = ttab[-1]
Masse = 1e16  # N'influe pas sur le résultat final de la surdensité


# Fonction pour déterminer t lorsque Rviriel est atteint
def t_viriel(borne_min, borne_max, tolerance, r, r_vir):
    for i in r[borne_min:borne_max]:
        inter = abs(i - r_vir)
        t = ttab[borne_min]
        while inter > tolerance:
            milieu = (borne_min + borne_max) / 2
            if i > r_vir:
                borne_min = milieu
            else:
                borne_max = milieu
            inter = borne_max - borne_min
            t += 1
        return t


# Fonction qui fait évoluer la surdensité après la surdensité viriel
def surdensite(_delta, index):
    _delta[index:] = rho_vir / rho_m(ttab[index:]) - 1
    return _delta


def vitesse(r_final):
    r = r_final
    _v = H(ttab) * r * (1 - (p * r ** 3) / (9 * Masse * G * ttab))
    v_test = H(ttab) * r - (r * p) / (3 * (1 + delta))
    return _v, v_test


def Vitesse_alter(r, dt=1e-5):  # VITESSE CALCULEE DIRECTEMENT A PARTIR DE R ET dt
    dr = np.diff(r)
    _v = dr / dt
    return np.pad(_v, (0, 1), mode='constant')  # numpy pad pour avoir len(v) == len(r)


def energie_cin(_v):
    return (3 * Masse * _v ** 2) / 10


def energie_pot(r_vir):
    return -(3 * G * Masse ** 2) / (5 * r_vir)
# Calcul du rayon viriel en prenant la moitié de son rayon max (méthode 1)


R = ((3 * Masse) / (4 * np.pi * rho_m(ttab) * (1 + delta))) ** (1 / 3)
Rmax, imax = R.max(), R.argmax()
tmax = ttab[imax]
Rvir = Rmax / 2

print("Rayon maximum:", Rmax, "atteint à:", tmax)
print("Le rayon viriel est:", Rvir)
# On ajuste les valeurs de la surdensité une fois R_vir atteint

tvir1 = t_viriel(0, len(delta), 1e-6, R, Rvir)
rho_vir = (3 * Masse) / (4 * np.pi * Rvir ** 3)  # Densité volumique de virialisation
delta_vir = rho_vir / rho_m(tvir1) - 1  # Surdensité viriel
delta[delta > delta_vir] = delta_vir  # Obsolète
indice = delta.argmax()
print("Indice de la virialisation :", indice)
print("Temps viriel calculé avec la fonction", tvir1)

# Rayon et surdensité finaux
delta_final = surdensite(delta, indice)
R_final = ((3 * Masse) / (4 * np.pi * rho_m(ttab) * (1 + delta_final))) ** (1 / 3)

# Calcul des vitesses via méthodes numériques et formule analytique
v, vtest = vitesse(R_final)
v_num = Vitesse_alter(R_final)

# Calcul des énergies - On regarde où est-ce qu'elles se croisent pour déterminer le rayon viriel (méthode 2)
Epp = energie_pot(R_final)
Ec = energie_cin(v)
Ectest = energie_cin(v_num)
Epvir = energie_pot(Rvir)

surd_viriel_finale = 1 + delta_vir * (a(teff) / a(tvir1)) ** 3
print("Surdensité viriel finale :", surd_viriel_finale)
print("Surdensité viriel finale attendue:", 18 * np.pi ** 2)
# print(delta_final[410] + 1)  # Ici l'indice n'est valable que pour une certaine valeur de N (pour N=5000)

if PLOT:
    plt.figure()
    plt.plot(ttab, 2 * Ectest, label="Energie cinétique en fonction de t")
    plt.plot(ttab, abs(Epp), label="Energie potentielle en fonction de t")
    plt.plot(tvir1, Epvir, "x", label='Energie potentielle au rayon viriel')

    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    plt.figure()
    # plt.plot(tmax,Rmax,'x',label="Rayon max")
    plt.plot(tvir1, Rvir, 'x', label="Rayon Viriel")
    plt.plot(tvir1, delta_vir, "+", label="Surdensité viriel")
    plt.plot(ttab, R_final, label="Rayon")
    plt.plot(ttab, delta_final, label="Surdensité")
    plt.plot(tmax, Rmax, 'x', label="Rayon max")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    plt.figure()
    plt.plot(ttab, abs(v), label="Vitesse théorique")
    plt.plot(ttab, abs(v_num), label="Vitesse numérique")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    plt.show()
