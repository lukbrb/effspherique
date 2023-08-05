import numpy as np


Mp = 3.08567758 * 1e22  # Valeur d'un mégaparsec
H0 = 7 * 1e4 / Mp  # Constante d'Hubble Actuelle
G = 6.6742 * 1e-11  # Constante de gravitation
rho_crit = (3 * H0 ** 2) / (8 * np.pi * G)  # Densité critique
sig = 1  # Densité de l'univers
milliard_annee = 3600 * 24 * 365 * 1e9
age_univ = (2 / (3 * H0)) / milliard_annee  # L'âge de l'univers
surd_mini = 0.0017104343414306644  # Surdensité mini pour que teff <= age_univers
# CONDITIONS INITIALES
# On veut ti très petit devant 1/H0, or 1/H0= 4.4e17s (=13 952 308 472 ans)
# On peut donc prendre ti=300000 ans par exemple

ti = (300000 * 365 * 24 * 3600) / milliard_annee
tf = (1 / H0) / milliard_annee


def H(t):
    return 2 / (3 * t)


def a(t):
    return (1.5 * H0 * t) ** (2 / 3)


def rho_m(t):
    return (sig * rho_crit) / (a(t) ** 3)


def rho(r, m):
    return (3 * m) / (4 * np.pi * r**3)


def eq_diff_lin(d, p, t):
    return -2 * H(t) * p + 1.5 * (H(t) ** 2) * d


def eq_diff(d, p, t):
    return -2 * H(t) * p + 1.5 * (H(t) ** 2) * d * (1 + d) + (4 / 3) * ((p ** 2) / (1 + d))
