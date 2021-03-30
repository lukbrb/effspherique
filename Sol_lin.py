import numpy as np
from main import *


def fonction(delta_et_p, ttab):
    return -2 * H(ttab) * delta_et_p[1] + 1.5 * (H(ttab) ** 2) * delta_et_p[0]


def Euler_lin(delta_i, delta_eff):
    p_i = H(ti) * delta_i
    delta_p = np.zeros((len(ttab), 2))
    delta_p[0, 0], delta_p[0, 1] = delta_i, p_i
    for n in range(1, len(ttab)):
        t = ttab[n]
        if (delta_p[n - 1, 0] and delta_p[n - 1, 1]) < delta_eff:
            pas_p = (-2 * H(t) * delta_p[n - 1, 1] + 1.5 * (H(t) ** 2) * delta_p[n - 1, 0]) * dt
            delta_p[n, 1] = delta_p[n - 1, 1] + pas_p
            delta_p[n, 0] = delta_p[n - 1, 0] + delta_p[n - 1, 1] * dt
        else:
            print("Le temps d'effondrement pour δi =", delta_i, "est de ", ttab[n] / milliard_annee,
                  "milliards d'années", "(indice:", n, ")")
    return delta_p


def rk2_lin(delta_i, delta_eff):
    pi = delta_i * H(ti)
    delta_p = np.zeros((len(ttab), 2))
    delta_p[0, 0], delta_p[0, 1] = delta_i, pi
    for n in range(1, len(ttab)):
        t = ttab[n]
        if (delta_p[n - 1, 0] and delta_p[n - 1, 1]) < delta_eff:
            k1 = -2 * H(t) * delta_p[n - 1, 1] + 1.5 * (H(t) ** 2) * delta_p[n - 1, 0]
            k2 = fonction([delta_p[n - 1, 0], delta_p[n - 1, 1] + dt * k1], t + dt)
            delta_p[n, 1] = delta_p[n - 1, 1] + 0.5 * dt * (k1 + k2)
            delta_p[n, 0] = delta_p[n - 1, 0] + delta_p[n - 1, 1] * dt
        else:
            print("Le temps d'effondrement pour δi =", delta_i, "est de ", ttab[n] / milliard_annee,
                  "milliards d'années", "indice:", n)
    return delta_p


def rk4_lin(delta_i, delta_eff):
    pi = delta_i * H(ti)
    delta_p = np.zeros((len(ttab), 2))
    delta_p[0, 0], delta_p[0, 1] = delta_i, pi
    for n in range(1, len(ttab)):
        t = ttab[n]
        if (delta_p[n - 1, 0] and delta_p[n - 1, 1]) < delta_eff:
            k1 = -2 * H(t) * delta_p[n - 1, 1] + 1.5 * (H(t) ** 2) * delta_p[n - 1, 0]
            k2 = fonction([delta_p[n - 1, 0], delta_p[n - 1, 1] + (dt * k1) / 2], t + dt / 2)
            k3 = fonction([delta_p[n - 1, 0], delta_p[n - 1, 1] + (dt * k2) / 2], t + dt / 2)
            k4 = fonction([delta_p[n - 1, 0], delta_p[n - 1, 1] + (dt * k3) / 2], t + dt)
            delta_p[n, 1] = delta_p[n - 1, 1] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            delta_p[n, 0] = delta_p[n - 1, 0] + delta_p[n - 1, 1] * dt
        else:
            temps_eff = ttab[n]
    if delta_p[-1, 0] != 0:
        temps_eff = 0
        print("La surdensité ne s'est pas encore effondrée ")

    return delta_p


def equation(delta_et_p, ttab):
    return [delta_et_p[1], -2 * H(ttab) * delta_et_p[1] + 1.5 * (H(ttab) ** 2) * delta_et_p[0]]
