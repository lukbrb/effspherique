from cosmofunc import H, ti, ttab, eq_diff, dt, milliard_annee, tf
import numpy as np
import time

surd_mini = 0.001777656936645508

# TODO: Avoir des fonctions RK2, RK4, etc indépendantes des deltas, tq:
# def solver(func, ordre, points), puis une autre fonction qui résout le problème spécifique


def rk2(_delta_i, _delta_eff):
    tps1 = time.time()
    pi = _delta_i * H(ti)
    delta_p = np.zeros((len(ttab), 2))
    delta_p[0, 0], delta_p[0, 1] = _delta_i, pi
    for n in range(1, len(ttab)):
        t = ttab[n]
        if (delta_p[n - 1, 0] and delta_p[n - 1, 1]) < _delta_eff:  # Condition pour s'assurer que delta ne diverge pas
            k1 = -2 * H(t) * delta_p[n - 1, 1] + 1.5 * (H(t) ** 2) * delta_p[n - 1, 0] * (1 + delta_p[n - 1, 0]) \
                 + (4 / 3) * ((delta_p[n - 1, 1] ** 2) / (1 + delta_p[n - 1, 0]))
            k2 = eq_diff([delta_p[n - 1, 0], delta_p[n - 1, 1] + dt * k1], t + dt)
            delta_p[n, 1] = delta_p[n - 1, 1] + 0.5 * dt * (k1 + k2)
            delta_p[n, 0] = delta_p[n - 1, 0] + delta_p[n - 1, 1] * dt
        else:
            delta_p[n, 0] = _delta_eff + 1
            delta_p[n, 1] = delta_p[n - 1, 1]
    delta_p[delta_p[:, 0] > _delta_eff, 0] = _delta_eff
    argmax = delta_p[:, 1].argmax()
    delta_p[delta_p[:, 1] > delta_p[argmax - 1, 1], 1] = delta_p[argmax - 1, 1]
    if delta_p[-1, 0] != _delta_eff:
        t_eff = 0
        print("La surdensité ne s'est pas encore effondrée ")
    else:
        argmax = delta_p[:, 0].argmax()
        t_eff = ttab[argmax]
        print("Effondrement à ", t_eff / milliard_annee, "milliards d'années, indice", argmax)
    tps2 = time.time()
    print("Temps de calcul RK2 :", tps2 - tps1, "secondes")
    return delta_p, t_eff


# Implémentation de Runge-Kutta 4

def rk4(_delta_i, _delta_eff):
    tps1 = time.time()
    pi = _delta_i * H(ti)
    delta_p = np.zeros((len(ttab), 2))
    delta_p[0, 0], delta_p[0, 1] = _delta_i, pi
    for n in range(1, len(ttab)):
        t = ttab[n]
        if (delta_p[n - 1, 0] and delta_p[n - 1, 1]) < _delta_eff:  # Condition pour s'assurer que delta ne diverge pas
            k1 = -2 * H(t) * delta_p[n - 1, 1] + 1.5 * (H(t) ** 2) * delta_p[n - 1, 0] * (1 + delta_p[n - 1, 0]) \
                 + (4 / 3) * ((delta_p[n - 1, 1] ** 2) / (1 + delta_p[n - 1, 0]))
            k2 = eq_diff([delta_p[n - 1, 0], delta_p[n - 1, 1] + (dt * k1) / 2], t + dt / 2)
            k3 = eq_diff([delta_p[n - 1, 0], delta_p[n - 1, 1] + (dt * k2) / 2], t + dt / 2)
            k4 = eq_diff([delta_p[n - 1, 0], delta_p[n - 1, 1] + (dt * k3) / 2], t + dt)
            delta_p[n, 1] = delta_p[n - 1, 1] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            delta_p[n, 0] = delta_p[n - 1, 0] + delta_p[n - 1, 1] * dt
        else:
            delta_p[n, 0] = _delta_eff + 1
            delta_p[n, 1] = delta_p[n - 1, 1]
    delta_p[delta_p[:, 0] > _delta_eff, 0] = _delta_eff
    argmax = delta_p[:, 1].argmax()
    delta_p[delta_p[:, 1] > delta_p[argmax - 1, 1], 1] = delta_p[argmax - 1, 1]
    if delta_p[-1, 0] != _delta_eff:
        t_eff = 0
        print("La surdensité ne s'est pas encore effondrée ")
    else:
        argmax = delta_p[:, 0].argmax()
        t_eff = ttab[argmax]
        print("Effondrement à ", t_eff / milliard_annee, "milliards d'années, indice", argmax)
    tps2 = time.time()
    print("Temps de calcul RK4 :", tps2 - tps1, "secondes")
    return delta_p, t_eff


def Euler(_delta_i, _delta_eff):
    p_i = H(ti) * _delta_i
    delta_p = np.zeros((len(ttab), 2))
    delta_p[0, 0], delta_p[0, 1] = p_i, _delta_i
    for k in range(1, len(ttab)):
        t = ttab[k]
        if (delta_p[k - 1, 0] and delta_p[k - 1, 1]) < _delta_eff:
            pas_p = (-2 * H(t) * delta_p[k - 1, 0] + 1.5 * (H(t) ** 2) * delta_p[k - 1, 1] * (1 + delta_p[k - 1, 1])
                     + (4 / 3) * ((delta_p[k - 1, 0] ** 2) / (1 + delta_p[k - 1, 1]))) * dt
            delta_p[k, 0] = delta_p[k - 1, 0] + pas_p
            delta_p[k, 1] = delta_p[k - 1, 1] + delta_p[k - 1, 0] * dt

        else:
            print("Le temps d'effondrement pour δi =", _delta_i, "est de ", ttab[k] / milliard_annee,
                  "milliards d'années", "arg(t)=", k)
    return delta_p


if __name__ == '__main__':
    x, teff = rk4(4 * surd_mini, 3 * 1e4)
    # ITERATIONS SUR LES SURDENSITÉS INITIALES

    temps_eff = open("resultats/temps_eff.txt", "w")
    surdens_init = open("resultats/surdensité_initiale.txt", "w")
    # On part de la surdensité qui s'effondre à l'âge de l'univers
    surd_init = np.linspace(0.001777656936645508, 15 * 0.001777656936645508, 100)
    for delta_i in surd_init:
        z, eff = rk4(delta_i, 1e4)
        temps_eff.write("{}\n".format(str(eff)))
        surdens_init.write("{}\n".format(str(delta_i)))
    temps_eff.close()
    surdens_init.close()

    # ITERATIONS SUR SURDENSITÉ D'EFFONDREMENT

    temps_eff = open("resultats/temps_deltaeff.txt", "w")
    surdens_eff = open("resultats/Delta_eff.txt", "w")
    surd_eff = np.linspace(1e3, 1e6, 1000)
    for delta_eff in surd_eff:
        z, eff = rk4(4 * surd_mini, delta_eff)
        temps_eff.write("{}\n".format(str(eff)))
        surdens_eff.write("{}\n".format(str(delta_eff)))
    temps_eff.close()
    surdens_eff.close()

    # ITERATIONS SUR LE PAS DE TEMPS

    temps_fctN = open("resultats/temps_fctN.txt", "w")
    pas_temps = open("resultats/pas_temps.txt", "w")
    for N in range(5000, 10000):
        temps = np.linspace(ti, tf, round(N))
        dt1 = (tf - ti) / N
        z, eff = rk4(4 * surd_mini, 3 * 1e4)
        temps_fctN.write("{}\n".format(str(eff)))
        pas_temps.write("{}\n".format(str(dt1)))
    temps_fctN.close()
    pas_temps.close()  # NE FONCTIONNAIT PAS, J'AI DONC REGARDE QUELQUES VALEURS A LA MAIN
