import numpy as np


def euler(init_conds, function, t_max, dt=1e-5, max_density=np.inf):
    """
    Fonction qui calcule l'algorithme d'Euler explicit.
    -------------
    Arguments:
        - init_conds: tuple (delta_0, p_0, t_0) des conditions initiales
        - function : fonction directrice du système. Supposée ici f(x,t)
        - t_max : Borne supérieure du domaine temporel (d'intégration)
        - dt : le pas d'intégration
        - max_density : valeur maximale autorisée pour la densité
    :return delta, t
    """

    results = []
    ttab = []
    delta, p, t = init_conds
    while t <= t_max and delta <= max_density:
        step = function([delta, p], t) * dt
        p += step
        delta += p * dt

        t += dt
        results.append(delta)
        ttab.append(t)

    return results, ttab


def rk2(init_conds, function, t_max, dt=1e-5, max_density=np.inf):
    """
        Fonction qui calcule l'algorithme de Runge-Kutta 2.
        -------------
        Arguments:
            - init_conds: tuple (delta_0, p_0, t_0) des conditions initiales
            - function : fonction directrice du système. Supposée ici f(x,t)
            - t_max : Borne supérieure du domaine temporel (d'intégration)
            - dt : le pas d'intégration
            - max_density : valeur maximale autorisée pour la densité
        :return delta_f, t_f
        """
    results = []
    ttab = []
    delta, p, t = init_conds
    while t <= t_max and delta <= max_density:
        k1 = function([delta, p], t)
        k2 = function([delta, p + dt * k1], t + dt)
        p += 0.5 * dt * (k1 + k2)
        delta += p * dt

        t += dt
        results.append(delta)
        ttab.append(t)

    return results, ttab


def rk4(init_conds, function, t_max, dt=1e-5, max_density=np.inf):
    """
        Fonction qui calcule l'algorithme de Runge-Kutta 4.
        -------------
        :argument
            - init_conds: tuple (delta_0, p_0, t_0) des conditions initiales
            - function : fonction directrice du système. Supposée ici f(x,t)
            - t_max : Borne supérieure du domaine temporel (d'intégration)
            - dt : le pas d'intégration
            - max_density : valeur maximale autorisée pour la densité
        :return delta_f, t_f
        """
    results = []
    derivee = []
    ttab = []
    delta, p, t = init_conds

    while t <= t_max and delta <= max_density:
        k1 = function([delta, p], t)
        k2 = function([delta, p + (dt * k1) / 2], t + dt / 2)
        k3 = function([delta, p + (dt * k2) / 2], t + dt / 2)
        k4 = function([delta, p + (dt * k3) / 2], t + dt)
        p += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        delta += p * dt

        t += dt
        results.append(delta)
        derivee.append(p)
        ttab.append(t)

    return np.array(results), np.array(derivee), np.array(ttab)


if __name__ == '__main__':
    # TODO: Changer code pour prendre en compte modif rk4
    import matplotlib.pyplot as plt

    from cosmofunc import H, tf, ti, milliard_annee, eq_diff, eq_diff_lin

    # from nonlinear import rk4 as rk4_nlin
    # from linear import rk4 as rk4_lin

    surd_mini = 0.001777656936645508
    tf /= milliard_annee
    ti /= milliard_annee
    print(f"Temps initial: {ti : .2f} milliards d'années")
    print(f"Temps initial: {tf : .2f} milliards d'années")
    init = (4 * surd_mini, 4 * surd_mini * H(ti), ti)
    res1 = euler(init, eq_diff, tf, dt=1e-5, max_density=1e4)
    res2 = rk2(init, eq_diff, tf, dt=1e-5, max_density=1e4)
    res3 = rk4(init, eq_diff, tf, dt=1e-5, max_density=1e4)

    # Comparaison avec les vieilles fonctions
    # res4 = rk4_nlin(4 * surd_mini, 1e4)

    plt.figure()
    plt.title("Comparaison pour les solutions non-linéaires")
    plt.plot(res1[1], res1[0], label="Méthode d'Euler")
    plt.plot(res2[1], res2[0], label="Méthode RK2")
    plt.plot(res3[1], res3[0], label="Méthode RK4")
    plt.yscale('log')
    # plt.plot(res4[1]/milliard_annee, res4[0][:, 0], '--r', label="Ancienne Méthode RK4")
    plt.legend()
    plt.show()

    # res5 = rk4_lin(4 * surd_mini, 1e4)
    res6 = rk4(init, eq_diff_lin, tf, dt=1e-5, max_density=1e4)

    plt.figure()
    plt.title("Comparaison pour les solutions linéaires")
    plt.plot(res6[1], res6[0], label="Méthode RK4")
    # plt.plot(res5[1]/milliard_annee, res5[0][:, 0], label="Ancienne Méthode RK4")
    # plt.plot(res6[1], np.array(res6[1])**(2/3), label="t^{2/3}")
    # plt.yscale('log')
    plt.legend()
    plt.show()
