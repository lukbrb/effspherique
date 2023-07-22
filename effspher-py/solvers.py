import numpy as np


def euler(init_conds, function, t_max, dt=1e-3, max_density=np.inf):
    """
    Fonction qui calcule l'algorithme d'Euler explicit (ou implicit ?).
    -------------
    Arguments:
        - init_conds: tuple (delta_0, p_0, t_0) des conditions initiales
        - function : fonction directrice du système. Supposée ici f(x,t)
        - t_max : Borne supérieure du domaine temporel (d'intégration)
        - dt : le pas d'intégration
        - max_density : valeur maximale autorisée pour la densité
    :return delta_f, t_f
    """

    delta, p, t = init_conds
    while t <= t_max and delta <= max_density:
        step = function([delta, p], t) * dt
        p += step
        delta += p * dt

        t += dt

    return delta, t


def rk2(init_conds, function, t_max, dt=1e-3, max_density=np.inf):
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

    delta, p, t = init_conds
    while t <= t_max and delta <= max_density:
        k1 = function([delta, p], t)
        k2 = function([delta, p + dt * k1], t + dt)
        p += 0.5 * dt * (k1 + k2)
        delta += p * dt

        t += dt

    return delta, t


def rk4(init_conds, function, t_max, dt=1e-3, max_density=np.inf):
    """
        Fonction qui calcule l'algorithme de Runge-Kutta 4.
        -------------
        Arguments:
            - init_conds: tuple (delta_0, p_0, t_0) des conditions initiales
            - function : fonction directrice du système. Supposée ici f(x,t)
            - t_max : Borne supérieure du domaine temporel (d'intégration)
            - dt : le pas d'intégration
            - max_density : valeur maximale autorisée pour la densité
        :return delta_f, t_f
        """

    delta, p, t = init_conds
    while t <= t_max and delta <= max_density:
        k1 = function([delta, p], t)
        k2 = function([delta, p + (dt * k1) / 2], t + dt / 2)
        k3 = function([delta, p + (dt * k2) / 2], t + dt / 2)
        k4 = function([delta, p + (dt * k3) / 2], t + dt)
        p += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        delta += p * dt

        t += dt

    return delta, t


if __name__ == '__main__':
    from cosmofunc import H, tf, ti, milliard_annee, eq_diff, eq_diff_lin
    from nonlinear import rk4 as rk4_nlin
    from linear import rk4 as rk4_lin

    surd_mini = 0.001777656936645508
    tf /= milliard_annee
    ti /= milliard_annee
    print(f"Temps initial: {ti : .2f} milliards d'années")
    print(f"Temps initial: {tf : .2f} milliards d'années")
    init = (4 * surd_mini, 4 * surd_mini * H(ti), ti)
    res1 = euler(init, eq_diff, tf, dt=1e-5, max_density=1e4)
    print(res1)
    res2 = rk2(init, eq_diff, tf, dt=1e-5, max_density=1e4)
    print(res2)
    res3 = rk4(init, eq_diff, tf, dt=1e-5, max_density=1e4)
    print(res3)
    # Comparaison avec les vieilles fonctions
    res4 = rk4_nlin(4 * surd_mini, 1e4)
    res5 = rk4_lin(4 * surd_mini, 1e4)
    print("-" * 50)
    print("Ancienne RK4 non-linéaire:", res4[0][-1, 0])
    print("Ancienne RK4 linéaire:", res5[-1, 0])
    print("-" * 50)
    print("Nouvelle RK4 non-linéaire:", res3)
    res6 = rk4(init, eq_diff_lin, tf, dt=1e-5, max_density=1e4)
    print("Nouvelle RK4 linéaire:", res6)
