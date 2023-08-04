import numpy as np
from cosmofunc import H


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


def rk4_array(delta_i, _ti, function, t_max, dt=1e-5, max_density=np.inf):
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
    delta, t = delta_i, _ti
    p = delta_i * H(_ti)

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


def rk4(delta_i, _ti, function, t_max, dt=1e-5, max_density=1e4):
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

    delta, t = delta_i, _ti
    p = delta_i * H(_ti)

    while t <= t_max and delta <= max_density:
        k1 = function([delta, p], t)
        k2 = function([delta, p + (dt * k1) / 2], t + dt / 2)
        k3 = function([delta, p + (dt * k2) / 2], t + dt / 2)
        k4 = function([delta, p + (dt * k3) / 2], t + dt)
        p += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        delta += p * dt

        t += dt

    return delta, p, t


def integrateur(delta_i, _ti, function, t_max, dt=1e-5, max_density=1e4, kind='rk4', array=False):
    if kind.lower() != 'rk4':
        raise NotImplementedError(f'Intégrateur général non implémenté pour {kind}')
    if array:
        solv = rk4_array
    else:
        solv = rk4

    return solv(delta_i, _ti, function, t_max, dt, max_density)
