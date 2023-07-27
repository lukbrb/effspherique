module solvers
    implicit none
    private

    public :: rk4, H

contains
    function H(t) result(hubble_cst)
        real, intent(in) :: t
        real :: hubble_cst
        hubble_cst = 2 / (3 * t)
    end function H

    function eq_diff(d, p, t) result(val_eq)
        real, intent(in) :: d, p, t
        real :: val_eq
        val_eq = -2 * H(t) * p + 1.5 * (H(t) ** 2) * d
    end function eq_diff

    function rk4(di, pi, ti, t_max, dt, max_density) result(teff)
        implicit none
        real, intent(in) :: di, pi, ti, t_max, dt, max_density
        real :: delta, p, t, teff
        real :: k1, k2, k3, k4

        delta = di
        p = pi 
        t = ti
        do while(t <= t_max .and. delta <= max_density)
            k1 = eq_diff(delta, p, t)
            k2 = eq_diff(delta, p + (dt * k1) / 2, t + dt / 2)
            k3 = eq_diff(delta, p + (dt * k2) / 2, t + dt / 2)
            k4 = eq_diff(delta, p + (dt * k3) / 2, t + dt)
            p = p + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            delta = delta + p * dt

            t = t + dt
            print *, 't=', t, 'delta=', delta
        end do

        teff = t
    end function rk4  
end module solvers