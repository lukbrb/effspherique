module solvers
    use cosmofunc, only: H, eq_diff
    implicit none
    private

    public :: euler, rk2, rk4, rk4_write
! TODO: Transformer les solveurs en subroutines
contains
    function euler(di, ti, t_max, dt, max_density) result(t_and_delta)
    ! Fonction qui calcule l'algorithme d'Euler explicit.
    !   -------------
    !   Arguments:
    !       - init_conds: tuple (delta_0, t_0) des conditions initiales
    !       - function : fonction directrice du système. Supposée ici f(x,t)
    !       - t_max : Borne supérieure du domaine temporel (d'intégration)
    !       - dt : le pas d'intégration
    !       - max_density : valeur maximale autorisée pour la densité
    !   :return t, delta
        real(kind = 8), intent(in) :: di, ti, t_max, dt, max_density
        real(kind = 8) :: delta, p, t, step
        real(kind = 8), dimension(2) :: t_and_delta

        delta = di
        p = H(ti) * di
        t = ti
        do while(t <= t_max .and. delta <= max_density)
            step = eq_diff(delta, p, t) * dt
            p = p + step
            delta = delta + p * dt

            t = t + dt
        end do
        t_and_delta = [t, delta]
    end function euler

    function rk2(di, ti, t_max, dt, max_density) result(t_and_delta)
        ! Fonction qui calcule l'algorithme de Runge-Kutta 2.
        ! -------------
        ! Arguments:
        !     - init_conds: tuple (delta_0, t_0) des conditions initiales
        !     - function : fonction directrice du système. Supposée ici f(x,t)
        !     - t_max : Borne supérieure du domaine temporel (d'intégration)
        !     - dt : le pas d'intégration
        !     - max_density : valeur maximale autorisée pour la densité
        ! :return t, delta
        real(kind = 8), intent(in) :: di, ti, t_max, dt, max_density
        real(kind = 8) :: delta, p, t
        real(kind = 8):: k1, k2
        real(kind = 8), dimension(2) :: t_and_delta

        delta = di
        p = H(ti) * di
        t = ti
        do while(t <= t_max .and. delta <= max_density)
            k1 = eq_diff(delta, p, t)
            k2 = eq_diff(delta, p + dt * k1, t + dt)
            p = p + 0.5 * dt * (k1 + k2)
            delta = delta + p * dt

            t = t + dt
            ! print *, 't=', t, 'delta=', delta
        end do
        t_and_delta = [t, delta]
    end function rk2

    function rk4(di, ti, t_max, dt, max_density) result(t_and_delta)
        ! Fonction qui calcule l'algorithme de Runge-Kutta 2.
        ! -------------
        ! Arguments:
        !     - init_conds: tuple (delta_0, t_0) des conditions initiales
        !     - function : fonction directrice du système. Supposée ici f(x,t)
        !     - t_max : Borne supérieure du domaine temporel (d'intégration)
        !     - dt : le pas d'intégration
        !     - max_density : valeur maximale autorisée pour la densité
        ! :return t, delta
        implicit none
        real(kind = 8), intent(in) :: di, ti, t_max, dt, max_density
        real(kind = 8) :: delta, p, t
        real(kind = 8) :: k1, k2, k3, k4
        real(kind = 8), dimension(2) :: t_and_delta

        delta = di
        p = H(ti) * di
        t = ti
        do while(t <= t_max .and. delta <= max_density)
            k1 = eq_diff(delta, p, t)
            k2 = eq_diff(delta, p + (dt * k1) / 2, t + dt / 2)
            k3 = eq_diff(delta, p + (dt * k2) / 2, t + dt / 2)
            k4 = eq_diff(delta, p + (dt * k3) / 2, t + dt)
            p = p + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            delta = delta + p * dt

            t = t + dt
            ! print *, 't=', t, 'delta=', delta
        end do
        t_and_delta = [t, delta]
    end function rk4

    subroutine rk4_write(di, ti, t_max, dt, max_density)
        ! Fonction qui calcule l'algorithme de Runge-Kutta 4.
        ! -------------
        ! Arguments:
        !     - init_conds: tuple (delta_0, t_0) des conditions initiales
        !     - function : fonction directrice du système. Supposée ici f(x,t)
        !     - t_max : Borne supérieure du domaine temporel (d'intégration)
        !     - dt : le pas d'intégration
        !     - max_density : valeur maximale autorisée pour la densité
        ! :return t, delta
        real(kind = 8), intent(in) :: di, ti, t_max, dt, max_density
        real(kind = 8) :: delta, ddelta, t
        real(kind = 8) :: k1, k2, k3, k4

        integer :: io
        open(newunit=io, file="results/rk4-results.dat", status="replace", action="write")

        delta = di
        ddelta = H(ti) * di
        t = ti
        do while(t <= t_max .and. delta <= max_density)
            k1 = eq_diff(delta, ddelta, t)
            k2 = eq_diff(delta, ddelta + (dt * k1) / 2, t + dt / 2)
            k3 = eq_diff(delta, ddelta + (dt * k2) / 2, t + dt / 2)
            k4 = eq_diff(delta, ddelta + (dt * k3) / 2, t + dt)

            ddelta = ddelta + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            delta = delta + ddelta * dt
            t = t + dt

            write(io, *) t, ddelta, delta
        end do
        close(io)
    end subroutine rk4_write 
end module solvers