module paramtune
    use solvers, only: rk4
    use cosmofunc, only: t_max, ti, eq_diff
    implicit none

    private
    public :: cond_init, iterate_on_param
contains

    real function cond_init(d_min, d_max, tolerance, delta_eff, age_univers) result(surd_mini)
        ! On cherche ici la surdensité minimale telle que le temps d'effondrement soit égal à celui de l'âge de l'univers,
        ! dans un univers d'Einstein - De Sitter.
        real, intent(in) :: d_min, d_max, tolerance, delta_eff, age_univers
        real, dimension(2) :: t_and_delta
        real :: inter, milieu, temps_eff_min, temps_eff_max, temps_eff_moy
        real :: d_i_min, d_i_max 

        d_i_min = d_min
        d_i_max = d_max
        ! Chercher le temps d'effondrement pour d_i_min
        t_and_delta = rk4(d_i_min, ti, t_max, 1e-5, delta_eff)
        temps_eff_max = t_and_delta(1)

        if (temps_eff_max < age_univers) then
            print*, "Surdensité minimum trop grande"
        end if
        ! Chercher le temps d'effondrement pour d_i_max
        t_and_delta = rk4(d_i_max, ti, t_max, 1e-5, max_density=delta_eff)
        temps_eff_min = t_and_delta(1)
        if (age_univers < temps_eff_min) then
            print *, "Surdensité maximum trop petite"
        end if
        inter = abs(d_i_max - d_i_min)
        do while (inter > tolerance)
            milieu = (d_i_min + d_i_max) / 2
            t_and_delta = rk4(milieu, ti, t_max, 1e-5, max_density=delta_eff)
            temps_eff_moy = t_and_delta(1)
            if (temps_eff_moy > age_univers) then
                d_i_min = milieu
            else
                d_i_max = milieu
            end if
            inter = d_i_max - d_i_min
        end do
        surd_mini = (d_i_max + d_i_min) / 2

    end function cond_init

    function iterate_on_param(type, array, save_res) result(res)
        integer, intent(in) :: type
        logical, intent(in) :: save_res
        real, dimension(:), intent(in) :: array
        
        real, dimension(1000) :: res
        res(:) = 1
    end function iterate_on_param
end module paramtune
