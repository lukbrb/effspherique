module paramtune
    implicit none

    private
    public :: cond_init, iterate_on_param
contains

    real function cond_init(d_i_min, d_i_max, tolerance, delta_eff, t_eff) result(surd_mini)
    real, intent(in) :: d_i_min, d_i_max, tolerance, delta_eff, t_eff
    surd_mini = 0.1
    end function cond_init

    function iterate_on_param(type, array, save_res) result(res)
        integer, intent(in) :: type
        logical, intent(in) :: save_res
        real, dimension(:), intent(in) :: array
        
        real, dimension(1000) :: res
        res(:) = 1
    end function iterate_on_param
end module paramtune
