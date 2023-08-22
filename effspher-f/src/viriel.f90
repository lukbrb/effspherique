module viriel
<<<<<<< HEAD
    use cosmofunc, only: H
    implicit none
    
contains
    function vitesse(r) result(v)
        real, intent(in) :: r(:)
        real :: v
        v = 1.
    end function
=======
    implicit none
    use cosmofunc, only: G
    private

    real, parameter :: masse = 1E16

contains
    
    function vitesse(r) result(v)
        real, intent(in) :: r
        real :: v
        v = 1.
    end function vitesse

    function energie_cin(r) result(ec)
        real, intent(in) :: r
        real :: ec, v
        v = vitesse(r)
        ec = (3 * masse * v ** 2) / 10
    end function energie_cin

    function energie_pot(r) result(ep)
        real, intent(in) :: r
        real :: ep
        ep = - (3 * G * masse ** 2) / (5 * r)
    end function energie_pot
>>>>>>> b40897af6c5e58a1d0b28c9b474b92ff975bc41b
end module viriel