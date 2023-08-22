module viriel
    use cosmofunc, only: H
    implicit none
    
contains
    function vitesse(r) result(v)
        real, intent(in) :: r(:)
        real :: v
        v = 1.
    end function
end module viriel