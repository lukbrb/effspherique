module effspher_f
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, effspher-f!"
  end subroutine say_hello
end module effspher_f
