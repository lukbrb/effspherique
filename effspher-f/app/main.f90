program main
  use solvers, only: rk2, rk4, euler, rk4_write
  use cosmofunc, only: di, ti, t_max, age_univers, PI
  use paramtune, only: cond_init, iterate_on_param
  use viriel, only: surd_vir, surd_finale, surd_ta
  
  implicit none
  ! Régime non-linéaire
  ! ----------
  ! Surdensité volte-face: 4.552053651856895
  ! Surdensité volte-face théorique: 4.551652475612764
  ! Différence relative:  0.008814%
  ! ----------
  ! Surdensité viriel finale: 177.6523156109369
  ! Surdensité viriel théorique: 177.65287921960845
  ! Différence relative:  0.000317%
  ! ----------
  
  real :: dt = 1e-5
  real :: max_density = 1e4
  real :: surdensite_mini
  real, dimension(2) :: surdensite_vir, sta
  real, dimension(1) :: surd_fin

  ! resulte = euler(4. * di, ti, t_max, dt, max_density)
  ! result2 = rk2(4. * di, ti, t_max, dt, max_density)
  ! result4 = rk4_write(4. * di, ti, t_max, dt, max_density)
  ! print *, 'Résultat Euler explicit :', resulte
  ! print *, 'Résultat RK2 :', result2
  ! print *, 'Résultat RK4 :', result4

  ! ----------------------------
  ! Itération sur les paramètres
  ! ----------------------------

  
  ! surdensite_mini = cond_init(1e-4, 1e-2, 1e-7, 1e4, age_univers)
  ! write(*, '(A, F12.6)') 'Surdensité minimum pour t_eff=âge univers:', surdensite_mini

  surdensite_vir = surd_vir()
  write(*, '(A, F12.6)') 'Surdensité viriel calculée:', surdensite_vir(2)

  sta = surd_ta()
  write(*, '(A, F12.6)') 'Surdensité volte-face calculée:', sta(2)
  write(*, '(A, F12.6)') 'Surdensité volte-face théorique:', (9./16.) * PI**2 - 1.

  surd_fin = surd_finale()
  write(*, '(A, F12.6)') 'Surdensité finale calculée:', surd_fin
  write(*, '(A, F12.6)') 'Surdensité finale théorique:', 18. * PI**2
end program main
