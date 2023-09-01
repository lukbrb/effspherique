program main
  use solvers, only: rk2, rk4, euler, rk4_write
  use cosmofunc, only: di, ti, t_max, age_univers
  use paramtune, only: cond_init, iterate_on_param
  use viriel, only: surd_vir, surd_finale
  
  implicit none

  real :: dt = 1e-5
  real :: max_density = 1e4
  real :: surdensite_mini
  real, dimension(2) :: resulte, result2, result4, surdensite_vir
  real, dimension(1) :: surd_fin

  resulte = euler(4. * di, ti, t_max, dt, max_density)
  result2 = rk2(4. * di, ti, t_max, dt, max_density)
  result4 = rk4_write(4. * di, ti, t_max, dt, max_density)
  print *, 'Résultat Euler explicit :', resulte
  print *, 'Résultat RK2 :', result2
  print *, 'Résultat RK4 :', result4

  ! ----------------------------
  ! Itération sur les paramètres
  ! ----------------------------

  
  surdensite_mini = cond_init(1e-4, 1e-2, 1e-7, 1e4, age_univers)
  print *, 'Surdensité minimum pour t_eff=âge univers:', surdensite_mini
  surdensite_vir = surd_vir()
  print *, surdensite_vir
  surd_fin = surd_finale()
  write(*, *), 'surd_fin', surd_fin
end program main
