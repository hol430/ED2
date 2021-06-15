Module mend_exchange
  implicit none

  real, parameter :: pom_burial_rate = 1./(365.*2.*24)  ! 1/hr
  real, parameter :: dom_burial_rate = 1./240.  ! 1/hr
  real, parameter :: som_depth = 0.2 ! m

Contains

  subroutine mend_plant2som_exchange(npom, input_pom_c, input_pom_n, input_pom_p, &
       input_dom_c, input_dom_n, input_dom_p, &
       input_nh4, input_no3, input_psol, plant_input_C, plant_input_N, &
       plant_input_P)
    use decomp_coms, only: litter_partition
    implicit none
    
    integer, intent(in) :: npom
    real, dimension(npom), intent(out) :: input_pom_c
    real, intent(out) :: input_dom_c
    real, dimension(npom), intent(out) :: input_pom_n
    real, intent(out) :: input_dom_n
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, dimension(npom), intent(out) :: input_pom_p
    real, intent(out) :: input_dom_p
    real, intent(out) :: input_psol
    real, dimension(4), intent(in) :: plant_input_C
    real, dimension(4), intent(in) :: plant_input_N
    real, dimension(4), intent(in) :: plant_input_P
    integer :: iplant, ipom

    input_dom_c=0;input_dom_n=0;input_dom_p=0
    input_pom_c(:)=0;input_pom_n(:)=0;input_pom_p(:)=0

    do iplant = 1, 4
       input_dom_c = input_dom_c + litter_partition(1,iplant) * plant_input_C(iplant)
       input_dom_n = input_dom_n + litter_partition(1,iplant) * plant_input_N(iplant)
       input_dom_p = input_dom_p + litter_partition(1,iplant) * plant_input_P(iplant)
    enddo

    do ipom = 1, npom
       do iplant = 1,4
          input_pom_c(ipom) = input_pom_c(ipom) + litter_partition(1+ipom,iplant) * &
               plant_input_C(iplant)
          input_pom_n(ipom) = input_pom_n(ipom) + litter_partition(1+ipom,iplant) * &
               plant_input_N(iplant)
          input_pom_p(ipom) = input_pom_p(ipom) + litter_partition(1+ipom,iplant) * &
               plant_input_P(iplant)
       enddo
    enddo
       
    ! all in gC,N,P/m2/s

    input_nh4 = 0. 
    input_no3 = 0.
    input_psol = 0.

    return
  end subroutine mend_plant2som_exchange

  subroutine zero_exchange_vars(evars)
    use mend_state_vars, only: exchange_vars, npom
    implicit none
    type(exchange_vars) :: evars
    integer :: ipom

    do ipom = 1, npom
       evars%pom_c(ipom) = 0.
       evars%pom_n(ipom) = 0.
       evars%pom_p(ipom) = 0.
    enddo

    evars%dom_c = 0.
    evars%dom_n = 0.
    evars%dom_p = 0.
    
    evars%nh4 = 0.
    evars%no3 = 0.
    evars%psol = 0.

    return
  end subroutine zero_exchange_vars

  subroutine inc_exchange_vars(esum, einc)
    use mend_state_vars, only: exchange_vars, npom
    implicit none
    type(exchange_vars) :: esum
    type(exchange_vars) :: einc
    integer :: ipom

    do ipom = 1, npom
       esum%pom_c(ipom) = esum%pom_c(ipom) + einc%pom_c(ipom)
       esum%pom_n(ipom) = esum%pom_n(ipom) + einc%pom_n(ipom)
       esum%pom_p(ipom) = esum%pom_p(ipom) + einc%pom_p(ipom)
    enddo

    esum%dom_c = esum%dom_c + einc%dom_c
    esum%dom_n = esum%dom_n + einc%dom_n
    esum%dom_p = esum%dom_p + einc%dom_p
    
    esum%nh4 = esum%nh4 + einc%nh4
    esum%no3 = esum%no3 + einc%no3
    esum%psol = esum%psol + einc%psol

    return
  end subroutine inc_exchange_vars

  subroutine mend_som2canopy_exchange(d_co2_lost, slden, consts, d_can_co2, &
       d_co2budget_storage,ccapcani)
    use mend_consts_coms, only: decomp_consts
    implicit none

    type(decomp_consts) :: consts
    real, intent(in) :: slden
    real, intent(in) :: d_co2_lost
    real :: d_co2_lost_units
    real(kind=8), intent(inout) :: d_co2budget_storage
    real(kind=8), intent(inout) :: d_can_co2
    real(kind=8), intent(in) :: ccapcani

    ! gC/kgSoil/s
    d_co2_lost_units = d_co2_lost
    ! gC/m3Soil/s
    d_co2_lost_units = d_co2_lost_units * slden
    ! gC/m2Soil/s
    d_co2_lost_units = d_co2_lost_units * consts%eff_soil_depth
    ! molC/m2Soil/s
    d_co2_lost_units = d_co2_lost_units / 12.
    ! umolC/m2Soil/s
    d_co2_lost_units = d_co2_lost_units * 1.0e6

    d_can_co2 = d_can_co2 + d_co2_lost_units * ccapcani
    d_co2budget_storage = d_co2budget_storage + d_co2_lost_units

    return
  end subroutine mend_som2canopy_exchange

end Module mend_exchange