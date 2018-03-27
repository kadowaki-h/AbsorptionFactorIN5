! test_AbsorptionFactorIN5.f90
!   H Kadowaki
!
! 6-CPU parallel computation using OpenMP
! ifort -qopenmp -static -warn all -O -o test_AbsorptionFactorIN5static.ex mod_AbsorptionFactorIN5.f90 test_AbsorptionFactorIN5.f90
! export OMP_NUM_THREADS=6
! time ./test_AbsorptionFactorIN5static.ex < AbsorptionFactorIN5_test_in > AbsorptionFactorIN5_test_out &
!
! single CPU, no OpenMP
! ifort -warn all -O -o test_AbsorptionFactorIN5.ex mod_AbsorptionFactorIN5.f90 test_AbsorptionFactorIN5.f90
!
! OpenMP
! gfortran -fopenmp -Wall -O -o test_AbsorptionFactorIN5.ex mod_AbsorptionFactorIN5.f90 test_AbsorptionFactorIN5.f90
! single CPU, no OpenMP
! gfortran -Wall -O -o test_AbsorptionFactorIN5.ex mod_AbsorptionFactorIN5.f90 test_AbsorptionFactorIN5.f90
!
!-----
program main
!
!$ use omp_lib
  use mod_absorp, only : absorp_init, absorp_calc
  use mod_in5, only : in5_init, read_par_file, read_spe_file, absorp_set_drct, abs_corr_SQE_err_SQE &
  & , write_spe_file, dealloc_in5, Delta_Omega, nABS, non_zero_SQE_tab, UKF_tab, mu_out_tab, ABSFAC_tab &
  & , UKI, mu_in
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
  character(256) :: TITLE
  integer(4) :: nf_spe_in, nf_spe_out, nf_par
  character(128) :: spe_file_name_in, spe_file_name_out, par_file_name
  real(DP) :: Omega, Offset
  real(DP) :: UKF(3),mu_out,ABSFAC
  integer(4) :: i, i_begin, i_end, ie, id, j
  integer(4), parameter :: parallel_do_unit=1000
!  parallel_do_unit=processing unit of  !$OMP parallel do  e.g. 200 * OMP_NUM_THREADS
!
!$ call omp_set_dynamic(.false.)
!$ write(*,'(A)')'  omp_set_dynamic(.false.)'
!$ write(*,'(A,I10)')'  omp_get_max_threads() = ',omp_get_max_threads()
!$ write(*,'(A,I10)')'  omp_get_num_procs() = ',omp_get_num_procs()
!$ write(*,'(A,L10)')'  omp_get_nested() = ',omp_get_nested()
!$ write(*,'(A,L10)')'  omp_get_dynamic() = ',omp_get_dynamic()
!
  write(*,'(A,I10)')'  parallel_do_unit =',parallel_do_unit
  READ(*,'(A256)') TITLE
  WRITE(*,'(/,2x,A)') TRIM(TITLE)
!
  call in5_init
  call absorp_init    ! read and initialize parameters of spectrometer
!
  read(*,*)par_file_name
  write(*,'(/,A,A)')'  par file name = par_file_name= ',TRIM(par_file_name)
  nf_par=31
  open(nf_par,file=par_file_name)
    call read_par_file(nf_par)  ! get nd, phi(nd), azim(nd)
  close(nf_par)
!
 read(*,*) Offset
 write(*,'(/,A,G15.7)')'  Offset = Omega_0 = ',Offset
 
 do j=1,90000
!
! see vectors.pdf
  read(*,*,end=990) Omega, spe_file_name_in, spe_file_name_out
    if(Omega < -3600.0d0) exit
    if(spe_file_name_in(1:3) == 'end' .or. spe_file_name_out(1:3) == 'end' ) exit
  Delta_Omega = Omega - Offset
  write(*,'(/,A,3G15.7)')'  Omega, Delta_Omega, psi=', Omega, Delta_Omega, -Delta_Omega
  write(*,'(A,A)')'  spe_file_name_in  =',TRIM(spe_file_name_in)
  write(*,'(A,A)')'  spe_file_name_out =',TRIM(spe_file_name_out)
!
  nf_spe_in  = 31
  nf_spe_out = 32
  open(nf_spe_in ,file=spe_file_name_in )
  open(nf_spe_out,file=spe_file_name_out)
!
  call read_spe_file(nf_spe_in,nf_spe_out)
  close(nf_spe_in)
!
  call absorp_set_drct
!
  do i_begin = 1, nABS, parallel_do_unit
    i_end   = MIN( i_begin + parallel_do_unit -1 , nABS )
!
!$OMP parallel do private(i,id,ie,UKF,mu_out,ABSFAC)
!!$ write(*,*)'  omp_get_num_threads(), omp_get_thread_num() = ',omp_get_num_threads(),omp_get_thread_num()
    do i = i_begin, i_end
      ie = non_zero_SQE_tab(1,i)
      id = non_zero_SQE_tab(2,i)
      UKF(1:3) = UKF_tab(1:3,id)
      mu_out = mu_out_tab(ie)
      call absorp_calc(UKI,UKF,mu_in,mu_out,ABSFAC)
      ABSFAC_tab(i) = ABSFAC
!!$ write(*,*)'  i = ',i,'    omp_get_thread_num() = ',omp_get_thread_num()
    end do
!$OMP end parallel do
!
  end do
!
  call abs_corr_SQE_err_SQE
  call write_spe_file(nf_spe_out)
  close(nf_spe_out)
  call dealloc_in5
!
 end do
 990 continue
!
end program main
!

