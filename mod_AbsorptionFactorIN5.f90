! mod_AbsorptionFactorIN5.f90
!**********************************************************************
! Copyright (c) 2018, Hiroaki Kadowaki (kadowaki@tmu.ac.jp)
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, 
! this list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation and/or 
! other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors may 
! be used to endorse or promote products derived from this software without specific 
! prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
! POSSIBILITY OF SUCH DAMAGE.
! 
!    (See https://opensource.org/licenses/BSD-3-Clause )
!
!**********************************************************************
! A fortran90 program of calculating absorption factors
! single crystal measurements at ILL-IN5
!
! mod_AbsorptionFactorIN5.f90
!    614 lines
!
!  calculate absorption factors
!      neutron scattering    IN5 and AMATERAS
!      SPE file --> absorption corrected SPE file
!
!  sample shape is convex and approximated by surfaces of planes and cylinders
!      convex region is a region where, for every pair of points within the region, 
!      every point on the straight line segment that joins the pair of points is also within the region
!
!  ABSFAC= absorption factor
!        =(  \int_{X1}^{X2} dx \int_{Y1}^{Y2} dy \int_{Z1}^{Z2} dz
!             exp( - (\mu_{in absorp}*(k_0/k_in) + \mu_{in incoh}) * l_in )
!           * exp( - (\mu_{out absorp}*(k_0/k_out) + \mu_{out incoh}) * l_out )
!           * f(x,y,z)  )
!         / (  \int_{X1}^{X2} dx \int_{Y1}^{Y2} dy \int_{Z1}^{Z2} dz  f(x,y,z)  )
!
!  xyz = orthogonal coordinates are defined by user (cm)
!             and are used in the calculation of ABSFAC
!
!  \mu_{in absorp}*(k_0/k_in) + \mu_{in incoh}
!       = linear absorption coefficient of incoming neutron (cm^{-1})
!  \mu_{out absorp}*(k_0/k_out) + \mu_{out incoh}
!       = linear absorption coefficient of outgoing neutron (cm^{-1})
!  k_0 = 3.494158 (A^{-1})   ( v_0 = 2200 (m/s) )
!  k_in = wavevector of incoming neutron (A^{-1})
!  k_out = wavevector of outgoing neutron (A^{-1})
!  l_in = path length of incoming neutron (cm)
!  l_out = path length of outgoing neutron (cm)
!    see https://www.ncnr.nist.gov/resources/n-lengths/list.html
!
!  f(x,y,z) = 1 inside the sample
!           = 0 outside the sample
!
!  the sample shape is expressed by surfaces of planes and cyclinders
!     a point (x,y,z) is inside the sample,
!     if
!       DOT_PRODUCT(PLA(1:3,I),(/x,y,z/))  >=  PLB(I)    (I=1,NPL)
!         PLA(1:3,I) = a direction orthogonal to a plane  (no dimension)
!         PLB(I)=   (cm)
!
!       DOT_PRODUCT(delx(1:3),delx(1:3)) - DOT_PRODUCT(delx(1:3),CYLD(1:3,J))**2  <=  CYLR(J)**2    (J=1,NCYL)
!         CYLC(1:3,J) = a point on the center axis of the cylinder   (xyz coordinates (cm))
!         CYLD(1:3,J) = the unit vector along the the cylinder axis  (xyz coordinates (no dimension))
!         CYLR(J) = radius of cylinder  (xyz coordinates (cm))
!         delx(1:3) = (/x,y,z/) - CYLC(1:3,J)
!
!   Delta_E(ie) = energy transfer  (meV)    ie-th Ebin
!
! input
!
!  READ(*,*) wave_length
!    wave_length =  wave_length (A)
!            k_in=(2.0d0*pi)/wave_length (A**^{-1}) 
!
!  READ(*,*) UT(1:3)
!  READ(*,*) VT(1:3)
!     UT(1:3) = a vector // (k_{in} for Delta_Omega = 0)
!           UT(1:3) is converted to the unit vector in the program
!           xyz coordinates (no dimension)
!     VT(1:3) = a vector perpendicular to (k_{in} for Delta_Omega = 0) and in the horizontal scattering plane
!           VT(1:3) is converted to exactly orthogonal to UT(1:3) and then is normalized to the unit vector.
!           xyz coordinates (no dimension)
!     WT(1:3) = UT(1:3) x VT(1:3)  (cross product)  // vertical direction
!
!  READ(*,*) mu_in_absorp, mu_in_incoh, mu_out_absorp, mu_out_incoh
!     mu_in_absorp = \mu_{in absorp}  (cm^{-1})
!           = linear absorption coefficient by absorption of incoming neutron
!             a value for 2200m/s neutrons (k_0 = 3.494158 A^{-1})
!              https://www.ncnr.nist.gov/resources/n-lengths/
!     mu_in_incoh = \mu_{in incoh}  (cm^{-1})
!           = linear absorption coefficient by incoherent scattering of incoming neutron
!             (no dependent on the neutron velocity)
!     mu_out_absorp = \mu_{out absorp}  (cm^{-1})
!           = linear absorption coefficient by absorption of outgoing neutron
!     mu_out_incoh = \mu_{out incoh}  (cm^{-1})
!           = linear absorption coefficient by incoherent scattering of outgoing neutron
!
!  READ(*,*) X1,X2,NX
!     X1= lower limit of integration along x-axis (cm)
!     X2= upper limit of integration along x-axis (cm)
!     NX= division number of integration along x-axis ( >= 1)
!     x = X1 + (X2-X1)*(i-0.5)/NX  i=1,NX
!
!  READ(*,*) Y1,Y2,NY
!     Y1= lower limit of integration along y-axis (cm)
!     Y2= upper limit of integration along y-axis (cm)
!     NY= division number of integration along y-axis ( >= 1)
!     y = Y1 + (Y2-Y1)*(j-0.5)/NY  j=1,NY
!
!  READ(*,*) Z1,Z2,NZ
!     Z1= lower limit of integration along z-axis (cm)
!     Z2= upper limit of integration along z-axis (cm)
!     NZ= division number of integration along z-axis ( >= 1)
!     z = Z1 + (Z2-Z1)*(k-0.5)/NZ  k=1,NZ
!
!  READ(*,*) NPL,NCYL
!     NPL = number of planes (>= 0)
!     NCYL = number of cylinder (>= 0)
!
!  DO I=1,NPL
!    READ(*,*) PLA(1:3,I),PLB(I)
!  END DO
!
!  DO J=1,NCYL
!    READ(*,*) CYLC(1:3,J),CYLD(1:3,J),CYLR(J)
!  END DO
!
!  par file
!  Offset (deg)
!  Omega (deg), SPE file
!
! output
!     SPE file
!     S(Q,E)/ABSFAC  and  error_S(Q,E)/ABSFAC
!
!-----
MODULE mod_absorp
  IMPLICIT NONE
!
  private                                   ! set all names private
  public :: absorp_init, absorp_calc        ! public names
  public :: mu_in_absorp, mu_in_incoh, mu_out_absorp, mu_out_incoh
  public :: UT, VT, WT
!
  integer, parameter :: DP = kind(1.0D0)
!
  INTEGER(4), save :: NPL  ! number of planes
  INTEGER(4), save :: NCYL ! number of cylinder
  real(DP), save , ALLOCATABLE :: PLA(:,:),PLB(:)
  real(DP), save , ALLOCATABLE :: CYLC(:,:),CYLD(:,:),CYLR(:)
  real(DP), save :: X1,X2, Y1,Y2, Z1,Z2
  INTEGER(4), save :: NX, NY, NZ
  real(DP), save :: UT(3),VT(3),WT(3)
  real(DP), save :: mu_in_absorp,mu_in_incoh,mu_out_absorp,mu_out_incoh
!
contains
!
!-----
  subroutine absorp_init
!      read and initialize parameters of spectrometer
!
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    real(DP), parameter :: pi = 4.0d0*atan(1.0d0)
    integer(4) :: I,J
!
    READ(*,*) UT(1:3)
    WRITE(*,'(/,A,3G15.7)') '  direction parallel to u = UT(1:3) =',UT(1:3)
!
    READ(*,*) VT(1:3)
    WRITE(*,'(A,3G15.7)')   '  direction parallel to v = VT(1:3) =',VT(1:3)
!
    UT(1:3) = UT(1:3) / sqrt( dot_product( UT(1:3), UT(1:3) ) )    !  normalization
    VT(1:3) = VT(1:3) - UT(1:3) * dot_product( UT(1:3), VT(1:3) )  !  make VT(1:3) orthogonal to UT(1:3)
    VT(1:3) = VT(1:3) / sqrt( dot_product( VT(1:3), VT(1:3) ) )    !  normalization
    WT(1) = UT(2)*VT(3) - UT(3)*VT(2)                              !  vector product UT(1:3) x VT(1:3)
    WT(2) = UT(3)*VT(1) - UT(1)*VT(3)
    WT(3) = UT(1)*VT(2) - UT(2)*VT(1)
    WT(1:3) = WT(1:3) / sqrt( dot_product( WT(1:3), WT(1:3) ) )    !  normalization
    WRITE(*,'(A,3G15.7)')'  UT(1:3)=',UT(1:3)
    WRITE(*,'(A,3G15.7)')'  VT(1:3)=',VT(1:3)
    WRITE(*,'(A,3G15.7)')'  WT(1:3)=',WT(1:3)
!
    READ(*,*) mu_in_absorp, mu_in_incoh, mu_out_absorp, mu_out_incoh
    WRITE(*,'(/,a)') '  linear absorption coefficients for k_0 = 3.494158 (A^{-1})'
    WRITE(*,'(a)') '   mu_in_absorp   mu_in_incoh  (cm^{-1}) incoming'
    WRITE(*,'(4G15.7)') mu_in_absorp, mu_in_incoh
    WRITE(*,'(a)') '   mu_out_absorp   mu_out_incoh  (cm^{-1}) outgoing'
    WRITE(*,'(4G15.7)') mu_out_absorp, mu_out_incoh
!
    READ(*,*) X1,X2,NX
    READ(*,*) Y1,Y2,NY
    READ(*,*) Z1,Z2,NZ
    WRITE(*,'(/,a)') '   integration limits (cm) and divisions'
    WRITE(*,'(''  X1,X2,NX='',2G15.7,I5)')X1,X2,NX
    WRITE(*,'(''  Y1,Y2,NY='',2G15.7,I5)')Y1,Y2,NY
    WRITE(*,'(''  Z1,Z2,NZ='',2G15.7,I5)')Z1,Z2,NZ
!
    READ(*,*) NPL,NCYL
    WRITE(*,'(/,a,I5)') '  number of planes = NPL=',NPL
    if(NPL  >= 1)  ALLOCATE( PLA(3,NPL) , PLB(NPL) )
    if(NCYL >= 1)  ALLOCATE( CYLC(3,NCYL) , CYLD(3,NCYL) , CYLR(NCYL) )
!
    IF(NPL > 0) THEN
      DO I=1,NPL
        READ(*,*)PLA(1:3,I),PLB(I)
      end do
      WRITE(*,'(a)') '  PLA(1,I)*x+PLA(2,I)*y+PLA(3,I)*z >= PLB(I) (cm) (inside sample)'
      do I=1,NPL
        WRITE(*,'(a,i5)') '   I=',I
        WRITE(*,'(3(g15.7,a),g15.7)') PLA(1,I),'*x +' ,PLA(2,I),'*y +' ,PLA(3,I),'*z >=',PLB(I)
      end do
    ENDIF
!
    WRITE(*,'(/,a,I5)') '  number of cylinders = NCYL=',NCYL
    IF(NCYL > 0) THEN
      DO J=1,NCYL
        READ(*,*)CYLC(1:3,J),CYLD(1:3,J),CYLR(J)
      end do
      WRITE(*,'(A)') '  a point on the cylinder axis (cm) = CYLC(1:3,J)'
      WRITE(*,'(A)') '  a vector along the cylinder axis (no dim) = CYLD(1:3,J)'
      WRITE(*,'(a,a)') '  radius of the cylinder  (cm) =',' CYLR(J)'
      do J=1,NCYL
        WRITE(*,'(a,i5)') '   J=',J
        WRITE(*,'(a,3G15.7)') '  CYLC=',CYLC(1:3,J),'  CYLD=',CYLD(1:3,J),'  CYLR=',CYLR(J)
        CYLD(1:3,J)=CYLD(1:3,J)/sqrt( dot_product( CYLD(1:3,J), CYLD(1:3,J) ) )  !  nomalization
      end do
    ENDIF
!
  end subroutine absorp_init
!
!-----
  subroutine absorp_calc(UKI,UKF,mu_in,mu_out,ABSFAC)
!      calc absorption factor
!
    implicit none
    integer, parameter :: DP = kind(1.0D0)
    real(DP), intent(in) :: UKI(3),UKF(3),mu_in,mu_out
    real(DP), intent(out) :: ABSFAC
    real(DP) :: SUMN,SUMD,DELZ,X(3),Sout,F,XMC(3)
    real(DP) :: ALI,ALO,S,SS,SSS,CI,AI,BI,AO,BO
    integer(4) :: Ix,Jy,Kz,I,J,IFL_outside
!
    SUMN=0.0D0
    SUMD=0.0D0
    DELZ=(Z2-Z1)/DBLE(NZ)
    do Ix=1,NX
      X(1)=X1+(X2-X1)*(Ix-0.5D0)/DBLE(NX)
      do Jy=1,NY
        X(2)=Y1+(Y2-Y1)*(Jy-0.5D0)/DBLE(NY)
        do Kz=1,NZ
          X(3)=Z1+DELZ*(Kz-0.5D0)
!
          ALI=1.0d30  ;  ALO=1.0d30   ! call CALCSF(X,Sout,F)
          IFL_outside=0        ! inside sample 
!
          if(NPL > 0) then
            do I=1,NPL
              S   =  PLA(1,I)*X(1)  +PLA(2,I)*X(2)  +PLA(3,I)*X(3) - PLB(I)
              if(S < 0.0d0) then
                IFL_outside=1  ;  exit  ! outside sample
              endif
              SS  =  PLA(1,I)*UKI(1)+PLA(2,I)*UKI(2)+PLA(3,I)*UKI(3)
              SSS =-(PLA(1,I)*UKF(1)+PLA(2,I)*UKF(2)+PLA(3,I)*UKF(3))
              if(SS  > 0.0d0) ALI=MIN(ALI,S/SS)
              if(SSS > 0.0d0) ALO=MIN(ALO,S/SSS)
            end do
          end if
!
          if(IFL_outside == 0) then
          if(NCYL > 0) then
            do J=1,NCYL
              XMC(1)=X(1)-CYLC(1,J)
              XMC(2)=X(2)-CYLC(2,J)
              XMC(3)=X(3)-CYLC(3,J)
              SS=XMC(1)*CYLD(1,J)+XMC(2)*CYLD(2,J)+XMC(3)*CYLD(3,J)
              CI=XMC(1)*XMC(1)+XMC(2)*XMC(2)+XMC(3)*XMC(3) -SS*SS-CYLR(J)*CYLR(J)
              IF(CI > 0.0D0) then
                IFL_outside=1  ;  exit  ! outside sample
              endif
              S=UKI(1)*CYLD(1,J)+UKI(2)*CYLD(2,J)+UKI(3)*CYLD(3,J)
              AI=1.0d0-S*S
              BI=-(XMC(1)*UKI(1)+XMC(2)*UKI(2)+XMC(3)*UKI(3))+SS*S
              S=UKF(1)*CYLD(1,J)+UKF(2)*CYLD(2,J)+UKF(3)*CYLD(3,J)
              AO=1.0d0-S*S
              BO=XMC(1)*UKF(1)+XMC(2)*UKF(2)+XMC(3)*UKF(3)-SS*S
              IF(AI > 0.0d0) ALI=MIN(ALI,(-BI+SQRT(BI*BI-AI*CI))/AI)
              IF(AO > 0.0d0) ALO=MIN(ALO,(-BO+SQRT(BO*BO-AO*CI))/AO)
            end do
          endif
          endif
!
          if(IFL_outside == 0) then                      !  inside sample
            F=1.0d0  ;  Sout=EXP(-(mu_in * ALI + mu_out * ALO))
          else                                           !  outside sample
            F=0.0d0  ;  Sout=0.0d0
          endif
! return from CALCSF(X,Sout,F)
!
          SUMN=SUMN+Sout
          SUMD=SUMD+F
        end do
      end do
    end do
!
    ABSFAC=SUMN/SUMD
!
  end subroutine absorp_calc
!
END MODULE mod_absorp
!
!-----
! mod_in5.f90
!  2017/1/6  H Kadowaki
!
!-----
MODULE mod_in5
  IMPLICIT NONE
!
  private    ! set all names private
!     public names
  public :: in5_init, read_par_file, read_spe_file, absorp_set_drct, abs_corr_SQE_err_SQE, write_spe_file
  public :: dealloc_in5
  public :: Delta_Omega, nABS
  public :: non_zero_SQE_tab, ABSFAC_tab
  public :: UKI, mu_in
  public :: UKF_tab, mu_out_tab
!
  integer, parameter :: DP = kind(1.0D0)
!
  INTEGER(4), save :: ne  ! number of Ebin
  INTEGER(4), save :: nd  ! number of detector elements
  real(DP), save , ALLOCATABLE :: phi(:),azim(:),UKF_tab(:,:)
  real(DP), save , ALLOCATABLE :: Delta_E(:),mu_out_tab(:)
  real(DP), save :: UKI(3),mu_in,Delta_Omega
  real(DP), save :: wave_length,k_in
!
  real(DP), save, ALLOCATABLE :: SQE(:,:)
  real(DP), save, ALLOCATABLE :: err_SQE(:,:)
  integer(4), save :: nABS
  integer(4), save, ALLOCATABLE :: non_zero_SQE_tab(:,:)   ! ie,id
  real(DP), save, ALLOCATABLE :: ABSFAC_tab(:)
!
contains
!
!-----
  subroutine in5_init
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    real(DP), parameter :: pi = 4.0d0*atan(1.0d0)
    READ(*,*) wave_length
    WRITE(*,'(/,A,G15.7)') '  wave_length (A) =',wave_length
    k_in=(2.0d0*pi)/wave_length
    WRITE(*,'(A,G15.7)') '  k_in (A**^{-1}) =',k_in
  end subroutine in5_init
!
!-----
  subroutine read_par_file(nf_par)
!       read par file
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    integer(4), intent(in) :: nf_par
    integer(4) :: id
    real(DP) :: xx
!
    read(nf_par,*) nd
    write(*,'(A,I10)')'  number of detector elements = nd =',nd
    if(nd < 1) then
      write(*,*)'  invalid nd<1;  nd =',nd  ;  stop
    endif
    ALLOCATE( phi(nd)       )
    ALLOCATE( azim(nd)      )
    do id=1,nd
      read(nf_par,*) xx,phi(id),azim(id)
    end do
!
  end subroutine read_par_file
!
!-----
  subroutine absorp_set_drct
!        set neutron directions and absorption coefficients
!
    use mod_absorp, only : UT, VT, WT, mu_in_absorp, mu_in_incoh, mu_out_absorp, mu_out_incoh
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    real(DP), parameter :: pi = 4.0d0*atan(1.0d0)
    real(DP), parameter :: d2r = pi/180.0d0        !  degree to radian
!
    real(DP), parameter :: k_0=3.494158d0  ! k_0 = 3.494158 A^{-1} (v_0 = 2200 m/s)
!                                      see  https://www.ncnr.nist.gov/resources/n-lengths/list.html
    real(DP) :: k_out,ux_s(3),uy_s(3),uz_s(3),ang
    integer(4) :: ie,id
!      input:
!         Delta_Omega (degrees)
!         Delta_E(ie) = energy transfer  (meV)    ie-th Ebin
!      output:
!         UKI(1:3) = direction of incoming neutron
!         mu_in = linear absorption coefficients of incoming neutron
!         mu_out_tab(1:ne) = linear absorption coefficients of outgoing neutron
!         UKF_tab(1:3,1:nd) = direction of outgoing neutron
!
    ang = Delta_Omega*d2r
    UKI(1:3)= cos(ang) * UT(1:3) + sin(ang) * VT(1:3)
    mu_in = mu_in_absorp * k_0 / k_in + mu_in_incoh
    do ie=1,ne
      k_out = SQRT(k_in*k_in - Delta_E(ie)/2.072125d0)             ! k_in fixed
      mu_out_tab(ie) = mu_out_absorp * k_0 / k_out + mu_out_incoh  ! linear absorption coefficients of outgoing neutron
    end do
!
    ang = 0.5d0*pi + Delta_Omega*d2r
    ux_s(1:3)= cos(ang) * UT(1:3) + sin(ang) * VT(1:3)
    uy_s(1:3)= WT(1:3)
    uz_s(1:3)= UKI(1:3)
    do id=1,nd
      UKF_tab(1:3,id)= cos( azim(id)*d2r ) * sin( phi(id)*d2r ) * ux_s(1:3) &
         &           + sin( azim(id)*d2r ) * sin( phi(id)*d2r ) * uy_s(1:3) &
         &           +                       cos( phi(id)*d2r ) * uz_s(1:3)    ! direction of outging neutron
    end do
!
  end subroutine absorp_set_drct
!
!-----
  subroutine read_spe_file(nf_spe_in,nf_spe_out)
!       read spe file
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    integer(4), intent(in) :: nf_spe_in,nf_spe_out
    real(DP), ALLOCATABLE :: Ebin_boundary(:)
    real(DP) :: small
    character(128) :: oneline
    integer(4) :: i,ie,id,ie_begin,ie_end
!
    read(nf_spe_in,'(A128)')oneline
    write(nf_spe_out,'(A)')TRIM(oneline)
    read(oneline,*) i,ne
    write(*,'(/,A,I10)')'  number of detector elements = nd =',i
    write(*,'(A,I10)')'  number of energy bins = ne =',ne
      if(nd /= i) then
        write(*,*)'  invalid  i /= nd;  nd,i=',nd,i  ;  stop
      endif
    if(ne < 1) then
      write(*,*)'  invalid ne<1;  ne =',ne  ;  stop
    endif
    ALLOCATE( UKF_tab(3,nd)       )
    ALLOCATE( Delta_E(ne)         )
    ALLOCATE( mu_out_tab(ne)      )
    ALLOCATE( SQE(ne,nd)          )
    ALLOCATE( err_SQE(ne,nd)      )
    ALLOCATE( Ebin_boundary(ne+1) )
!
    read(nf_spe_in,'(A128)')oneline
    write(nf_spe_out,'(A)')TRIM(oneline)
    if( oneline(1:12) /= '### Phi Grid' ) then
      write(*,*) " oneline(1:12) /= '### Phi Grid' "  ;  stop
    endif
    do id = 1, nd+1, 8
      read(nf_spe_in,'(A128)')oneline
      write(nf_spe_out,'(A)')TRIM(oneline)
    end do
!
    read(nf_spe_in,'(A128)')oneline
    write(nf_spe_out,'(A)')TRIM(oneline)
    if( oneline(1:15) /= '### Energy Grid' ) then
      write(*,*) " oneline(1:15) /= '### Energy Grid' "  ;  stop
    endif
!
    do ie = 1, ne+1, 8
      read(nf_spe_in,'(A128)')oneline
      write(nf_spe_out,'(A)')TRIM(oneline)
      ie_begin = ie
      ie_end = MIN( ie+7 , ne+1 )
      read(oneline,*) Ebin_boundary( ie_begin : ie_end )
    end do
    do ie=1,ne
      Delta_E(ie)=0.5d0*( Ebin_boundary(ie) + Ebin_boundary(ie+1) )
    end do
    DEALLOCATE( Ebin_boundary )
!
    do id=1,nd
      read(nf_spe_in,'(A128)')oneline
      if( oneline(1:12) /= '### S(Phi,w)' ) then
        write(*,*) " oneline(1:12) /= '### S(Phi,w)' "  ;  stop
      endif
      read(nf_spe_in,*) SQE( 1:ne , id )
!
      read(nf_spe_in,'(A128)')oneline
      if( oneline(1:10) /= '### Errors' ) then
        write(*,*) " oneline(1:10) /= '### Errors' "  ;  stop
      endif
      read(nf_spe_in,*) err_SQE( 1:ne , id )
    end do
!
    small = 1.0d-15
    i=0
    do id=1,nd
    do ie=1,ne
      if( ABS( SQE(ie,id) ) + ABS( err_SQE(ie,id) ) < small ) then
        cycle
      else
        i=i+1
      endif
    end do
    end do
    nABS=i
    write(*,'(A,I10)')'  number of nonzero S(Q,E) or error_S(Q,E) = nABS =',nABS
!
    ALLOCATE( non_zero_SQE_tab(2,nABS) )
    ALLOCATE( ABSFAC_tab(nABS)         )
    i=0
    do id=1,nd
    do ie=1,ne
      if( ABS( SQE(ie,id) ) + ABS( err_SQE(ie,id) ) < small ) then
        cycle
      else
        i=i+1
        non_zero_SQE_tab(1,i)=ie
        non_zero_SQE_tab(2,i)=id
      endif
    end do
    end do
!
  end subroutine read_spe_file
!-----
  subroutine abs_corr_SQE_err_SQE
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    real(DP) :: ABS_corr,ABS_corr_max,ABS_corr_min
    integer(4) :: i,ie,id
    ABS_corr_max = 1.0d-30
    ABS_corr_min = 1.0d+30
    do i = 1, nABS
      ie = non_zero_SQE_tab(1,i)
      id = non_zero_SQE_tab(2,i)
      ABS_corr = 1.0d0 / MAX( ABSFAC_tab(i), 1.0d-15 )
      ABS_corr_max = MAX( ABS_corr , ABS_corr_max )
      ABS_corr_min = MIN( ABS_corr , ABS_corr_min )
      SQE(ie,id) = SQE(ie,id) * ABS_corr
      err_SQE(ie,id) = err_SQE(ie,id) * ABS_corr
    end do
    write(*,'(/,A,2G15.7)')'  ABS_corr_min, ABS_corr_max =',ABS_corr_min,ABS_corr_max
  end subroutine abs_corr_SQE_err_SQE
!
!-----
  subroutine write_spe_file(nf_spe_out)
!       write spe file
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    integer(4), intent(in) :: nf_spe_out
    integer(4) :: id,ie_begin,ie_end
    character(128) :: oneline
!
    do id=1,nd
      write(nf_spe_out,'(A)')'### S(Phi,w)'
      do ie_begin = 1, ne, 8
        ie_end = MIN( ie_begin + 7 , ne )
        write(oneline,'(8ES11.3E2)') SQE( ie_begin : ie_end , id )
        call E2e(oneline)
        write(nf_spe_out,'(A)') TRIM(oneline)
      end do
      write(nf_spe_out,'(A)')'### Errors'
      do ie_begin = 1, ne, 8
        ie_end = MIN( ie_begin + 7 , ne )
        write(oneline,'(8ES11.3E2)') err_SQE( ie_begin : ie_end , id )
        call E2e(oneline)
        write(nf_spe_out,'(A)') TRIM(oneline)
      end do
    end do
!
  end subroutine write_spe_file
!
!-----
  subroutine E2e(oneline)
    IMPLICIT NONE
    character(128), intent(inout) :: oneline
    integer(4) :: i,j
    do j=1,128
      i=INDEX(oneline,'E')
      if(i == 0) exit
      oneline(i:i)='e'
    end do
  end subroutine E2e
!
!------
  subroutine dealloc_in5
    IMPLICIT NONE
    DEALLOCATE( UKF_tab          )
    DEALLOCATE( Delta_E          )
    DEALLOCATE( mu_out_tab       )
    DEALLOCATE( SQE              )
    DEALLOCATE( err_SQE          )
    DEALLOCATE( non_zero_SQE_tab )
    DEALLOCATE( ABSFAC_tab       )
  end subroutine dealloc_in5
!
END MODULE mod_in5
!
