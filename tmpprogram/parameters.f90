#include "../parallel.h"
 module parameters
!***************************************************************************
   implicit none
   save

   double precision            ::  d_Re         = 100d0
   
   !NUMBER OF MODES TO USE IE HIGHEST WAVENUMBER + 1
   integer,          parameter :: i_MM	        = 32!32!32!!192!64
   integer,          parameter :: i_NN          = 128!256!32!32!64!128
   
   integer,          parameter :: i_K0          = 4
   
   double precision, parameter :: d_PI          = 3.1415926535897931d0
   double precision, parameter :: d_Lx          = 10d0!480d0!120d0!!2d0*d_PI/d_alpha
   double precision, parameter :: d_Lz          = 40d0!240d0!240d0!240d0!120d0!!2d0*d_PI/d_gamma
   double precision, parameter :: d_alpha       = 2d0*d_PI/d_Lx!0.5d0
   double precision, parameter :: d_gamma       = 2d0*d_PI/d_Lz!0.5d0

   logical,          parameter :: s_reflect     = .false.!.TRUE.!.FALSE. 
   logical,          parameter :: s_uvreflect   = .FALSE.
   logical,          parameter :: s_wreflect    = .false.

   integer, parameter          :: i_maxtstep    = 50!1000000!1000
   integer, parameter          :: i_save_rate1  = 5000!5000
   integer, parameter          :: i_save_rate2  = 5!50!50 
   double precision, parameter :: d_maxt        = -1d0
   double precision, parameter :: d_cpuhours    = 9.75d0!1d99 !90d0
   double precision, parameter :: d_time        = -1d0 !0d0
   double precision, parameter :: d_thdeg       = 24d0
   double precision, parameter :: d_dt          = 0.01d0
   double precision, parameter :: d_minE        = 1d-8
   
   logical, parameter          :: s_HIS         = .FALSE. !HISTORY DATA
   integer, parameter          :: i_H           = 32
   logical, parameter          :: s_TUR         = .TRUE.!FALSE. !TURBULENT DATA
   integer, parameter          :: i_DM          = 6
   integer, parameter          :: i_DN          = 6
   integer, parameter          :: i_MT          = 10
   !---------------------------------------------------------------------------
   !  Fixed parameters
   !---------------------------------------------------------------------------
   integer,	     parameter :: i_N           = 2*(i_NN-1)
   integer,          parameter :: i_K           = 2*i_K0-1
   integer,          parameter :: i_KK          = 3*i_K0-1
   integer,          parameter :: i_M           = 2*(i_MM-1)
   integer,          parameter :: i_Ny          = 15
   double precision, parameter :: d_beta        = d_PI/2d0
   double precision, parameter :: d_theta       = d_thdeg/180d0*d_PI

   integer,          parameter :: i_3M          = 3*i_MM
   integer,          parameter :: i_3N          = 3*i_NN
   integer,          parameter :: i_SM          = i_3M/i_DM
   integer,          parameter :: i_SN          = i_3N/i_DN
   
   integer,          parameter :: i_Mp          = (_Np + i_3M - 1)/_Np
   integer,          parameter :: i_Np          = (_Np + i_NN- 1)/_Np
   integer,          parameter :: i_SMp         = (_Np + i_SM - 1)/_Np

   integer,          parameter :: i_NN1  = i_NN-1
   integer,          parameter :: i_N1  = i_N-1
   integer,          parameter :: i_M1  = i_M-1
   integer,          parameter :: i_MM1  = i_MM-1
   integer,          parameter :: i_KK1  = i_KK-1
   integer,          parameter :: i_K1  = i_K0-1

   double precision :: tim_t
   integer          :: tim_step
!---------------------------------------------------------------------------
!  check params are ok
!---------------------------------------------------------------------------

contains

   subroutine par_precompute()
      integer :: itmp
      if (d_time > 0d0) then
         tim_t=d_time
      else
         tim_t=0d0
      end if
      tim_step=0
      itmp=mod(i_M,_Np)
!      if(itmp /= 0) stop 'mpi_precompute: incorrect num procs, M'
      itmp=mod(i_NN,_Np)
!      if(itmp /= 0) stop 'mpi_precompute: incorrect num procs, N'
   end subroutine par_precompute
 

!***************************************************************************
 end module parameters
!***************************************************************************
