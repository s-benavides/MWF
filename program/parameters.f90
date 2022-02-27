#include "../parallel.h"
 module parameters
 use mpif
!***************************************************************************
   implicit none
   save

   double precision            ::  d_Re         != 180d0!!67.8d0!875d0
   
   !NUMBER OF MODES TO USE IE HIGHEST WAVENUMBER + 1
   integer,          parameter :: i_MM          = 512!4096!512!2048!4096 !Streamwise
   integer,          parameter :: i_NN          = 256!1024!512!2048!4096 !Spanwise
   integer,          parameter :: i_K0          = 4
   integer                     :: i_kICx        ! Parameter for random initial  condition
   integer                     :: i_KICz        ! Parameter for random initial condition
   
   double precision, parameter :: d_PI          = 3.1415926535897931d0
   double precision            :: d_Lx          != 180d0!4096d0!960d0 !Streamwise
   double precision            :: d_Lz          != 80d0!4096d0!1024d0!960d0! Spanwise
   double precision            :: d_E0          ! Initial condition KE
   double precision            :: d_alpha       != 2d0*d_PI/d_Lx!0.5d0
   double precision            :: d_gamma       != 2d0*d_PI/d_Lz!0.5d0

   logical,          parameter :: s_reflect     = .false.!.TRUE.!.FALSE. 
   logical,          parameter :: s_uvreflect   = .FALSE.
   logical,          parameter :: s_wreflect    = .false.

   integer                     :: i_maxtstep    != 1d8
   integer                     :: i_save_rate1  != 9000!5000!12500!1d8!50000!100000!5000
   integer                     :: i_save_rate2  != 230!25!25!50 
   double precision, parameter :: d_maxt        = -1d0
   double precision            :: d_cpuhours    != 11.8d0!1.7d0!19.6d0!0.98d0!1d99 !90d0
   double precision            :: d_time        != -1d0 !0d0  ! if d_time < 0, then it reads it from the statefile.
   integer                     :: i_tstep       != 1d8
   double precision, parameter :: d_thdeg       = 0d0!24d0
   double precision            :: d_dt          != 0.00025d0!0.01d0!0.02d0
   double precision, parameter :: d_minE        = 1d-10

   double precision            :: d_HYPO        != 0d0!1d-6 !1d-01
   integer                     :: i_PHYPO       != 2
   double precision            :: d_drag        != 1d-2 !1d-01
   double precision            :: d_vdrag       != 0d0
   logical,          parameter :: s_dragall     = .true.!.false.

   logical, parameter          :: s_HIS         = .FALSE. !HISTORY DATA
   integer, parameter          :: i_H           = 32
   logical, parameter          :: s_TUR         = .TRUE.! TURBULENT DATA
   integer, parameter          :: i_DM          = 6!60!!0!12!60!48!6!6!3!96 !Single point in M
   integer, parameter          :: i_DN          = 6!6!12!6!3!9  !equal to dz of 1.875
   integer, parameter          :: i_MT          = 10!Number of times to average over
   integer, parameter          :: i_WT          = 10!20  ! Gaps of MT to make between output
   !---------------------------------------------------------------------------
   !  Fixed parameters
   !---------------------------------------------------------------------------
   integer,          parameter :: i_N           = 2*(i_NN-1)
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

     ! Load parameters from file 'parameter.inp'
     NAMELIST / parameters / d_Re, d_Lx, d_Lz, d_E0,i_kICx,i_kICz, i_save_rate1, i_save_rate2, i_maxtstep, d_cpuhours, d_time,i_tstep, d_dt, d_HYPO, i_PHYPO, d_drag, d_vdrag
       if (mpi_rnk==0) then
          open(1,file='parameter.inp',status='unknown',form='formatted')
          read(1,NML=parameters)
          close(1)
      endif
#ifdef _MPI
      call mpi_bcast(d_Re,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_Lx,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_Lz,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_E0,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(i_kICx,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(i_kICz,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(i_save_rate1,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(i_save_rate2,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(i_maxtstep,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_cpuhours,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_time,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(i_tstep,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_dt,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_HYPO,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(i_PHYPO,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_drag,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
      call mpi_bcast(d_vdrag,1,mpi_double_precision,0,mpi_comm_world,mpi_er)
#endif

      d_alpha = 2d0*d_PI/d_Lx
      d_gamma = 2d0*d_PI/d_Lz 

      if (d_time.ge.0d0) then
         tim_t=d_time
         tim_step= i_tstep
      else
         tim_t=0d0
         tim_step=0
      end if

   end subroutine par_precompute
 

!***************************************************************************
 end module parameters
!***************************************************************************
