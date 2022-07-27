!**************************************************************************
!  IN/OUT 
!
!**************************************************************************
#include "../parallel.h"
 module io
!**************************************************************************
   use velocity
   use mpif
   use netcdf
   use turb
   implicit none
   save

   character(200)        :: io_statefile
   integer               :: io_save1
   integer,     private  :: io_KE, io_HI, io_uq
   type (mpt),  private  :: c1
   type (spec),  private  :: s1
   type (spec_xavg_odd), private :: x1o
   type (spec_xavg_even), private :: x1e
   character(200) :: oloc='RUNNING'
 contains
 
!--------------------------------------------------------------------------
!  initialiser fn
!--------------------------------------------------------------------------
   subroutine io_precompute()
     character(4) :: cnum
     type (phys)   :: vp
     logical       :: exist
     double precision :: E
     _loop_kmn_vars

     ! Check
     ! Namelist for input files
     NAMELIST / status / io_save1, turb_save2, i_restress_save
       if (mpi_rnk==0) then
          open(1,file='parameter.inp',status='unknown',form='formatted')
          read(1,NML=status)
          close(1)
      endif
#ifdef _MPI
      call mpi_bcast(io_save1,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(turb_save2,1,mpi_integer,0,mpi_comm_world,mpi_er)
      call mpi_bcast(i_restress_save,1,mpi_integer,0,mpi_comm_world,mpi_er)
#endif
      i_restress_save = i_restress_save + 1 ! So that it won't overwrite
      io_KE    = 20
      io_hi    = 0      
      io_uq    = 0
      if (s_HIS) &
           io_hi    = 30
      if (s_uq)  &
           io_uq = 40
           if (d_theta > 0) print*, "WARNING, you are outputting (u,q), but this &
           was only intended to be used in the theta=0 case. It is not ouputting &
            the correct quantity now."

       !!! INITIAL CONDITIONS

       if (io_save1==0) then 
        ! Make initial condition and save it, then load it as state.cdf.in 
        call var_randphasempt(vel_c)
        !Calculate energy of this IC
        call vel_mpt2phys(vel_c,vp)
        call vel_energy(vp,E)
        ! Now renormalize and multiply by energy we want.
        _loop_kmn_begin
            if ((n+var_N%pH0)==0) cycle
            if (m==0) cycle
            vel_c%Re(k,m,n) = dsqrt(d_E0)*vel_c%Re(k,m,n)/dsqrt(E)
            vel_c%Im(k,m,n) = dsqrt(d_E0)*vel_c%Im(k,m,n)/dsqrt(E)
        _loop_kmn_end
        call vel_mpt2phys(vel_c,vp)
        call vel_energy(vp,E)

        call io_save_IC()
        
#ifdef _MPI
        call mpi_barrier(mpi_comm_world, mpi_er)
#endif

        io_statefile = 'state.cdf.in'
      else
        write(cnum,'(I4.4)') io_save1
        if (mpi_rnk==0)then
            print*, 'starting from state'// cnum
        endif
        io_statefile = 'state'//cnum//'.cdf.dat'
      endif

      !!! SPEC_IN 
      inquire(file='spec.cdf.in',exist=exist)
      if (exist) call io_load_spec('spec.cdf.in',spec_in)
   end subroutine io_precompute 
 

!--------------------------------------------------------------------------
!  Open files written to every ... steps runtime
!--------------------------------------------------------------------------

   subroutine io_openfiles()
      logical :: exist
      character(10), save :: s = 'replace', a = 'sequential', p='append'
      if(mpi_rnk/=0) return
      ! Check for currently existing vel_energy.dat
      inquire(file='vel_energy.dat',exist=exist)
      if (exist) then
          s = 'old'
          if(io_KE/=0)  open(io_KE,status=s,position=p, file='vel_energy.dat')
      else 
          s = 'new'
          if(io_KE/=0)  open(io_KE,status=s,position=p, file='vel_energy.dat')
      endif
      ! Check for currently existing vel_history.dat
      inquire(file='vel_history.dat',exist=exist)
      if (exist) then
          s = 'old'
          if(io_hi/=0)  open(io_hi,status=s,position=p, file='vel_history.dat')
      else 
          s = 'new'
          if(io_hi/=0)  open(io_hi,status=s,position=p, file='vel_history.dat')
      endif
      ! Check for currently existing uq.dat
      inquire(file='uq.dat',exist=exist)
      if (exist) then
          s = 'old'
          if(io_uq/=0)  open(io_uq,status=s,position=p, file='uq.dat')
      else 
          s = 'new'
          if(io_uq/=0)  open(io_uq,status=s,position=p, file='uq.dat')
      endif
      s = 'old'
!      a = 'stream'
   end subroutine io_openfiles


!--------------------------------------------------------------------------
!  Close files written to during runtime
!--------------------------------------------------------------------------
   subroutine io_closefiles()
      if(mpi_rnk/=0) return
      if(io_KE/=0) close(io_KE)
      if(io_hi/=0) close(io_hi)
      if(io_uq/=0) close(io_uq)
   end subroutine io_closefiles


   !-------------------------------------------------------------------------
   subroutine io_clk_time(t)
     real, intent(out) :: t
#ifdef _CPUTIME
     call cpu_time(t)
#else
     integer, save :: ct,ctrt,ctmx,ct_=0,wrap=0
     call system_clock(ct,ctrt,ctmx)
     if(ct<ct_) wrap = wrap + 1
     ct_ = ct
     t = (real(ct)+real(ctmx)*wrap)/real(ctrt)
#endif
   end subroutine io_clk_time
    
!-------------------------------------------------------------------------


!--------------------------------------------------------------------------
!  Write to files
!--------------------------------------------------------------------------
   subroutine io_write2files()
      type (phys) :: vp
      if(modulo(tim_step,i_save_rate1)==0) then
!         if (tim_step > 0) then
            call io_save_state()
            call io_save_Uspec()
!         end if
         io_save1 = io_save1+1
      endif

      if(modulo(tim_step,i_save_rate2)==0) then
         call vel_mpt2phys(vel_c,vp)
         if(io_KE/=0) call io_write_energy(vp)
         if(io_hi/=0) call io_write_history(vp)
         if(io_uq/=0) call io_write_uq(vel_c)
         if(s_tur) then
            call turb_sample(vp)
            if (modulo(tim_step,i_WT*i_MT*i_save_rate2)==0) &
                 call turb_save_state()
         end if
      end if

     ! Reynolds stresses
     if ((d_avg_time > d_avg_window).and.((s_restress_xavg).or.(s_restress_2d).or.(s_restress_filt))) then
       ! Exports Reynolds averages and stresses
       if (s_restress_xavg) call io_save_xavg()
       if (s_restress_2d) call io_save_restress_2d()
       if (s_restress_filt) call io_save_restress_filt()
       
       ! Updates counters
       i_restress_save = i_restress_save + 1 
       i_count = 0
       d_avg_time = 0.0

       ! Resets average values
       if (s_restress_xavg) call vel_restress_reset()
       if (s_restress_2d) call vel_restress_2d_reset()
     end if

! Every 50 i_save_rate2 reopen vel_energy to
! overcome buffering and force output
      if(modulo(tim_step,i_save_rate2*50)==0) then
         call io_closefiles()
         call io_openfiles()
      end if

   end subroutine io_write2files


!--------------------------------------------------------------------------
!  Load state - start from previous solution
!--------------------------------------------------------------------------
   subroutine io_load_state()
      integer :: e, f, i, rd
      double precision :: d

      if(mpi_rnk==0)      print*, "Loading: ",io_statefile
      e=nf90_open(io_statefile, nf90_nowrite, f)      
      if(e/=nf90_noerr) then
         if(mpi_rnk==0) open(99, file='PRECOMPUTING')
         if(mpi_rnk==0) close(99, status='delete')
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'state file not found!'
      end if

      e=nf90_get_att(f,nf90_global,'t', d)
      if(d_time<0d0) tim_t = d
      if(mpi_rnk==0 .and. dabs(tim_t-d)>1d-8)  &
         print*,' t    :',d,' --> ', tim_t

      e=nf90_get_att(f,nf90_global,'tstep', d)
      if(d_time<0d0) tim_step = d
      if(mpi_rnk==0 .and. dabs(tim_step-d)>1d-8)  &
         print*,' tstep    :',d,' --> ', tim_step

      e=nf90_get_att(f,nf90_global,'Re', d)
      if(mpi_rnk==0 .and. dabs(d_Re-d)>1d-8)  &
         print*,' Re   :',d,' --> ', d_Re
      e=nf90_get_att(f,nf90_global,'alpha', d)
      if(mpi_rnk==0 .and. dabs(d_alpha-d)>1d-8)  &
         print*,' alpha:',d,' --> ', d_alpha
     e=nf90_get_att(f,nf90_global,'gamma', d)
      if(mpi_rnk==0 .and. dabs(d_gamma-d)>1d-8)  &
         print*,' gamma:',d,' --> ', d_gamma

      call io_load_mpt(f,'mpt',vel_c)

      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

   end subroutine io_load_state

!--------------------------------------------------------------------------
!  Load spec file
!--------------------------------------------------------------------------
   subroutine io_load_spec(io_specfile,u)
      integer :: e, f, i
      character(*),  intent(in) :: io_specfile
      type (spec),   intent(out) :: u
      integer :: K__, M__, N__
      logical :: error=.false.

      if(mpi_rnk==0)      print*, "Loading: ",io_specfile
      e=nf90_open(io_specfile, nf90_nowrite, f)      
      if(e/=nf90_noerr) then
         if(mpi_rnk==0) open(99, file='PRECOMPUTING')
         if(mpi_rnk==0) close(99, status='delete')
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'spec file not found!'
      end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      e=nf90_inq_varid(f,'spec', i)
      if(e/=nf90_noerr)  print*, 'Field spec not found!'
      if(e/=nf90_noerr)  stop 'io_load_spec'
      e=nf90_get_att(f,i, 'KK',  K__)
      e=nf90_get_att(f,i, 'MM',  M__)
      e=nf90_get_att(f,i, 'NN',  N__)
      if(K__ /=i_KK)  error=.true.
      if(M__ /=i_MM) error=.true.
      if(N__ /=i_NN) error=.true.

      if(mpi_rnk==0) then
         if(K__ /=i_KK)  print*,'spec  KK :', K__, ' --> ',i_KK 
         if(M__ /=i_MM)  print*,'spec  MM :', M__, ' --> ',i_MM 
         if(N__ /=i_NN)  print*,'spec  NN :', N__,' --> ',i_NN 
      end if
      if (error) then
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'wrong state size'
      end if
      
      u%Re = 0d0
      u%Im = 0d0
      
      e=nf90_get_var(f,i, u%Re(:,:,:),start=(/1,1,var_N%pH0+1,1/))
      e=nf90_get_var(f,i, u%Im(:,:,:),start=(/1,1,var_N%pH0+1,2/))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

   end subroutine io_load_spec

!--------------------------------------------------------------------------
!  Load coll variable
!--------------------------------------------------------------------------
   subroutine io_load_mpt(f,nm, a)
      integer,          intent(in)  :: f
      character(*),     intent(in)  :: nm
      type (mpt),      intent(out) :: a
      integer :: N1,M1, mm_,mm,nn,e,i
      integer :: K__, M__, N__,K1,K2
      integer :: m,n
      logical :: error=.false.
          
      e=nf90_inq_varid(f,nm, i)
      if(e/=nf90_noerr)  print*, 'Field '//nm//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_get_att(f,i, 'K',  K__)
      e=nf90_get_att(f,i, 'MM',  M__)
      e=nf90_get_att(f,i, 'NN',  N__)
      if(K__ /=i_K)  error=.true.
      if(M__ /=i_MM) error=.true.
      if(N__ /=i_NN) error=.true.

      if(mpi_rnk==0) then
         if(K__ /=i_K)  print*, nm, ' K :', K__, ' --> ',i_K 
         if(M__ /=i_MM)  print*, nm, ' MM :', M__, ' --> ',i_MM 
         if(N__ /=i_NN)  print*, nm, ' NN :', N__,' --> ',i_NN 
      end if
      if (error) then
#ifdef _MPI
         call mpi_barrier(mpi_comm_world, mpi_er)
         call mpi_finalize(mpi_er)
#endif
         stop 'wrong state size'
      end if
      
      a%Re = 0d0
      a%Im = 0d0
      
      e=nf90_get_var(f,i, a%Re(:,:,:),start=(/1,1,var_N%pH0+1,1/))
      e=nf90_get_var(f,i, a%Im(:,:,:),start=(/1,1,var_N%pH0+1,2/))
      
    end subroutine io_load_mpt

   subroutine io_load_mpt_old(f,nm, a)
      integer,          intent(in)  :: f
      character(*),     intent(in)  :: nm
      type (mpt),      intent(out) :: a
      integer :: N1,M1, mm_,mm,nn,e,i
      integer :: K__, M__, N__,K1,K2
      integer :: m,n
          
      e=nf90_inq_varid(f,nm, i)
      if(e/=nf90_noerr)  print*, 'Field '//nm//' not found!'
      if(e/=nf90_noerr)  stop 'io_load_coll'
      e=nf90_get_att(f,i, 'K',  K__)
      e=nf90_get_att(f,i, 'MM',  M__)
      e=nf90_get_att(f,i, 'NN',  N__)
      if(mpi_rnk==0) then
         if(K__ /=i_K)  print*, nm, ' K :', K__, ' --> ',i_K 
         if(M__ /=i_MM)  print*, nm, ' MM :', M__, ' --> ',i_MM 
         if(N__ /=i_NN)  print*, nm, ' NN :', N__,' --> ',i_NN 
      end if

      a%Re = 0d0
      a%Im = 0d0

      N1 = min(N__,i_NN)-1
      M1 = min(M__,i_MM)-1
      K1 = min(i_K0,(K__+1)/2)
      K2 = (K__+1)/2+1
      do n = 0, N__-1 
         if(n<var_N%pH0 .or. n>var_N%pH0+var_N%pH1)  cycle
         nn = n - var_N%pH0
!         print*,n,mpi_rnk,nn
         do m = -M__,M__
            if (abs(m) > i_MM1) cycle
            mm_ = modulo(m,2*(M__-1))
            mm  = modulo(m,i_M)
            if (i_K == K__) then
               e=nf90_get_var(f,i, a%Re(1:i_K,mm,nn),start=(/1,mm_+1,n+1,1/))
               e=nf90_get_var(f,i, a%Im(1:i_K,mm,nn),start=(/1,mm_+1,n+1,2/))
            else 
               e=nf90_get_var(f,i, a%Re(1:K1,mm,nn)            ,start=(/1,mm_+1,n+1,1/))
               e=nf90_get_var(f,i, a%Re(i_K0+1:i_K0+K1-1,mm,nn),start=(/K2,mm_+1,n+1,1/))
               e=nf90_get_var(f,i, a%Im(1:K1,mm,nn)            ,start=(/1,mm_+1,n+1,2/))
               e=nf90_get_var(f,i, a%Im(i_K0+1:i_K0+K1-1,mm,nn),start=(/K2,mm_+1,n+1,2/))
            end if
         end do
      end do
    end subroutine io_load_mpt_old

!--------------------------------------------------------------------------
!  Save state
!--------------------------------------------------------------------------
   subroutine io_save_IC()
      character(4) :: cnum
      integer :: e, f
      integer :: md,nd,kd, ReImd, dims(4)
      integer :: ss

      write(cnum,'(I4.4)') io_save1

      if(mpi_rnk==0) then
         print*, ' saving initial condition'
         e=nf90_create('state.cdf.in', nf90_clobber, f)
         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'tstep', tim_step)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         e=nf90_put_att(f, nf90_global, 'gamma', d_gamma)

         e=nf90_def_dim(f, 'N', i_NN, nd)
         e=nf90_def_dim(f, 'M', i_M, md)
         e=nf90_def_dim(f, 'K', i_K, kd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         dims = (/kd,md,nd,ReImd/)
         call io_define_mpt(f, 'mpt', dims, ss)

         e=nf90_enddef(f)
      end if

      call io_save_mpt(f,ss, vel_c)

      if(mpi_rnk==0)  &
         e=nf90_close(f)

   end subroutine io_save_IC

!--------------------------------------------------------------------------
!  Save state
!--------------------------------------------------------------------------
   subroutine io_save_state()
      character(4) :: cnum
      integer :: e, f
      integer :: md,nd,kd, ReImd, dims(4)
      integer :: ss 

      write(cnum,'(I4.4)') io_save1

      if(mpi_rnk==0) then
         print*, ' saving state'//cnum//'  t=', tim_t
         e=nf90_create('state'//cnum//'.cdf.dat', nf90_clobber, f)
         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'tstep', tim_step)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         e=nf90_put_att(f, nf90_global, 'gamma', d_gamma)

         e=nf90_def_dim(f, 'N', i_NN, nd)
         e=nf90_def_dim(f, 'M', i_M, md)
         e=nf90_def_dim(f, 'K', i_K, kd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         dims = (/kd,md,nd,ReImd/)
         call io_define_mpt(f, 'mpt', dims, ss)

         e=nf90_enddef(f)
      end if

      call io_save_mpt(f,ss, vel_c)

      if(mpi_rnk==0)  &
         e=nf90_close(f)

   end subroutine io_save_state

!--------------------------------------------------------------------------
!  Save Reynold averages and stresses
!--------------------------------------------------------------------------
   subroutine io_save_xavg()
      character(4) :: cnum
      integer :: e, f
      integer :: nd,kd, ReImd, dims(3)
      integer :: ssu,ssv,ssw,ssuu,ssuv,ssuw,ssww,sswv,ssvv 

      write(cnum,'(I4.4)') i_restress_save

      if(mpi_rnk==0) then
         print*, ' saving Reynolds averages and stresses (even)'//cnum//'  t=', tim_t
         e=nf90_create('restress_even'//cnum//'.cdf.dat', nf90_clobber, f)
         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'tstep', tim_step)
         e=nf90_put_att(f, nf90_global, 'avg_window', d_avg_window)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         e=nf90_put_att(f, nf90_global, 'gamma', d_gamma)

         e=nf90_def_dim(f, 'N', i_NN, nd)
         e=nf90_def_dim(f, 'K', i_K0, kd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)
         dims = (/kd,nd,ReImd/)

         call io_define_restress(f, 'umean', dims, ssu,i_K0)
         call io_define_restress(f, 'wmean', dims, ssw,i_K0)
         call io_define_restress(f, 'uumean', dims, ssuu,i_K0)
         call io_define_restress(f, 'uwmean', dims, ssuw,i_K0)
         call io_define_restress(f, 'wwmean', dims, ssww,i_K0)
         call io_define_restress(f, 'vvmean', dims, ssvv,i_K0)
         
         e=nf90_enddef(f)
      end if

      call io_save_restress_even(f,ssu,umean)
      call io_save_restress_even(f,ssw,wmean)
      call io_save_restress_even(f,ssuu,uumean)
      call io_save_restress_even(f,ssuw,uwmean)
      call io_save_restress_even(f,ssww,wwmean)
      call io_save_restress_even(f,ssvv,vvmean)

      if(mpi_rnk==0)  &
         e=nf90_close(f)

      if(mpi_rnk==0) then
         print*, ' saving Reynolds averages and stresses (odd)'//cnum//'  t=', tim_t
         e=nf90_create('restress_odd'//cnum//'.cdf.dat', nf90_clobber, f)
         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'tstep', tim_step)
         e=nf90_put_att(f, nf90_global, 'avg_window', d_avg_window)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         e=nf90_put_att(f, nf90_global, 'gamma', d_gamma)

         e=nf90_def_dim(f, 'N', i_NN, nd)
         e=nf90_def_dim(f, 'K', i_K1, kd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)
         dims = (/kd,nd,ReImd/)

         ! Odd parity functions
         call io_define_restress(f, 'vmean', dims, ssv,i_K1)
         call io_define_restress(f, 'uvmean', dims, ssuv,i_K1)
         call io_define_restress(f, 'wvmean', dims, sswv,i_K1)
         
         e=nf90_enddef(f)
      end if

      call io_save_restress_odd(f,ssv,vmean)
      call io_save_restress_odd(f,sswv,wvmean)
      call io_save_restress_odd(f,ssuv,uvmean)
      
      if(mpi_rnk==0)  &
         e=nf90_close(f)

   end subroutine io_save_xavg

   subroutine io_define_restress(f,nm,dims, id,KK)
      integer,      intent(in) :: f, dims(3)
      character(*), intent(in) :: nm
      integer,      intent(in) :: KK
      integer, intent(out) :: id
      integer :: e
      e=nf90_def_var(f, nm, nf90_double, dims, id)
      e=nf90_put_att(f, id,  'K', KK)      
      e=nf90_put_att(f, id,  'N', i_NN)
   end subroutine io_define_restress

   subroutine io_save_restress_even(f,id,a)
      integer,     intent(in) :: f, id
      type (spec_xavg_even), intent(in) :: a
      integer :: e
      
#ifndef _MPI
      e=nf90_put_var(f,id,a%Re(1:i_K0,0:i_NN1), start=(/1,1,1/))
      e=nf90_put_var(f,id,a%Im(1:i_K0,0:i_NN1), start=(/1,1,2/))

#else
      integer :: r, pN0,pN1

      if(mpi_rnk==0) then
         e=nf90_put_var(f,id,a%Re(1:i_K0,0:var_N%pH1), start=(/1,1,1/))
         e=nf90_put_var(f,id,a%Im(1:i_K0,0:var_N%pH1), start=(/1,1,2/))         
         do r = 1, mpi_sze-1
            pN0 = var_N%pH0_(r)
            pN1 = var_N%pH1_(r)
            mpi_tg = r
            call mpi_recv( x1e%Re(1,0), (pN1+1)*i_K0, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            call mpi_recv( x1e%Im(1,0), (pN1+1)*i_K0, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            e=nf90_put_var(f,id,x1e%Re(1:i_K0,0:pN1), start=(/1,pN0+1,1/))
            e=nf90_put_var(f,id,x1e%Im(1:i_K0,0:pN1), start=(/1,pN0+1,2/))
         end do
      else
         mpi_tg = mpi_rnk
         call mpi_send( a%Re(1,0), (var_N%pH1+1)*i_K0, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%Im(1,0), (var_N%pH1+1)*i_K0, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
      end if
#endif      
   end subroutine io_save_restress_even

   subroutine io_save_restress_odd(f,id,a)
      integer,     intent(in) :: f, id
      type (spec_xavg_odd), intent(in) :: a
      integer :: e
      
#ifndef _MPI
      e=nf90_put_var(f,id,a%Re(1:i_K1,0:i_NN1), start=(/1,1,1/))
      e=nf90_put_var(f,id,a%Im(1:i_K1,0:i_NN1), start=(/1,1,2/))

#else
      integer :: r, pN0,pN1

      if(mpi_rnk==0) then
         e=nf90_put_var(f,id,a%Re(1:i_K1,0:var_N%pH1), start=(/1,1,1/))
         e=nf90_put_var(f,id,a%Im(1:i_K1,0:var_N%pH1), start=(/1,1,2/))         
         do r = 1, mpi_sze-1
            pN0 = var_N%pH0_(r)
            pN1 = var_N%pH1_(r)
            mpi_tg = r
            call mpi_recv( x1o%Re(1,0), (pN1+1)*i_K1, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            call mpi_recv( x1o%Im(1,0), (pN1+1)*i_K1, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            e=nf90_put_var(f,id,x1o%Re(1:i_K1,0:pN1), start=(/1,pN0+1,1/))
            e=nf90_put_var(f,id,x1o%Im(1:i_K1,0:pN1), start=(/1,pN0+1,2/))
         end do
      else
         mpi_tg = mpi_rnk
         call mpi_send( a%Re(1,0), (var_N%pH1+1)*i_K1, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%Im(1,0), (var_N%pH1+1)*i_K1, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
      end if
#endif      
   end subroutine io_save_restress_odd
   
   subroutine io_save_restress_2d()
      character(4) :: cnum
      integer :: e, f
      integer :: md,nd,kd, ReImd, dims(4)
      integer :: mss,rss1,rss2

      write(cnum,'(I4.4)') i_restress_save

      !! Mean velocities
      if(mpi_rnk==0) then
         print*, ' Saving Reynolds averages and stresses (2D)'//cnum//'  t=', tim_t
         e=nf90_create('restress_2d'//cnum//'.cdf.dat', nf90_clobber, f)
         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'tstep', tim_step)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         e=nf90_put_att(f, nf90_global, 'gamma', d_gamma)

         e=nf90_def_dim(f, 'N', i_NN, nd)
         e=nf90_def_dim(f, 'M', i_M, md)
         e=nf90_def_dim(f, 'KK', i_KK, kd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         dims = (/kd,md,nd,ReImd/)
         call io_define_spec(f, 'umean_2d', dims, mss)
         call io_define_spec(f, 'remean1_2d', dims, rss1)
         call io_define_spec(f, 'remean2_2d', dims, rss2)

         e=nf90_enddef(f)
      end if

      call io_save_spec(f,mss, umean_2d)
      call io_save_spec(f,rss1,remean1_2d)
      call io_save_spec(f,rss2, remean2_2d)

      if(mpi_rnk==0)  &
         e=nf90_close(f)

   end subroutine io_save_restress_2d
   
   subroutine io_save_restress_filt()
      character(4) :: cnum
      integer :: e, f
      integer :: md,nd,kd, ReImd, dims(4)
      integer :: mss,rss1,rss2
      type(spec) :: umean_filt, remean1_filt, remean2_filt

      write(cnum,'(I4.4)') i_restress_save

      if(mpi_rnk==0) then
         print*, ' Saving Reynolds averages and stresses (filt)'//cnum//'  t=', tim_t
         e=nf90_create('restress_filt'//cnum//'.cdf.dat', nf90_clobber, f)
         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'tstep', tim_step)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         e=nf90_put_att(f, nf90_global, 'gamma', d_gamma)
         e=nf90_put_att(f, nf90_global, 'nx_c', i_nx_c)
         e=nf90_put_att(f, nf90_global, 'nz_c', i_nz_c)

         e=nf90_def_dim(f, 'N', i_NN, nd)
         e=nf90_def_dim(f, 'M', i_M, md)
         e=nf90_def_dim(f, 'KK', i_KK, kd)
         e=nf90_def_dim(f, 'ReIm', 2, ReImd)

         dims = (/kd,md,nd,ReImd/)
         call io_define_spec(f, 'umean_filt', dims, mss)
         call io_define_spec(f, 'remean1_filt', dims, rss1)
         call io_define_spec(f, 'remean2_filt', dims, rss2)

         e=nf90_enddef(f)
      end if

      call vel_restress_filt_calc(vel_c,i_nx_c,i_nz_c,umean_filt,remean1_filt,remean2_filt)

      call io_save_spec(f,mss, umean_filt)
      call io_save_spec(f,rss1,remean1_filt)
      call io_save_spec(f,rss2, remean2_filt)

      if(mpi_rnk==0)  &
         e=nf90_close(f)

   end subroutine io_save_restress_filt


!--------------------------------------------------------------------------
!  Save coll variable
!--------------------------------------------------------------------------
   subroutine io_define_mpt(f,nm,dims, id)
      integer,      intent(in) :: f, dims(4)
      character(*), intent(in) :: nm
      integer, intent(out) :: id
      integer :: e
      e=nf90_def_var(f, nm, nf90_double, dims, id)
      e=nf90_put_att(f, id,  'K', i_K)      
      e=nf90_put_att(f, id,  'MM', i_MM)
      e=nf90_put_att(f, id,  'NN', i_NN)
   end subroutine io_define_mpt

   subroutine io_save_mpt(f,id,a)
      integer,     intent(in) :: f, id
      type (mpt), intent(in) :: a
      integer :: e
      
#ifndef _MPI
      e=nf90_put_var(f,id,a%Re(1:i_K,0:i_M1,0:i_NN1), start=(/1,1,1,1/))
      e=nf90_put_var(f,id,a%Im(1:i_K,0:i_M1,0:i_NN1), start=(/1,1,1,2/))

#else
      integer :: r, pN0,pN1

      if(mpi_rnk==0) then
         e=nf90_put_var(f,id,a%Re(1:i_K,0:i_M1,0:var_N%pH1), start=(/1,1,1,1/))
         e=nf90_put_var(f,id,a%Im(1:i_K,0:i_M1,0:var_N%pH1), start=(/1,1,1,2/))         
         do r = 1, mpi_sze-1
            pN0 = var_N%pH0_(r)
            pN1 = var_N%pH1_(r)
            mpi_tg = r
            call mpi_recv( c1%Re(1,0,0), i_M*(pN1+1)*i_K, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            call mpi_recv( c1%Im(1,0,0), i_M*(pN1+1)*i_K, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            e=nf90_put_var(f,id,c1%Re(1:i_K,0:i_M1,0:pN1), start=(/1,1,pN0+1,1/))
            e=nf90_put_var(f,id,c1%Im(1:i_K,0:i_M1,0:pN1), start=(/1,1,pN0+1,2/))
         end do
      else
         mpi_tg = mpi_rnk
         call mpi_send( a%Re(1,0,0), i_M*(var_N%pH1+1)*i_K, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%Im(1,0,0), i_M*(var_N%pH1+1)*i_K, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
      end if
#endif      
   end subroutine io_save_mpt


   subroutine io_define_spec(f,nm,dims, id)
      integer,      intent(in) :: f, dims(4)
      character(*), intent(in) :: nm
      integer, intent(out) :: id
      integer :: e
      e=nf90_def_var(f, nm, nf90_double, dims, id)
      e=nf90_put_att(f, id,  'KK', i_KK)      
      e=nf90_put_att(f, id,  'MM', i_MM)
      e=nf90_put_att(f, id,  'NN', i_NN)
   end subroutine io_define_spec

   subroutine io_save_spec(f,id,a)
      integer,     intent(in) :: f, id
      type (spec), intent(in) :: a
      integer :: e
      
#ifndef _MPI
      e=nf90_put_var(f,id,a%Re(1:i_KK,0:i_M1,0:i_NN1), start=(/1,1,1,1/))
      e=nf90_put_var(f,id,a%Im(1:i_KK,0:i_M1,0:i_NN1), start=(/1,1,1,2/))

#else
      integer :: r, pN0,pN1

      if(mpi_rnk==0) then
         e=nf90_put_var(f,id,a%Re(1:i_KK,0:i_M1,0:var_N%pH1), start=(/1,1,1,1/))
         e=nf90_put_var(f,id,a%Im(1:i_KK,0:i_M1,0:var_N%pH1), start=(/1,1,1,2/))         
         do r = 1, mpi_sze-1
            pN0 = var_N%pH0_(r)
            pN1 = var_N%pH1_(r)
            mpi_tg = r
            call mpi_recv( s1%Re(1,0,0), i_M*(pN1+1)*i_KK, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            call mpi_recv( s1%Im(1,0,0), i_M*(pN1+1)*i_KK, mpi_double_precision,  &
               r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
            e=nf90_put_var(f,id,s1%Re(1:i_KK,0:i_M1,0:pN1), start=(/1,1,pN0+1,1/))
            e=nf90_put_var(f,id,s1%Im(1:i_KK,0:i_M1,0:pN1), start=(/1,1,pN0+1,2/))
         end do
      else
         mpi_tg = mpi_rnk
         call mpi_send( a%Re(1,0,0), i_M*(var_N%pH1+1)*i_KK, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%Im(1,0,0), i_M*(var_N%pH1+1)*i_KK, mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
      end if
#endif      
   end subroutine io_save_spec


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  save spectrum !var_mpt2spec
!--------------------------------------------------------------------------
   subroutine io_save_Uspec()
     
     type(spec) :: tmp
     double precision :: n_(0:i_NN1), m_(0:i_MM1)
     double precision :: n__(0:i_NN1), m__(0:i_MM1)
     double precision :: d, dRe, dIm
     character(4) :: cnum
     integer :: i
     _loop_kmn_vars
10   format(i4,1e20.12)
     
     call var_mpt2spec(vel_c,tmp)
     n_ = 0d0
     m_ = 0d0
     _loop_mn_begin
     mm = m
     if (m > i_MM1) mm = i_M - m
     do k=1,i_KK
        dRe=tmp%Re(k,m,n)
        dIm=tmp%Im(k,m,n) 
        d = dsqrt(dRe*dRe+dIm*dIm)
        n_(nn)  = max(d, n_(nn))
        m_(mm)  = max(d, m_(mm))
     end do
     _loop_mn_end
     
#ifdef _MPI
     call mpi_barrier(mpi_comm_world, mpi_er)
     call mpi_allreduce(n_(0), n__(0), i_NN, mpi_double_precision,  &
          mpi_max, mpi_comm_world, mpi_er)
     n_ = n__
     call mpi_allreduce(m_(0), m__(0), i_MM, mpi_double_precision,  &
          mpi_max, mpi_comm_world, mpi_er)
     m_ = m__
#endif
     if(mpi_rnk/=0) return
     write(cnum,'(I4.4)') io_save1
     open(11, status='unknown', file='vel_spec'//cnum//'.dat')
     write(11,*) '# t = ', tim_t
     write(11,*) '# m'
     do i = 0, i_MM1
        write(11,10) i, m_(i)      
     end do
     write(11,*)
     write(11,*) '# n'
     do i = 0, i_NN1
        write(11,10) i, n_(i)      
     end do
     close(11)
     
   end subroutine io_save_Uspec


!--------------------------------------------------------------------------
!  save spectrum !var_mpt2spec
!--------------------------------------------------------------------------
   subroutine io_save_spectrum()
      double precision :: n_(0:i_NN1), m_(0:i_MM1)
      double precision :: n__(0:i_NN1), m__(0:i_MM1)
      double precision :: d, dRe, dIm
      character(4) :: cnum
      integer :: i
      _loop_kmn_vars
   10 format(i4,1e20.12)
      
      n_ = 0d0
      m_ = 0d0
      _loop_kmn_begin
      mm = m
      if (m > i_MM1) mm = i_M - m
!      print*,mpi_rnk,k,mm,nn
      dRe=vel_c%Re(k,m,n)
      dIm=vel_c%Im(k,m,n) 
      d = dsqrt(dRe*dRe+dIm*dIm)
      n_(nn)  = max(d, n_(nn))
      m_(mm)  = max(d, m_(mm))
      _loop_kmn_end
      
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
      call mpi_allreduce(n_(0), n__(0), i_NN, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      n_ = n__
      call mpi_allreduce(m_(0), m__(0), i_MM, mpi_double_precision,  &
         mpi_max, mpi_comm_world, mpi_er)
      m_ = m__
#endif
      if(mpi_rnk/=0) return
      write(cnum,'(I4.4)') io_save1
      open(11, status='unknown', file='vel_spec'//cnum//'.dat')
      write(11,*) '# t = ', tim_t
      write(11,*) '# m'
      do i = 0, i_MM1
         write(11,10) i, m_(i)      
      end do
      write(11,*)
      write(11,*) '# n'
      do i = 0, i_NN1
         write(11,10) i, n_(i)      
      end do
      close(11)

   end subroutine io_save_spectrum

 
!--------------------------------------------------------------------------
!  write to energy file
!--------------------------------------------------------------------------
   subroutine io_write_energy(sp)
      type (phys), intent(in) :: sp
      double precision :: E

      call vel_energy(sp,E)
      
      if(mpi_rnk/=0) return
      write(io_KE,'(2e20.12)')  tim_t, E
      
      if(E>d_minE .or. tim_t<20d0) return
      print*, 'io_write_energy: Relaminarised!'
      open(99,file=oloc)
      close(99, status='delete')

   end subroutine io_write_energy

!--------------------------------------------------------------------------
!  write to uq file
!--------------------------------------------------------------------------
   subroutine io_write_uq(mpt_in)
      type (mpt), intent(in) :: mpt_in
      type (mpt)      :: fluc
      type (phys)      :: fluc_phy
      double precision :: u,q

      ! First copy mpt array
      fluc%Re(:,:,:) = mpt_in%Re(:,:,:)
      fluc%Im(:,:,:) = mpt_in%Im(:,:,:)
      
      ! Extract u and remove it from fluc.
      ! Note that (in first cpu) u[2,0,0] = mpt[2,0,0]=f(y), so no need to call mpt2spec
      if (var_N%pH0 == 0) then
        u = mpt_in%Re(2,0,0) + 1 ! Adding laminar flow
        
        ! Remove in fluc
        fluc%Re(2,0,0) = 0.0d0
        fluc%Im(2,0,0) = 0.0d0
      endif

#ifdef _MPI
        call mpi_barrier(mpi_comm_world, mpi_er)
#endif

      call vel_mpt2phys(fluc,fluc_phy)
      call vel_energy(fluc_phy,q)

      if(mpi_rnk/=0) return
      write(io_uq,'(3e20.12)')  tim_t, u, q
      
   end subroutine io_write_uq

!--------------------------------------------------------------------------
!  write to history file
!--------------------------------------------------------------------------
   subroutine io_write_history(sp)
      type (phys), intent(in) :: sp
      double precision :: H(4,i_H)
      integer :: i

      if(mpi_rnk/=0) return
      call vel_history(sp,0d0,H)
      do i=1,i_H
         write(io_HI,'(6e20.12)')  tim_t, 0d0,H(1,i),H(2,i),H(3,i),H(4,i)
      end do

   end subroutine io_write_history


   Subroutine io_writeVTK_xz(V,y,cnum,xs,zs)
  
     type(phys),intent(in) :: V
     double precision, intent(in) :: y
     integer, optional :: xs,zs
     character(4) :: cnum
     integer :: io,ix,iz,M_,N_
     character(10) :: s = 'unknown', a = 'sequential'
     
     if (.not. present(zs)) then
        xs=1
        zs=1
     end if

     M_=i_3M/xs
     N_=i_3N/zs

     if(mpi_rnk==0) then
     
        open(io,status=s,access=a, file='XZ'//cnum//'.vtk')
        
        write(io,'(A)') "# vtk DataFile Version 3.0"
        write(io,'(A)') "M x 2 x N"
        write(io,'(A)') "ASCII"
        write(io,'(A)') "DATASET RECTILINEAR_GRID"
        write(io,'(A,3i5)') "DIMENSIONS",M_+1,1,N_+1
        
        write(io,'(A,i5,A)') "X_COORDINATES ",M_+1, " FLOAT"
        do ix=0,i_3M,xs
           write(io,'(f9.5)') d_Lx*dble(ix)/(i_3M)
        end do
        
        write(io,'(A,i5,A)') "Y_COORDINATES ",1, " FLOAT"
        write(io,'(f9.5)') 0d0
        
        write(io,'(A,i5,A)') "Z_COORDINATES ",N_+1, " FLOAT"
        do ix=0,i_3N,zs
           write(io,'(f9.5)') d_Lz*dble(ix)/(i_3N)
        end do
        
        write(io,*) 
        
        write(io,'(A,i10)') "POINT_DATA",(N_+1)*(M_+1) 
        write(io,'(A)') "VECTORS velocity FLOAT"

        do iz=0,i_3N-1,zs
           do ix=0,i_3M-1,xs
              write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
           end do
           ix=0
           write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
        end do
        iz=0
        do ix=0,i_3M-1,xs
           write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
        end do
        ix=0
        write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
        
        close(io)
     END if
  
   end subroutine io_writeVTK_xz

   Subroutine io_writeVTK_txz(V,cnum,xs,zs,ct)
  
     type(phys),intent(in) :: V
     double precision, intent(in) :: ct
     integer, optional :: xs,zs
     character(4) :: cnum
     integer :: io,ix,iz,M_,N_,t
     character(10) :: s = 'unknown', a = 'sequential'
     double precision :: ke

     if (.not. present(zs)) then
        xs=1
        zs=1
     end if

     M_=i_3M/xs
     N_=i_3N/zs
     t=0

     if(mpi_rnk==0) then
     
        open(io,status=s,access=a, file='TXZ'//cnum//'.vtk')
        
        write(io,'(A)') "# vtk DataFile Version 3.0"
        write(io,'(A)') "M x 2 x N"
        write(io,'(A)') "ASCII"
        write(io,'(A)') "DATASET RECTILINEAR_GRID"
        write(io,'(A,3i5)') "DIMENSIONS",M_+1,1,N_+1
        
        write(io,'(A,i5,A)') "X_COORDINATES ",M_+1, " FLOAT"
        do ix=0,i_3M,xs
           write(io,'(f9.5)') d_Lx*dble(ix)/(i_3M)
        end do
        
        write(io,'(A,i5,A)') "Y_COORDINATES ",1, " FLOAT"
        write(io,'(f9.5)') 0d0
        
        write(io,'(A,i5,A)') "Z_COORDINATES ",N_+1, " FLOAT"
        do ix=0,i_3N,zs
           write(io,'(f9.5)') d_Lz*dble(ix)/(i_3N)
        end do
        
        write(io,*) 
        
        write(io,'(A,i10)') "POINT_DATA",(N_+1)*(M_+1) 
        write(io,'(A)') "SCALAR turb FLOAT"

        do iz=0,i_3N-1,zs
           do ix=0,i_3M-1,xs
              ke = epos(V%Re(:,iz,ix))
              if (ke > ct) t = t + 1
              write(io,'(3e13.5)') ke
           end do
           ix=0
           write(io,'(3e13.5)') epos(V%Re(:,iz,ix))
        end do
        iz=0
        do ix=0,i_3M-1,xs
           write(io,'(3e13.5)') epos(V%Re(:,iz,ix))
        end do
        ix=0
        write(io,'(3e13.5)') epos(V%Re(:,iz,ix))
        
        close(io)
     END if
     print*,"Turbulent fraction :", dble(t)/((N_+1)*(M_+1))
  
   end subroutine io_writeVTK_txz


   Subroutine io_writeVTK(V,cnum,xs,zs)
  
     type(phys),intent(in) :: V
     integer, optional :: xs,zs
     character(4) :: cnum
     integer :: io,ix,iy,iz,M_,N_
     double precision :: y
     character(10) :: s = 'unknown', a = 'sequential'
     
     if (.not. present(zs)) then
        xs=1
        zs=1
     end if
     
     M_=i_3M/xs
     N_=i_3N/zs
     
     if(mpi_rnk==0) then
        
        open(io,status=s,access=a, file='XYZ'//cnum//'.vtk')
        
        write(io,'(A)') "# vtk DataFile Version 3.0"
        write(io,'(A)') "Lx x 2 x Lz"
        write(io,'(A)') "ASCII"
        write(io,'(A)') "DATASET RECTILINEAR_GRID"
        write(io,'(A,3i5)') "DIMENSIONS",M_+1,i_Ny,N_+1
        
        write(io,'(A,i5,A)') "X_COORDINATES ",M_+1, " FLOAT"
        do ix=0,i_3M,xs
           write(io,'(f9.5)') d_Lx*dble(ix)/(i_3M)
        end do
        
        write(io,'(A,i5,A)') "Y_COORDINATES ",i_Ny, " FLOAT"
        do iy=0,i_Ny-1
           write(io,'(f9.5)') -1d0 + 2d0*dble(iy)/(i_Ny-1)
        end do
        
        
        write(io,'(A,i5,A)') "Z_COORDINATES ",N_+1, " FLOAT"
        do ix=0,i_3N,zs
           write(io,'(f9.5)') d_Lz*dble(ix)/(i_3N)
        end do
        
        write(io,*) 
        
        write(io,'(A,i10)') "POINT_DATA",(N_+1)*(M_+1)*i_Ny
        write(io,'(A)') "VECTORS velocity FLOAT"
        
        do iz=0,i_3N-1,zs
           do iy=0,i_Ny-1
              y=-1d0 + 2d0*dble(iy)/(i_Ny-1)
              do ix=0,i_3M-1,xs
                 write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
              end do
              ix=0
              write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
           end do
        end do
        iz=0
        do iy=0,i_Ny-1
           y=-1d0 + 2d0*dble(iy)/(i_Ny-1)
           do ix=0,i_3M-1,xs
              write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
           end do
           ix=0
           write(io,'(3e13.5)') velU(V%Re(:,iz,ix),y),velV(V%Re(:,iz,ix),y),velW(V%Re(:,iz,ix),y)
        end do
        close(io)
     END if
  
   end subroutine io_writeVTK
   

!**************************************************************************
 end module io
!**************************************************************************

