!*************************************************************************
! main.f90 (executable)
!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
PROGRAM MAIN
!*************************************************************************
  use variables
  use mpif
  use transform
  use velocity
  use io
  implicit none
  real :: d_start, d_stop
   
   call initialise()
   
   if (mpi_rnk==0) print*, 'time',tim_t,' tstep', tim_step,'  io_save1', io_save1
   do while(.not.terminate())
      call var_null()
      call vel_imposesym()
      ! Update Reynolds average and stress calculations
      i_count       = i_count + 1
      d_avg_time    = d_avg_time + d_dt
      if (s_restress_xavg) call vel_restress_calc(i_count) 
      if (s_restress_2d) call vel_restress_2d_calc(i_count) 
      call io_write2files()
      call vel_TS()
      tim_t         = tim_t    + d_dt
      tim_step      = tim_step + 1
   end do
   call cleanup()

   contains

   subroutine initialise()
      character(4) :: cnum
      integer :: i_tmp,row
      logical :: file_exist
      character(200) :: loadfile

      call mpi_precompute()
      if(mpi_rnk==0) then
         call system('touch PRECOMPUTING')
         call system('echo $HOSTNAME > HOST')
      end if

      if(mpi_rnk==0)  print*, 'precomputing function requisites...'
      call par_precompute()
      call var_precompute()
      call tra_precompute()
      call vel_precompute()
      call  io_precompute()

      if(mpi_rnk==0)  print*, 'loading state...'

       call clk_time(d_start)
       call io_load_state()
       if ((i_restress_save.gt.0).and.((s_restress_xavg).or.(s_restress_2d))) then
           ! Load count, and avg_time
           if (mpi_rnk==0) then
             open(1,file='count_avg_time.txt',status='unknown')
             do row=1,i_restress_save
                 read(1,*) i_tmp, i_count, d_avg_time
             end do
             close(1)

             print*, 'Continuing previous averages. count = ',i_count,', avg_time = ',d_avg_time
           endif
#ifdef _MPI
           call mpi_bcast(i_count,1,mpi_integer,0,mpi_comm_world,mpi_er)
           call mpi_bcast(d_avg_time,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,mpi_er)
#endif
          ! Load averages into arrays
           if (s_restress_xavg) call io_load_xavg()
           if (s_restress_2d) then 
               write(cnum,'(I4.4)') i_restress_save
               loadfile = 'restress_2d'//cnum//'.cdf.dat'
               call io_load_spec(loadfile,'umean_2d',umean_2d)
               call io_load_spec(loadfile,'KE_mean',KE_mean)
               call io_load_spec(loadfile,'remean1_2d',remean1_2d)
               call io_load_spec(loadfile,'remean2_2d',remean2_2d)
           endif 
       else ! Load averages
          ! If no previous save, then set Reynolds averaging counts to zero 
          i_count = 0 
          d_avg_time = 0.0
       endif ! Load averages
       call clk_time(d_stop)
#ifdef _CPUTIME
       if(mpi_rnk==0)       print*, ' CPU loadtime  = ', int((d_stop-d_start)/6d1), ' mins.'
#else
       if(mpi_rnk==0)       print*, ' WALL loadtime = ', int((d_stop-d_start)/6d1), ' mins.'
#endif 
      

      if(mpi_rnk==0)  print*, 'initialising output files...'
      call io_openfiles()

      if(mpi_rnk==0) then
         open (99, file='PRECOMPUTING')
         close(99, status='delete')
         inquire(file='OLOCATION', exist=file_exist)
         if (file_exist) then
            open (99, file='OLOCATION')
            read (99,'(A)') oloc
            close(99, status='delete')
         end if
         open (99, file=oloc)
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
         close(99)
      end if

      call clk_time(d_start)

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

      if (mpi_rnk==0)   print*, 'timestepping.....'
	end subroutine initialise

!-------------------------------------------------------------------------
 subroutine clk_time(t)
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
   end subroutine clk_time

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!  Termaination conditions
!-------------------------------------------------------------------------
   logical function terminate()
      logical :: file_exist
            
      if(mpi_rnk==0) then
         terminate = .false.
      
         if(tim_step==i_maxtstep) then
            terminate = .true.
            print*, 'maxtstep reached!'
         end if

         if(d_maxt>0d0 .and. tim_t>=d_maxt) then
            terminate = .true.
            print*, 'maxt reached!'
         end if

         call clk_time(d_stop)
         if(dble(d_stop-d_start)/36d2 >= d_cpuhours) then
            terminate = .true.
            print*, 'cpuhours reached!'
         end if

         if( modulo(tim_step,i_save_rate2)==0) then
            inquire(file=oloc, exist=file_exist)
            if(.not. file_exist) then
               terminate = .true.
               print*, 'RUNNING deleted !'
            end if
         end if
      end if

#ifdef _MPI
      call mpi_bcast(terminate,1,mpi_logical, 0,mpi_comm_world,mpi_er)
#endif

   end function terminate


    subroutine cleanup()
      character(4) :: cnum
      logical :: file_exist
   
      if(mpi_rnk==0) then
         print*, 'cleanup...'
         call clk_time(d_stop)
         print*, ' sec/step  = ', (d_stop-d_start)/real(tim_step)
#ifdef _CPUTIME
         print*, ' CPU time  = ', int((d_stop-d_start)/6d1), ' mins.'
#else
         print*, ' WALL time = ', int((d_stop-d_start)/6d1), ' mins.'
#endif 
      end if
      
      call io_closefiles()
      call clk_time(d_start)
      call io_save_state()
        ! Reynolds stresses
        if ((s_restress_xavg).or.(s_restress_2d).or.(s_restress_filt)) then
           i_restress_save = i_restress_save + 1 
           ! Exports Reynolds averages and stresses
           if (s_restress_xavg) call io_save_xavg()
           if (s_restress_2d) call io_save_restress_2d()
           if (s_restress_filt) call io_save_restress_filt()

            ! Save i_count and d_avg_time
            if (mpi_rnk==0) then
            open(1,file='count_avg_time.txt',position='append')
             write(1,10) i_restress_save, i_count, d_avg_time
       10    format( I4.4, I20.12,E20.12 )
             close(1)
            endif
        endif ! Re stresses
       call clk_time(d_stop)
#ifdef _CPUTIME
       if(mpi_rnk==0)       print*, ' CPU savetime  = ', int((d_stop-d_start)/6d1), ' mins.'
#else
       if(mpi_rnk==0)       print*, ' WALL savetime = ', int((d_stop-d_start)/6d1), ' mins.'
#endif 

      if(mpi_rnk==0) then
         inquire(file=oloc, exist=file_exist)
         if(file_exist) open(99, file='RUNNING')
         if(file_exist) close(99, status='delete')
      end if      

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
      call mpi_finalize(mpi_er)
#endif
      if(mpi_rnk==0) print*, '...done!'

   end subroutine cleanup
         
!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************

