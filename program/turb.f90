#include "../parallel.h"
 module turb
!*************************************************************************
! this module calculates turbulent frac on a reduced grid
! possibly averaging over small length scales
! definitely averaging in time 
!*************************************************************************
   use mpif
   use velocity
   use netcdf
   implicit none
   save
   type tst
      double precision     :: Re(3,0:i_SN-1, 0:i_SMp-1,i_MT)
    end type tst
   
   integer :: turb_save2
   type(tst) :: tdata

 contains


   subroutine turb_sample_t(ind)
     !Samples energy onto grid
     type(phys), intent(in) :: ind
     _loop_mn_vars
     tdata%Re(:,:,:,2:i_MT)=tdata%Re(:,:,:,1:i_MT-1)
     tdata%Re(:,:,:,1)=0d0
     do m=0,i_SMp-1
        do n=0,i_3N-1
           nn=FLOOR(dble(n)/dble(i_DN))
           tdata%Re(:,nn,m,1) = tdata%Re(:,nn,m,1) + epos(ind%Re(:,n,m*i_DM))/i_DN
     _loop_mn_end
           
  end subroutine turb_sample_t
         
   subroutine turb_sample(ind)
     !Samples energy onto grid
     type(phys), intent(in) :: ind
     _loop_mn_vars
     tdata%Re(:,:,:,2:i_MT)=tdata%Re(:,:,:,1:i_MT-1)
     
     do m=0,i_SMp-1
        do n=0,i_SN-1
           tdata%Re(:,n,m,1) = eposC(ind%Re(:,n*i_DN,m*i_DM))
     _loop_mn_end
   end subroutine turb_sample

   subroutine turb_prepare()
     integer :: i
     do i=1,i_MT-1
        tdata%Re(:,:,:,i_MT)=tdata%Re(:,:,:,i_MT)+tdata%Re(:,:,:,i)
     end do
     tdata%Re(:,:,:,i_MT)=tdata%Re(:,:,:,i_MT)/i_MT
!     print*, tdata%Re(0:10,0,i_MT)
   end subroutine turb_prepare

!--------------------------------------------------------------------------
!  Save state
!--------------------------------------------------------------------------
   subroutine turb_save_state()
      character(4) :: cnum
      integer :: e, f
      integer :: md,nd, dims(3),cd
      integer :: ss 

      write(cnum,'(I4.4)') turb_save2!io_save2

      if(mpi_rnk==0) then
         print*, ' saving turb'//cnum//'  t=', tim_t
         e=nf90_create('turb'//cnum//'.cdf.dat', nf90_clobber, f)

         e=nf90_put_att(f, nf90_global, 't', tim_t)
         e=nf90_put_att(f, nf90_global, 'Re', d_Re)
         e=nf90_put_att(f, nf90_global, 'alpha', d_alpha)
         e=nf90_put_att(f, nf90_global, 'gamma', d_gamma)

         e=nf90_def_dim(f, 'N', i_SN, nd)
         e=nf90_def_dim(f, 'M', i_SM, md)
         e=nf90_def_dim(f, 'C', 3, cd)

         dims = (/cd,nd,md/)
         call turb_define_turb(f, 'uvw', dims, ss)

         e=nf90_enddef(f)
      end if
      call turb_prepare()
      call turb_save_turb(f,ss, tdata)

      if(mpi_rnk==0)  &
         e=nf90_close(f)
      turb_save2 = turb_save2+1

    end subroutine turb_save_state

!--------------------------------------------------------------------------
!  Save coll variable
!--------------------------------------------------------------------------
   subroutine turb_define_turb(f,nm,dims, id)
     integer,      intent(in) :: f, dims(3)
     character(*), intent(in) :: nm
     integer, intent(out) :: id
     integer :: e
     e=nf90_def_var(f, nm, nf90_double, dims, id)
     e=nf90_put_att(f, id,  'C', 3)
     e=nf90_put_att(f, id,  'MM', i_SM)
     e=nf90_put_att(f, id,  'NN', i_SN)
   end subroutine turb_define_turb
    
   subroutine turb_save_turb(f,id,a)
     integer,     intent(in) :: f, id
     type (tst), intent(in) :: a
     integer :: e
     type (tst) :: c1
     
#ifndef _MPI
      e=nf90_put_var(f,id,a%Re(1:3,0:i_SN-1,0:i_SM-1,10), start=(/1,1,1/))
#else
      integer :: r, pT0,pT1
      
      if(mpi_rnk==0) then
         pT1 = max((var_M%pH1+1)/i_DM - 1,0) 
         e=nf90_put_var(f,id,a%Re(1:3,0:i_SN-1, 0:pT1,10), start=(/1,1,1/))
      end if
      if (i_SM > 1) then
         if (mpi_rnk==0) then
            do r = 1, mpi_sze-1
               pT0 = (var_M%pH0_(r))/i_DM
               pT1 = (var_M%pH1_(r)+1)/i_DM - 1 
               mpi_tg = r
               call mpi_recv( c1%Re(1,0,0,1), 3*i_SM*(pT1+1), mpi_double_precision,  &
                    r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
               e=nf90_put_var(f,id,c1%Re(1:3,0:i_SN-1,0:pT1,1), start=(/1,1,pT0+1/))
            end do
         else
            mpi_tg = mpi_rnk
            pT1 = (var_M%pH1+1)/i_DM - 1
           call mpi_send( a%Re(1,0,0,10), 3*i_SM*(pT1+1), mpi_double_precision,  &
                 0, mpi_tg, mpi_comm_world, mpi_er)
         end if
      end if
#endif      
    end subroutine turb_save_turb

  end module turb
