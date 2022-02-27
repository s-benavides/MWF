!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
module velocity
  !*************************************************************************
  use mpif
  use variables
  use transform
  use modes
  implicit none
  save
  
  type (phys) :: vel_p
  type (mpt)  :: vel_c,vel_c2,vel_onl,vel_nl,lhs,rhs
  logical :: nst
  
contains
  
  subroutine vel_TS()
    _loop_kmn_vars
    call vel_nonlinear()
    if (nst) then
       call var_mpt_copy(vel_nl,vel_onl)
       nst = .false.
    end if
    _loop_kmn_begin
    vel_c2%Re(k,m,n)=rhs%Re(k,m,n)*vel_c%Re(k,m,n)
    vel_c2%Im(k,m,n)=rhs%Re(k,m,n)*vel_c%Im(k,m,n)
    vel_c2%Re(k,m,n)=vel_c2%Re(k,m,n) + d_dt/2d0*(3d0*vel_nl%Re(k,m,n)-vel_onl%Re(k,m,n))
    vel_c2%Im(k,m,n)=vel_c2%Im(k,m,n) + d_dt/2d0*(3d0*vel_nl%Im(k,m,n)-vel_onl%Im(k,m,n))
    vel_c2%Re(k,m,n)=vel_c2%Re(k,m,n) / lhs%Re(k,m,n)
    vel_c2%Im(k,m,n)=vel_c2%Im(k,m,n) / lhs%Re(k,m,n)
    _loop_kmn_end
    if (var_N%pH0 == 0) then
       vel_c2%Re(1,0,0)=0d0
    end if
    call var_mpt_copy(vel_nl,vel_onl)
    call var_mpt_copy(vel_c2,vel_c)
  end subroutine vel_TS
 
  subroutine vel_precompute()
   
    call var_mpt_init(vel_c)
    call var_mpt_init(vel_nl)
    call var_LHSRHS(lhs,rhs)
    nst = .true.

  end subroutine vel_precompute

   !-------------------------------------------------------------------------
   subroutine vel_clk_time(t)
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
   end subroutine vel_clk_time
    
!-------------------------------------------------------------------------

  subroutine vel_imposesym()

    if (s_reflect)   call var_reflect(vel_c)
    if (s_uvreflect) call var_uvreflect(vel_c)
    if (s_wreflect)  call var_wreflect(vel_c)

  end subroutine vel_imposesym

  subroutine vel_addlam(u)
    type (spec), intent(inout) :: u
    if (var_N%pH0 == 0) then 
       u%Re(2,0,0) = u%Re(2,0,0) + dcos(d_theta)
       u%Re(2*i_K0+1,0,0) = u%Re(2*i_K0+1,0,0) + dsin(d_theta)
    end if
  end subroutine vel_addlam
  
  subroutine vel_mpt2phys(u,p)
    type(mpt), intent(in) :: u
    type(phys), intent(out) :: p
    type(spec) :: tmp
    call var_mpt2spec(u,tmp)
    call tra_spec2phys(tmp,p)
  end subroutine vel_mpt2phys
  
  subroutine vel_nonlinear()
    real :: ci,co
    type (spec) :: u,ux,uz
    type (phys) :: p,px,pz,ans
    _loop_mn_vars
    
    call var_mpt2spec(vel_c,u)
    call vel_addlam(u)
    call var_spec_grad(u,ux,uz)

    call tra_spec2phys(u,p)
    call tra_spec2phys(ux,px)
    call tra_spec2phys(uz,pz)

    _loop_phy_begin
    ans%Re(:,n,m)=udotgradu(p%Re(:,n,m),px%Re(:,n,m),pz%Re(:,n,m))
    _loop_mn_end
    
    call tra_phys2spec(ans,u)
    call var_spec2mpt(u,vel_nl)

  end subroutine vel_nonlinear
  
  function epos(u) result(e)
    double precision :: e
    double precision,intent(in) :: u(i_KK)
    double precision :: udotu(i_K0)
    
    udotu = nluw(u(1:i_K0),u(1:i_K0))
    e = udotu(1)
    udotu = nlvv(u(i_K0+1:2*i_K0-1))
    e = e + udotu(1)
    udotu= nluw(u(2*i_K0:i_KK),u(2*i_K0:i_KK))
    e = e +udotu(1)
  end function epos

  function eposC(u) result(e)
    double precision :: e(3)
    double precision,intent(in) :: u(i_KK)
    double precision :: udotu(i_K0)
    double precision :: t(i_KK)
    t=u
!    t(1) = 0d0
!    t(2*i_K0) = 0d0
    udotu = nluw(t(1:i_K0),t(1:i_K0))
    e(1) = udotu(1)
    udotu = nlvv(t(i_K0+1:2*i_K0-1))
    e(2) = udotu(1)
    udotu= nluw(t(2*i_K0:i_KK),t(2*i_K0:i_KK))
    e(3) = udotu(1)
  end function eposC

  subroutine vel_energy(a,e)
    type (phys), intent(in) :: a
    double precision, intent(out) :: e
    double precision :: e_
    _loop_mn_vars
    e_=0
    _loop_phy_begin
    e_ = e_ + epos(a%Re(:,n,m))
    _loop_mn_end
    e_ = e_ /( i_3M * i_3N) 
#ifdef _MPI
    call mpi_allreduce( e_, e, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
#endif

  end subroutine vel_energy

    subroutine vel_history(V,y,ans)
    
    type(phys), intent(in) :: V
    double precision, intent(out) :: ans(4,i_H)
    integer :: j,ji
    double precision,intent(in) :: y

    if(mpi_rnk/=0) return 
    do j=1,i_H
       ji=int((j-1)*i_3N/i_H)+1
       ans(1,j)=dble(ji)*d_Lz/i_3N
       ans(2,j)=velU(V%Re(:,ji,0),y)
       ans(3,j)=velV(V%Re(:,ji,0),y)
       ans(4,j)=velW(V%Re(:,ji,0),y)
    end do
    
  end subroutine vel_history

end module velocity

