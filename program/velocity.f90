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
  type(spec) :: umean_2d,remean1_2d,remean2_2d,spec_in
  !umean_2d = (umean,vmean,wmean), remean1_2d = (uumean,uvmean,uwmean), remean2_2d = (vvmean,wvmean,wwmean)
  type(spec_xavg_even) :: umean,wmean,uumean,uwmean,wwmean,vvmean
  type(spec_xavg_odd) :: vmean,uvmean,wvmean
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
      ! If u1_fix == .true. then we set vel_c2(2,0,0) = f(2) here
      if (s_u1_fixed) then
       vel_c2%Re(2,0,0)= d_u1_in
       vel_c2%Re(4,0,0)= 0d0
      endif
    end if

    call var_mpt_copy(vel_nl,vel_onl)
    call var_mpt_copy(vel_c2,vel_c)
  end subroutine vel_TS
 
  subroutine vel_precompute()
    call var_spec_init(spec_in) 
    call var_spec_init(umean_2d) 
    call var_spec_init(remean1_2d) 
    call var_spec_init(remean2_2d) 
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
    logical :: exist
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
      
    !! If (spec.cdf.in) file is present, spec_in =/= 0, otherwise spec_in = 0
    u%Re = u%Re + spec_in%Re
    u%Im = u%Im + spec_in%Im

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

  subroutine vel_restress_calc(i_count)
  !! This subroutine takes the current fields and updates the average values and the Reynolds stresses
  !! These are x-averaged, and the time-mean is updated.
  !! Uses global variables: i_count, vel_c as inputs. and umean,vmean,wmean,
  !! uumean,uvmean,uwmean,vvmean,wvmean,wwmean as outputs. The outputs are of the
  !! type 'spec_xavg_even' (shape = (i_K0,0:i_Np-1), uumean, etc.), 'spec_xavg_odd' (shape = (i_K1,0:i_Np-1), uvmean, etc.)

    type(spec) :: u,up,nl_tmp
    type(spec_xavg_even) :: ubar,wbar,uubar,uwbar,wwbar,vvbar
    type(spec_xavg_odd) :: vbar,uvbar,wvbar
    type(phys) :: p,mult_aux
    integer, intent(in) ::  i_count
    _loop_mn_vars

    ! Takes current mpt field and converts to spec
    call var_mpt2spec(vel_c,u)
    
    ! -------- MEAN VELOCITIES -------
    ! Calculate the x-avg
    ubar%Re = reshape(u%Re(1:i_K0,0,:),shape(ubar%Re))
    ubar%Im = reshape(u%Im(1:i_K0,0,:),shape(ubar%Im))
    vbar%Re = reshape(u%Re(i_K0+1:2*i_K0-1,0,:),shape(vbar%Re))
    vbar%Im = reshape(u%Im(i_K0+1:2*i_K0-1,0,:),shape(vbar%Im))
    wbar%Re = reshape(u%Re(2*i_K0:3*i_K0-1,0,:),shape(wbar%Re))
    wbar%Im = reshape(u%Im(2*i_K0:3*i_K0-1,0,:),shape(wbar%Im))

    ! Then update the mean
    umean%Re = umean%Re + (ubar%Re - umean%Re)/i_count  ! <u>_n
    umean%Im = umean%Im + (ubar%Im - umean%Im)/i_count  ! <u>_n
    vmean%Re = vmean%Re + (vbar%Re - vmean%Re)/i_count  ! <v>_n
    vmean%Im = vmean%Im + (vbar%Im - vmean%Im)/i_count  ! <v>_n
    wmean%Re = wmean%Re + (wbar%Re - wmean%Re)/i_count  ! <w>_n
    wmean%Im = wmean%Im + (wbar%Im - wmean%Im)/i_count  ! <w>_n

    ! -------- REYNOLDS STRESSES  -------
    ! <u'u'>, <u'v'> and <u'w'>
    ! Calculate u' = u-umean, etc
    up%Re(1:i_K0,0,:) = u%Re(1:i_K0,0,:) - umean%Re
    up%Im(1:i_K0,0,:) = u%Im(1:i_K0,0,:) - umean%Im
    up%Re(i_K0+1:2*i_K0-1,0,:) = u%Re(i_K0+1:2*i_K0-1,0,:) - vmean%Re 
    up%Im(i_K0+1:2*i_K0-1,0,:) = u%Im(i_K0+1:2*i_K0-1,0,:) - vmean%Im
    up%Re(2*i_K0:3*i_K0-1,0,:) = u%Re(2*i_K0:3*i_K0-1,0,:) - wmean%Re
    up%Im(2*i_K0:3*i_K0-1,0,:) = u%Im(2*i_K0:3*i_K0-1,0,:) - wmean%Im
    
    do m = 1,i_M1
        up%Re(1:i_K0,m,:) = u%Re(1:i_K0,m,:)
        up%Im(1:i_K0,m,:) = u%Im(1:i_K0,m,:)
        up%Re(i_K0+1:2*i_K0-1,m,:) = u%Re(i_K0+1:2*i_K0-1,m,:) 
        up%Im(i_K0+1:2*i_K0-1,m,:) = u%Im(i_K0+1:2*i_K0-1,m,:)
        up%Re(2*i_K0:3*i_K0-1,m,:) = u%Re(2*i_K0:3*i_K0-1,m,:)
        up%Im(2*i_K0:3*i_K0-1,m,:) = u%Im(2*i_K0:3*i_K0-1,m,:)
    end do
    
    ! Calculate u'u', u'v', u'w'
    call tra_spec2phys(up,p)

    _loop_phy_begin
    mult_aux%Re(1:i_K0,n,m)=nluw(p%Re(1:i_K0,n,m),p%Re(1:i_K0,n,m)) ! u'u'
    mult_aux%Re(i_K0+1:2*i_K0-1,n,m)=nluv(p%Re(1:i_K0,n,m),p%Re(i_K0+1:2*i_K0-1,n,m)) ! u'v'
    mult_aux%Re(2*i_K0:3*i_K0-1,n,m)=nluw(p%Re(1:i_K0,n,m),p%Re(2*i_K0:3*i_K0-1,n,m)) ! u'w'
    _loop_mn_end
    
    call tra_phys2spec(mult_aux,nl_tmp)

    ! x-average to get bar(u'u') etc.
    uubar%Re = reshape(nl_tmp%Re(1:i_K0,0,:),shape(uubar%Re))
    uubar%Im = reshape(nl_tmp%Im(1:i_K0,0,:),shape(uubar%Im))
    uvbar%Re = reshape(nl_tmp%Re(i_K0+1:2*i_K0-1,0,:),shape(uvbar%Re))
    uvbar%Im = reshape(nl_tmp%Im(i_K0+1:2*i_K0-1,0,:),shape(uvbar%Im))
    uwbar%Re = reshape(nl_tmp%Re(2*i_K0:3*i_K0-1,0,:),shape(uwbar%Re))
    uwbar%Im = reshape(nl_tmp%Im(2*i_K0:3*i_K0-1,0,:),shape(uwbar%Im))

    ! Update the mean
    uumean%Re = uumean%Re + (uubar%Re - uumean%Re)/i_count ! <u'u'>
    uumean%Im = uumean%Im + (uubar%Im - uumean%Im)/i_count ! <u'u'>
    uvmean%Re = uvmean%Re + (uvbar%Re - uvmean%Re)/i_count ! <u'v'>
    uvmean%Im = uvmean%Im + (uvbar%Im - uvmean%Im)/i_count ! <u'v'>
    uwmean%Re = uwmean%Re + (uwbar%Re - uwmean%Re)/i_count ! <u'w'>
    uwmean%Im = uwmean%Im + (uwbar%Im - uwmean%Im)/i_count ! <u'w'>

    ! Let's do the same, but for w'w', v'w'
    _loop_phy_begin
    mult_aux%Re(i_K0+1:2*i_K0-1,n,m)=nluv(p%Re(2*i_K0:3*i_K0-1,n,m),p%Re(i_K0+1:2*i_K0-1,n,m)) ! w'v'
    mult_aux%Re(2*i_K0:3*i_K0-1,n,m)=nluw(p%Re(2*i_K0:3*i_K0-1,n,m),p%Re(2*i_K0:3*i_K0-1,n,m)) ! w'w'
    _loop_mn_end
    
    call tra_phys2spec(mult_aux,nl_tmp)

    ! x-average to get bar(u'u') etc.
    wvbar%Re = reshape(nl_tmp%Re(i_K0+1:2*i_K0-1,0,:),shape(wvbar%Re))
    wvbar%Im = reshape(nl_tmp%Im(i_K0+1:2*i_K0-1,0,:),shape(wvbar%Im))
    wwbar%Re = reshape(nl_tmp%Re(2*i_K0:3*i_K0-1,0,:),shape(wwbar%Re))
    wwbar%Im = reshape(nl_tmp%Im(2*i_K0:3*i_K0-1,0,:),shape(wwbar%Im))

    ! Update the mean
    wvmean%Re = wvmean%Re + (wvbar%Re - wvmean%Re)/i_count ! <w'v'>
    wvmean%Im = wvmean%Im + (wvbar%Im - wvmean%Im)/i_count ! <w'v'>
    wwmean%Re = wwmean%Re + (wwbar%Re - wwmean%Re)/i_count ! <w'w'>
    wwmean%Im = wwmean%Im + (wwbar%Im - wwmean%Im)/i_count ! <w'w'>

    ! And finally for v'v'
    _loop_phy_begin
    mult_aux%Re(1:i_K0,n,m)=nlvv(p%Re(i_K0+1:2*i_K0-1,n,m)) ! v'v'
    _loop_mn_end
    
    call tra_phys2spec(mult_aux,nl_tmp)

    ! x-average to get bar(u'u') etc.
    vvbar%Re = reshape(nl_tmp%Re(1:i_K0,0,:),shape(vvbar%Re))
    vvbar%Im = reshape(nl_tmp%Im(1:i_K0,0,:),shape(vvbar%Im))

    ! Update the mean
    vvmean%Re = vvmean%Re + (vvbar%Re - vvmean%Re)/i_count ! <v'v'>
    vvmean%Im = vvmean%Im + (vvbar%Im - vvmean%Im)/i_count ! <v'v'>

  end subroutine vel_restress_calc

  subroutine vel_restress_reset()
    umean%Re = 0.0*umean%Re
    umean%Im = 0.0*umean%Im
    vmean%Re = 0.0*vmean%Re
    vmean%Im = 0.0*vmean%Im
    wmean%Re = 0.0*wmean%Re
    wmean%Im = 0.0*wmean%Im
    uumean%Re = 0.0*uumean%Re
    uumean%Im = 0.0*uumean%Im
    uvmean%Re = 0.0*uvmean%Re
    uvmean%Im = 0.0*uvmean%Im
    uwmean%Re = 0.0*uwmean%Re
    uwmean%Im = 0.0*uwmean%Im
    wvmean%Re = 0.0*wvmean%Re
    wvmean%Im = 0.0*wvmean%Im
    wwmean%Re = 0.0*wwmean%Re
    wwmean%Im = 0.0*wwmean%Im
    vvmean%Re = 0.0*vvmean%Re
    vvmean%Im = 0.0*vvmean%Im
  end subroutine vel_restress_reset

  subroutine vel_restress_2d_calc(i_count)
  !! This subroutine takes the current fields and updates the average values and the Reynolds stresses
  !! These are full 2d, and the time-mean is updated.
  !! Uses global variables: i_count, vel_c as inputs. and umean_2d,vmean_2d,wmean_2d
  !! uumean_2d,uvmean_2d,uwmean_2d,vvmean_2d,wvmean_2d,wwmean_2d as outputs. The outputs are of the
  !! type 'spec' 

    type(spec) :: u,up,nl_tmp
    type(phys) :: p,mult_aux
    integer, intent(in) :: i_count
    _loop_mn_vars

    ! Takes current mpt field and converts to spec
    call var_mpt2spec(vel_c,u)
    
    ! -------- MEAN VELOCITIES -------

    ! Then update the mean
    umean_2d%Re = umean_2d%Re + (u%Re - umean_2d%Re)/i_count 
    umean_2d%Im = umean_2d%Im + (u%Im - umean_2d%Im)/i_count  

    ! -------- REYNOLDS STRESSES  -------
    ! <u'u'>, <u'v'> and <u'w'>
    ! Calculate u' = u-umean, etc
    up%Re = u%Re - umean_2d%Re
    up%Im = u%Im - umean_2d%Im
    
    ! Calculate u'u', u'v', u'w'
    call tra_spec2phys(up,p)

    _loop_phy_begin
    mult_aux%Re(1:i_K0,n,m)=nluw(p%Re(1:i_K0,n,m),p%Re(1:i_K0,n,m)) ! u'u'
    mult_aux%Re(i_K0+1:2*i_K0-1,n,m)=nluv(p%Re(1:i_K0,n,m),p%Re(i_K0+1:2*i_K0-1,n,m)) ! u'v'
    mult_aux%Re(2*i_K0:3*i_K0-1,n,m)=nluw(p%Re(1:i_K0,n,m),p%Re(2*i_K0:3*i_K0-1,n,m)) ! u'w'
    _loop_mn_end
    
    call tra_phys2spec(mult_aux,nl_tmp)

    ! Update the mean
    remean1_2d%Re = remean1_2d%Re + (nl_tmp%Re - remean1_2d%Re)/i_count 
    remean1_2d%Im = remean1_2d%Im + (nl_tmp%Im - remean1_2d%Im)/i_count 

    ! Let's do the same, but for w'w', v'w'
    _loop_phy_begin
    mult_aux%Re(1:i_K0,n,m)=nlvv(p%Re(i_K0+1:2*i_K0-1,n,m)) ! v'v'
    mult_aux%Re(i_K0+1:2*i_K0-1,n,m)=nluv(p%Re(2*i_K0:3*i_K0-1,n,m),p%Re(i_K0+1:2*i_K0-1,n,m)) ! w'v'
    mult_aux%Re(2*i_K0:3*i_K0-1,n,m)=nluw(p%Re(2*i_K0:3*i_K0-1,n,m),p%Re(2*i_K0:3*i_K0-1,n,m)) ! w'w'
    _loop_mn_end
    
    call tra_phys2spec(mult_aux,nl_tmp)
    
    ! Update the mean
    remean2_2d%Re = remean2_2d%Re + (nl_tmp%Re - remean2_2d%Re)/i_count 
    remean2_2d%Im = remean2_2d%Im + (nl_tmp%Im - remean2_2d%Im)/i_count 

  end subroutine vel_restress_2d_calc

  subroutine vel_restress_2d_reset()
    umean_2d%Re = 0.0*umean_2d%Re
    umean_2d%Im = 0.0*umean_2d%Im
    remean1_2d%Re = 0.0*remean1_2d%Re
    remean1_2d%Im = 0.0*remean1_2d%Im
    remean2_2d%Re = 0.0*remean2_2d%Re
    remean2_2d%Im = 0.0*remean2_2d%Im
  end subroutine vel_restress_2d_reset
  
  subroutine vel_restress_filt_calc(vel_c,nx_c,nz_c,umean_filt,remean1_filt,remean2_filt)
  !! This subroutine takes the current fields and calculates the mean
  !! via lowpass filter. Then calculates Reynolds averages.
  !! These are full 2d in x-z plane. No x-averaging
  !! Inputs: vel_c (mpt), cutoff mode number nx_c in the x-direction and nz_c in
  !! the z direction
  !! Outputs: filtered mean velocity umean_filt, and reynolds stresses
  !! remean1_filt, remean2_filt (all of spec type)

    type(mpt), intent(in) :: vel_c
    type(spec) :: u,up,nl_tmp
    type(spec), intent(out) :: umean_filt,remean1_filt,remean2_filt
    type(phys) :: p,mult_aux
    integer, intent(in) :: nx_c,nz_c
    integer :: nx
    _loop_mn_vars

    ! Takes current mpt field and converts to spec
    call var_mpt2spec(vel_c,u)
    
    ! -------- MEAN VELOCITIES -------

    ! Filter the velocity field
    _loop_mn_begin
    if (m<i_MM) then
        nx = m
    else if (m.ge.i_MM) then
        nx = abs(m-i_M)
    end if 
    if ((nn.ge.nz_c).or.(nx.ge.nx_c)) then
        umean_filt%Re(:,m,n) = 0.0*u%Re(:,m,n) 
        umean_filt%Im(:,m,n) = 0.0*u%Im(:,m,n)
    else
        umean_filt%Re(:,m,n) = u%Re(:,m,n) 
        umean_filt%Im(:,m,n) = u%Im(:,m,n)
    end if
    _loop_mn_end

    ! -------- REYNOLDS STRESSES  -------
    ! <u'u'>, <u'v'> and <u'w'>
    ! Calculate u' = u-umean, etc
    up%Re = u%Re - umean_filt%Re
    up%Im = u%Im - umean_filt%Im
    
    ! Calculate u'u', u'v', u'w'
    call tra_spec2phys(up,p)

    _loop_phy_begin
    mult_aux%Re(1:i_K0,n,m)=nluw(p%Re(1:i_K0,n,m),p%Re(1:i_K0,n,m)) ! u'u'
    mult_aux%Re(i_K0+1:2*i_K0-1,n,m)=nluv(p%Re(1:i_K0,n,m),p%Re(i_K0+1:2*i_K0-1,n,m)) ! u'v'
    mult_aux%Re(2*i_K0:3*i_K0-1,n,m)=nluw(p%Re(1:i_K0,n,m),p%Re(2*i_K0:3*i_K0-1,n,m)) ! u'w'
    _loop_mn_end
    
    call tra_phys2spec(mult_aux,nl_tmp)

    ! Filter 
    _loop_mn_begin
    if (m<i_MM) then
        nx = m
    else if (m.ge.i_MM) then
        nx = abs(m-i_M)
    end if 
    if ((nn.ge.nz_c).or.(nx.ge.nx_c)) then
        remean1_filt%Re(:,m,n) = 0.0*nl_tmp%Re(:,m,n) 
        remean1_filt%Im(:,m,n) = 0.0*nl_tmp%Im(:,m,n)
    else
        remean1_filt%Re(:,m,n) = nl_tmp%Re(:,m,n) 
        remean1_filt%Im(:,m,n) = nl_tmp%Im(:,m,n)
    end if
    _loop_mn_end

    ! Let's do the same, but for w'w', v'w'
    _loop_phy_begin
    mult_aux%Re(1:i_K0,n,m)=nlvv(p%Re(i_K0+1:2*i_K0-1,n,m)) ! v'v'
    mult_aux%Re(i_K0+1:2*i_K0-1,n,m)=nluv(p%Re(2*i_K0:3*i_K0-1,n,m),p%Re(i_K0+1:2*i_K0-1,n,m)) ! w'v'
    mult_aux%Re(2*i_K0:3*i_K0-1,n,m)=nluw(p%Re(2*i_K0:3*i_K0-1,n,m),p%Re(2*i_K0:3*i_K0-1,n,m)) ! w'w'
    _loop_mn_end
    
    call tra_phys2spec(mult_aux,nl_tmp)
    
    ! Update the mean
    ! Filter 
    _loop_mn_begin
    if (m<i_MM) then
        nx = m
    else if (m.ge.i_MM) then
        nx = abs(m-i_M)
    end if 
    if ((nn.ge.nz_c).or.(nx.ge.nx_c)) then
        remean2_filt%Re(:,m,n) = 0.0*nl_tmp%Re(:,m,n) 
        remean2_filt%Im(:,m,n) = 0.0*nl_tmp%Im(:,m,n)
    else
        remean2_filt%Re(:,m,n) = nl_tmp%Re(:,m,n) 
        remean2_filt%Im(:,m,n) = nl_tmp%Im(:,m,n)
    end if
    _loop_mn_end

  end subroutine vel_restress_filt_calc

end module velocity

