#include "../parallel.h"
 module variables
!*************************************************************************
   use mpif
   use parameters
   implicit none
   save
   type phys
      double precision     :: Re(i_KK,0:i_3N-1, 0:i_Mp-1)
   end type phys
   
   type mpt
      double precision     :: Re(i_K,0:i_M1, 0:i_Np-1)
      double precision     :: Im(i_K,0:i_M1, 0:i_Np-1)
   end type mpt 

   type spec
      double precision     :: Re(i_KK,0:i_M1, 0:i_Np-1)
      double precision     :: Im(i_KK,0:i_M1, 0:i_Np-1)
   end type spec
   
   type spec_xavg_even
      double precision     :: Re(i_K0, 0:i_Np-1)
      double precision     :: Im(i_K0, 0:i_Np-1)
   end type spec_xavg_even
   
   type spec_xavg_odd
      double precision     :: Re(i_K0-1, 0:i_Np-1)
      double precision     :: Im(i_K0-1, 0:i_Np-1)
   end type spec_xavg_odd
   
   
   type tran
      double precision     :: Re(i_KK,0:i_NN1, 0:i_Mp-1)
      double precision     :: Im(i_KK,0:i_NN1, 0:i_Mp-1)
   end type tran

   type harm
      integer              :: pH0,pH1, pH0_(0:_Np-1),pH1_(0:_Np-1)
   end type harm


   type (harm)               :: var_M,var_N
   double precision, private :: ad_n2(0:i_NN1)
   double precision, private :: ad_n1(0:i_NN1)
   double precision, private :: ad_m1(0:i_M1)
   double precision, private :: ad_m2(0:i_M1)
   double precision, private :: ad_k1(0:i_K0-1)
   double precision, private :: ad_k2(0:i_K0-1)

 contains

   subroutine var_uvreflect(c)
     type (mpt), intent(inout) :: c
     _loop_kmn_vars
     _loop_mn_begin
     if (m == 0) then
        if (nn==0) then
           c%Re(1,m,n)= 0d0
           c%Re(3,m,n)= 0d0
           c%Re(5,m,n)= 0d0
           c%Re(7,m,n)= 0d0
        else
           c%Re(1,m,n)= 0d0
           c%Im(1,m,n)= 0d0
           c%Re(3,m,n)= 0d0
           c%Im(3,m,n)= 0d0
           c%Re(5,m,n)= 0d0
           c%Im(5,m,n)= 0d0
           c%Re(7,m,n)= 0d0
           c%Im(7,m,n)= 0d0
        end if
     else if (m > i_MM1) then
        mm = abs(m - i_M)
        c%Re(1,m,n)=-c%Re(1,mm,n)
        c%Im(1,m,n)=-c%Im(1,mm,n)
        c%Re(2,m,n)= c%Re(2,mm,n)
        c%Im(2,m,n)= c%Im(2,mm,n)
        c%Re(3,m,n)=-c%Re(3,mm,n)
        c%Im(3,m,n)=-c%Im(3,mm,n)
        c%Re(4,m,n)= c%Re(4,mm,n)
        c%Im(4,m,n)= c%Im(4,mm,n)
        c%Re(5,m,n)=-c%Re(5,mm,n)
        c%Im(5,m,n)=-c%Im(5,mm,n)
        c%Re(6,m,n)= c%Re(6,mm,n)
        c%Im(6,m,n)= c%Im(6,mm,n)
        c%Re(7,m,n)=-c%Re(7,mm,n)
        c%Im(7,m,n)=-c%Im(7,mm,n)
     end if
     _loop_mn_end
   end subroutine var_uvreflect

   subroutine var_wreflect(c)
     type (mpt), intent(inout) :: c
     _loop_kmn_vars
     _loop_mn_begin
     if (nn==0) then
        if (m==0) then
           c%Re(5:7,m,n)= 0d0
        else
           c%Re(1:4,m,n)= 0d0
           c%Im(1:4,m,n)= 0d0
        end if
     else if (m > i_MM1) then
        mm = abs(m - i_M)
        c%Re(1,m,n)=-c%Re(1,mm,n)
        c%Im(1,m,n)= c%Im(1,mm,n)
        c%Re(2,m,n)=-c%Re(2,mm,n)
        c%Im(2,m,n)= c%Im(2,mm,n)
        c%Re(3,m,n)=-c%Re(3,mm,n)
        c%Im(3,m,n)= c%Im(3,mm,n)
        c%Re(4,m,n)=-c%Re(4,mm,n)
        c%Im(4,m,n)= c%Im(4,mm,n)
        c%Re(5,m,n)= c%Re(5,mm,n)
        c%Im(5,m,n)=-c%Im(5,mm,n)
        c%Re(6,m,n)= c%Re(6,mm,n)
        c%Im(6,m,n)=-c%Im(6,mm,n)
        c%Re(7,m,n)= c%Re(7,mm,n)
        c%Im(7,m,n)=-c%Im(7,mm,n)
     end if
     _loop_mn_end
   end subroutine var_wreflect

   subroutine var_reflect(c)
     type (mpt), intent(inout) :: c
     _loop_kmn_vars
     _loop_mn_begin
     if (m == 0 .and. nn==0) then
        c%Re(1,m,n)= 0d0
        c%Re(3,m,n)= 0d0
        c%Re(6,m,n)= 0d0
     else
        c%Im(1,m,n)= 0d0
        c%Re(2,m,n)= 0d0
        c%Im(3,m,n)= 0d0
        c%Re(4,m,n)= 0d0
        c%Re(5,m,n)= 0d0
        c%Im(6,m,n)= 0d0
        c%Re(7,m,n)= 0d0
     end if
     _loop_mn_end
   end subroutine var_reflect

   subroutine var_maskmpt(c)
      type (mpt), intent(inout) :: c
      double precision :: mval,mn
   _loop_mn_vars
   _loop_mn_begin
      mm = m
      if (mm > i_MM1) mm = m - i_M
      mn = dsqrt(1d0*mm * mm + 1d0*nn * nn)
      mval = 0.005 * 10 ** (-mn/3d0)
      c%Re(:,m,n) = c%Re(:,m,n) * mval
      c%Im(:,m,n) = c%Im(:,m,n) * mval
   _loop_mn_end
   end subroutine var_maskmpt

   subroutine var_doubleX(c)
      type (mpt), intent(inout) :: c
      type (mpt) :: tmp
      double precision :: mval,mn
_loop_mn_vars
      call var_mpt_init(tmp)
_loop_mn_begin
      mm = 2*m
      if (m > i_MM1) then
         mm = i_M - (2*(m - i_M))
         if (mm < i_MM1) cycle
      else
         if (mm > i_MM1) cycle
      end if
      tmp%Re(:,mm,n) = c%Re(:,m,n)
      tmp%Im(:,mm,n) = c%Im(:,m,n)
_loop_mn_end
   call var_mpt_copy(tmp,c)
   end subroutine var_doubleX

   subroutine var_quadX(c)
      type (mpt), intent(inout) :: c
      type (mpt) :: tmp
      double precision :: mval,mn
_loop_mn_vars
      call var_mpt_init(tmp)
_loop_mn_begin
      mm = 4*m
      if (m > i_MM1) then
         mm = i_M - (4*(i_M - m))
         if (mm < i_MM1) cycle
      else
         if (mm > i_MM1) cycle
      end if
      tmp%Re(:,mm,n) = c%Re(:,m,n)
      tmp%Im(:,mm,n) = c%Im(:,m,n)
_loop_mn_end
   call var_mpt_copy(tmp,c)
   end subroutine var_quadX

   subroutine var_printmpt(c)
      type (mpt), intent(in) :: c
      integer :: n
         do n = 0, 4!var_N%pH1
            write(*,'(2I3,10F9.3)') mpi_rnk, n+var_N%pH0,c%Re(1,0:4,n),c%Im(1,0:4,n)
!            print*, mpi_rnk, c%Re(1,:,m)
         end do
   end subroutine var_printmpt

   subroutine var_printmpt_all(c)
      type (mpt), intent(in) :: c
      integer :: n
         do n = 0, var_N%pH1
            write(*,'(2I3,140F9.3)') mpi_rnk, n+var_N%pH0,c%Re(1,:,n)!,c%Im(1,:,n)
!            print*, mpi_rnk, c%Re(1,:,m)
         end do
   end subroutine var_printmpt_all


   subroutine var_printphys(c)
      type (phys), intent(in) :: c
      integer :: m
         do m = 0, 4! var_M%pH1
            write(*,'(2I3,30F7.3)') mpi_rnk, m+var_M%pH0,c%Re(1,0:4,m)
!            print*, mpi_rnk, c%Re(1,:,m)
         end do
   end subroutine var_printphys

   subroutine var_printtran(c)
      type (tran), intent(in) :: c
      integer :: m
         do m = 0, 4! var_M%pH1
            write(*,'(2I3,30F7.3)') mpi_rnk,m+var_M%pH0, c%Re(1,0:4,m)
!            print*, mpi_rnk, c%Re(1,:,m)
         end do
   end subroutine var_printtran

   subroutine var_printspec(c)
      type (spec), intent(in) :: c
      integer :: n
         do n = 0,  4!var_N%pH1
            write(*,'(2I3,30F7.3)') mpi_rnk,n+var_N%pH0 ,c%Re(1,0:4,n),c%Im(1,0:4,n)
!            print*, mpi_rnk, c%Re(1,:,n)
         end do
   end subroutine var_printspec  

   subroutine var_randspec(c)
     type (spec), intent(out) :: c
     integer,parameter :: seed = 86456
     double precision :: r
     _loop_kmn_vars
     !         call srand(seed+mpi_rnk)
     _loop_kmn_begin
     !         if ((m+var_M%pH0)/=1) cycle
     if ((n+var_N%pH0)==0) cycle
     if (nn==i_NN1) cycle
     if (m==0) cycle
     call RANDOM_NUMBER(r)
     c%Re(k,m,n)=r!RAND()
     if (m==0) cycle
     if (n+var_N%pH0==0) cycle
     call RANDOM_NUMBER(r)
     c%Im(k,m,n)=r!RAND()
     _loop_kmn_end
   end subroutine var_randspec

   subroutine var_randadd(c,scl)
     type (mpt), intent(inout) :: c
     double precision, intent(in) :: scl
     integer,parameter :: seed = 86456
     double precision :: r
     _loop_kmn_vars
     _loop_kmn_begin
     if ((n+var_N%pH0)==0) cycle
     if (m==0) cycle
     call RANDOM_NUMBER(r)
     c%Re(k,m,n)=c%Re(k,m,n)+scl*r
     if (m==0) cycle
     if (n+var_N%pH0==0) cycle
     call RANDOM_NUMBER(r)
     c%Im(k,m,n)=c%Im(k,m,n)+scl*r
     _loop_kmn_end
   end subroutine var_randadd
   
   subroutine var_randmpt(c)
     type (mpt), intent(out) :: c
!      integer,parameter :: seed = 82342
      double precision :: r
      _loop_kmn_vars
      _loop_kmn_begin
      if ((n+var_N%pH0)==0) cycle
      if (m==0) cycle
      call RANDOM_NUMBER(r)
      c%Re(k,m,n)=r!RAND()
      if (m==0) cycle
      if (n+var_N%pH0==0) cycle
      call RANDOM_NUMBER(r)
      c%Im(k,m,n)=r!RAND()
      _loop_kmn_end
    end subroutine var_randmpt

    subroutine var_randphasempt(c)
     type (mpt), intent(out) :: c
!      integer,parameter :: seed = 82342
      double precision :: r
      integer :: state_size
      integer, allocatable, dimension(:) :: state
      _loop_kmn_vars
      
      !! Calls a new random seed
      call random_seed(size=state_size)
      allocate(state(state_size))
      state=i_rand_seed
      call RANDOM_SEED(put=state)
      
      _loop_kmn_begin
      mm = m
      if (mm > i_MM1) mm = m - i_M 
      if ((abs(mm).lt.i_kICx).and.(abs(nn).lt.i_kICz)) then
          call RANDOM_NUMBER(r)
          c%Re(k,m,n)=cos(2*d_PI*r)
          c%Im(k,m,n)=sin(2*d_Pi*r)
       else
          c%Re(k,m,n)=0d0
          c%Im(k,m,n)=0d0
       end if
       if (nn.eq.0) then
           c%Re(k,m,n)=0d0
           c%Im(k,m,n)=0d0
       end if
       if (mm.eq.0) then
           c%Re(k,m,n)=0d0
           c%Im(k,m,n)=0d0
       end if
      _loop_kmn_end
    end subroutine var_randphasempt
    
    subroutine var_simplempt(c)
      type (mpt), intent(out) :: c
      if (var_N%pH0==0) then
         c%Re(:,1,1)=1d0
         c%Im(:,2,2)=-1d0
      end if
    end subroutine var_simplempt
    
!-------------------------------------------------------------------------
!  initialise variable
!-------------------------------------------------------------------------
   subroutine var_mpt_init(a)
      type (mpt), intent(out) :: a
      a%Re = 0d0
      a%Im = 0d0
   end subroutine var_mpt_init

   subroutine var_spec_init(a)
      type (spec), intent(out) :: a
      a%Re = 0d0
      a%Im = 0d0
   end subroutine var_spec_init

   subroutine var_phys_init(a)
      type (phys), intent(out) :: a
      a%Re = 0d0
   end subroutine var_phys_init
!-------------------------------------------------------------------------
!  Copy a mptocated variable
!-------------------------------------------------------------------------
   subroutine var_mpt_copy(in, out)
      type (mpt), intent(in)  :: in
      type (mpt), intent(out) :: out
      out%Re(:,:,0:var_N%pH1) = in%Re(:,:,0:var_N%pH1)
      out%Im(:,:,0:var_N%pH1) = in%Im(:,:,0:var_N%pH1)
   end subroutine var_mpt_copy

   subroutine var_mpt_nonzero(in,in2)
      type (mpt), intent(in)  :: in,in2
     _loop_kmn_vars
     _loop_kmn_begin
     if ( k==1 .and. m==0 .and. abs(in%Re(k,m,n)-in2%Re(k,m,n)) > 1d-5) print*, mpi_rnk,n,in%Re(k,m,n),in2%Re(k,m,n)
     _loop_kmn_end
   end subroutine var_mpt_nonzero


!------------------------------------------------------------------------
!     out := out + in
!------------------------------------------------------------------------
   subroutine var_mpt_add(ac, a)
      type (mpt), intent(in)    :: ac
      type (mpt), intent(inout) :: a
      a%Re(:,:,0:var_N%pH1) = a%Re(:,:,0:var_N%pH1) + ac%Re(:,:,0:var_N%pH1)
      a%Im(:,:,0:var_N%pH1) = a%Im(:,:,0:var_N%pH1) + ac%Im(:,:,0:var_N%pH1)
   end subroutine var_mpt_add


!------------------------------------------------------------------------
!     out := out - in
!------------------------------------------------------------------------
   subroutine var_mpt_sub(ac, a)
      type (mpt), intent(in)    :: ac
      type (mpt), intent(inout) :: a
      a%Re(:,:,0:var_N%pH1) = a%Re(:,:,0:var_N%pH1) - ac%Re(:,:,0:var_N%pH1)
      a%Im(:,:,0:var_N%pH1) = a%Im(:,:,0:var_N%pH1) - ac%Im(:,:,0:var_N%pH1)
   end subroutine var_mpt_sub

!-------------------------------------------------------------------------
!  Copy a specocated variable
!-------------------------------------------------------------------------
   subroutine var_spec_copy(in, out)
      type (spec), intent(in)  :: in
      type (spec), intent(out) :: out
      out%Re(:,:,0:var_N%pH1) = in%Re(:,:,0:var_N%pH1)
      out%Im(:,:,0:var_N%pH1) = in%Im(:,:,0:var_N%pH1)
   end subroutine var_spec_copy


!------------------------------------------------------------------------
!     out := out + in
!------------------------------------------------------------------------
   subroutine var_spec_add(ac, a)
      type (spec), intent(in)    :: ac
      type (spec), intent(inout) :: a
      a%Re(:,:,0:var_N%pH1) = a%Re(:,:,0:var_N%pH1) + ac%Re(:,:,0:var_N%pH1)
      a%Im(:,:,0:var_N%pH1) = a%Im(:,:,0:var_N%pH1) + ac%Im(:,:,0:var_N%pH1)
   end subroutine var_spec_add


!------------------------------------------------------------------------
!     out := out - in
!------------------------------------------------------------------------
   subroutine var_spec_sub(ac, a)
      type (spec), intent(in)    :: ac
      type (spec), intent(inout) :: a
      a%Re(:,:,0:var_N%pH1) = a%Re(:,:,0:var_N%pH1) - ac%Re(:,:,0:var_N%pH1)
      a%Im(:,:,0:var_N%pH1) = a%Im(:,:,0:var_N%pH1) - ac%Im(:,:,0:var_N%pH1)
   end subroutine var_spec_sub

!-------------------------------------------------------------------------
! Convert mean-poloidal-toroidal to u,v,w representation
!-------------------------------------------------------------------------
   subroutine var_mpt2spec(mp,s)
     !ad_k1 designed to act on phi/v not psi/u/w
     ! u = M - dz psi + dxdy phi
     ! v = - (dxdx + dzdz) phi
     ! w = M + dx psi + dydz phi
     type(mpt), intent(in) :: mp
     type(spec), intent(out) :: s
     integer :: k,k1,k2,k3
     _loop_mn_vars
     _loop_mn_begin
     do k=0,i_K0-1
        k1=k+1
        k2=k+i_K0
        k3=k+2*i_K0
        s%Re(k1,m,n)=-(-mp%Im(k1,m,n)*ad_n1(nn) &
             +mp%Im(k2,m,n)*ad_k1(k)*ad_m1(m))
        s%Im(k1,m,n)= (-mp%Re(k1,m,n)*ad_n1(nn) &
             +mp%Re(k2,m,n)*ad_k1(k)*ad_m1(m))
        
        s%Re(k2,m,n)=-mp%Re(k2,m,n)*(ad_m2(m)+ad_n2(nn))
        s%Im(k2,m,n)=-mp%Im(k2,m,n)*(ad_m2(m)+ad_n2(nn))
        
        s%Re(k3,m,n)=-(mp%Im(k1,m,n)*ad_m1(m) &
             +mp%Im(k2,m,n)*ad_k1(k)*ad_n1(nn))
        s%Im(k3,m,n)= (mp%Re(k1,m,n)*ad_m1(m) &
             +mp%Re(k2,m,n)*ad_k1(k)*ad_n1(nn))
     end do
     _loop_mn_end
     if (var_N%pH0 == 0) then
        !no flux
        s%Re(1,0,0)=0d0
        s%Im(1,0,0)=0d0
        s%Re(2*i_K0,0,0)=0d0
        s%Im(2*i_K0,0,0)=0d0
        do k=1,i_K0-1
           k1=k+1
           k2=k+i_K0
           k3=k+2*i_K0
           s%Re(k1,0,0)= mp%Re(k1,0,0)
           s%Im(k1,0,0)= mp%Im(k1,0,0) 
           s%Re(k3,0,0)= mp%Re(k2,0,0)
           s%Im(k3,0,0)= mp%Im(k2,0,0)            
        end do
     end if
   end subroutine var_mpt2spec
   
   !-------------------------------------------------------------------------
   ! Convert u,v,w representation to mean-poloidal-toroidal after NL
   !-------------------------------------------------------------------------
   
   subroutine var_spec2mpt(s,mp)
     !k1  = - curl b = dx bz - dz bx
     !k2  = - curl curl b = (dxx + dzz)by - dy (dx bx + dz bz)  
     !Mean flows contributions in x and z are kept
     !    #-dxdy(bx) #--=+ because ad_y1 designed to act on phi/v not psi/u/w
     
     type(mpt), intent(out) :: mp
     type(spec), intent(in) :: s
     integer :: k,k1,k2,k3
     _loop_mn_vars
     _loop_mn_begin
     do k=0,i_K0-1
        k1=k+1
        k2=k+i_K0
        k3=k+2*i_K0
        mp%Re(k2,m,n)=-( s%Im(k1,m,n)*ad_m1(m) *ad_k1(k) &
             + s%Im(k3,m,n)*ad_n1(nn)*ad_k1(k)) &
             + s%Re(k2,m,n)*(ad_m2(m)+ad_n2(nn))
        mp%Im(k2,m,n)= ( s%Re(k1,m,n)*ad_m1(m) *ad_k1(k) &
             + s%Re(k3,m,n)*ad_n1(nn)*ad_k1(k)) &
             + s%Im(k2,m,n)*(ad_m2(m)+ad_n2(nn))
        mp%Re(k1,m,n)= ( s%Im(k3,m,n)*ad_m1(m) &
           - s%Im(k1,m,n)*ad_n1(nn))
        mp%Im(k1,m,n)=-( s%Re(k3,m,n)*ad_m1(m) &
             - s%Re(k1,m,n)*ad_n1(nn))
     end do
     _loop_mn_end
     if (var_N%pH0 == 0) then 
        mp%Re(1,0,0)=0d0
        mp%Im(1,0,0)=0d0
        do k=1,i_K0-1
           k1=k+1
           k2=k+i_K0
           k3=k+2*i_K0
           mp%Re(k1,0,0)=-s%Re(k1,0,0)
           mp%Im(k1,0,0)=-s%Im(k1,0,0) 
           mp%Re(k2,0,0)=-s%Re(k3,0,0)
           mp%Im(k2,0,0)=-s%Im(k3,0,0)            
        end do
     end if
   end subroutine var_spec2mpt
   
!-------------------------------------------------------------------------
!  null function
!-------------------------------------------------------------------------
   subroutine var_null()
   end subroutine var_null

   subroutine var_precompute()
      integer :: m,n,r,k
      ! distribute modes across processors
      var_M%pH1_ = -1
      do n = 0, i_3M-1
         r = _Np - modulo(n,_Np) - 1
         var_M%pH1_(r) = var_M%pH1_(r) + 1
      end do
      var_M%pH0_(0) = 0
      do r = 1, _Np-1
         var_M%pH0_(r) = var_M%pH0_(r-1) + var_M%pH1_(r-1) + 1
      end do
      var_M%pH0 = var_M%pH0_(mpi_rnk)
      var_M%pH1 = var_M%pH1_(mpi_rnk)

      var_N%pH1_ = -1
      do n = 0, i_NN1
         r = _Np - modulo(n,_Np) - 1
         var_N%pH1_(r) = var_N%pH1_(r) + 1
      end do
      var_N%pH0_(0) = 0
      do r = 1, _Np-1
         var_N%pH0_(r) = var_N%pH0_(r-1) + var_N%pH1_(r-1) + 1
      end do
      var_N%pH0 = var_N%pH0_(mpi_rnk)
      var_N%pH1 = var_N%pH1_(mpi_rnk)
!      print*, mpi_rnk,var_M%pH0,var_M%pH1
!      print*, mpi_rnk,var_N%pH0,var_N%pH1
      do m = 0, i_MM1
         if (m > 0) then   
            ad_m1(i_M-m) = -d_alpha*m
            ad_m2(i_M-m) = -d_alpha*m*d_alpha*m
         end if
         ad_m1(m) =  d_alpha*m
         ad_m2(m) = -d_alpha*m*d_alpha*m
      end do
      do n = 0, i_NN1
         ad_n2(n) = -d_gamma*n * d_gamma*n
         ad_n1(n) =  d_gamma*n
      end do
      do k = 0 , i_K0-1
         ad_k1(k)=(-1d0)**k * k * d_beta
         ad_k2(k)=-k*k*d_beta*d_beta
      end do 


   end subroutine var_precompute

  subroutine old_var_LHSRHS(l,r) 
   type (mpt), intent(out) :: l,r
    double precision :: lap,hlap,tvl,ilap16
   _loop_kmn_vars
   _loop_k0mn_begin
   lap =(ad_m2(m)+ad_k2(k)+ad_n2(nn))
   hlap=(ad_m2(m)         +ad_n2(nn))
   ilap16=d_hypo/(ad_m2(m)+ad_k2(k)+ad_n2(nn))
   if (s_dragall .or. k1==1) then
      tvl=(1d0 - d_dt * (lap/d_Re + ilap16 - d_drag))
   else
      tvl=(1d0 - d_dt * (lap/d_Re + ilap16))
   end if
   if (m==0 .and. nn==0) then
      l%Re(k1,m,n)=tvl
      l%Re(k2,m,n)=tvl
      r%Re(k1,m,n)=1d0
      r%Re(k2,m,n)=1d0
      l%Re(1,m,n)=1d0
      r%Re(1,m,n)=0d0
   else
      l%Re(k1,m,n)=tvl * hlap
      l%Re(k2,m,n)=tvl * hlap * lap
      r%Re(k1,m,n)=hlap
      r%Re(k2,m,n)=hlap * lap
   end if
   _loop_kmn_end

   end subroutine old_var_LHSRHS

  subroutine var_LHSRHS(l,r) 
   type (mpt), intent(out) :: l,r
    double precision :: lap,hlap,tvl,ilap16
   _loop_kmn_vars
   _loop_k0mn_begin
   lap =(ad_m2(m)+ad_k2(k)+ad_n2(nn))
   hlap=(ad_m2(m)         +ad_n2(nn))
   ilap16=d_hypo/(ad_m2(m)+ad_k2(k)+ad_n2(nn))
   if (s_dragall .or. k1==1) then
      tvl=(1d0 - d_dt * (lap/d_Re + ilap16 - d_drag))
   else
      tvl=(1d0 - d_dt * (lap/d_Re + ilap16))
   end if
   if (m==0 .and. nn==0) then
      l%Re(k1,m,n)=tvl
      l%Re(k2,m,n)=tvl
      r%Re(k1,m,n)=1d0
      r%Re(k2,m,n)=1d0
      l%Re(1,m,n)=1d0
      r%Re(1,m,n)=0d0
   else
      l%Re(k1,m,n)=tvl * hlap
      if (s_dragall) then
         l%Re(k2,m,n)=tvl * hlap * lap + d_dt * (d_vdrag - d_drag) * hlap * hlap
      else
         l%Re(k2,m,n)=tvl * hlap * lap 
      end if
      r%Re(k1,m,n)=hlap
      r%Re(k2,m,n)=hlap * lap
   end if
   _loop_kmn_end

   end subroutine var_LHSRHS

   
   subroutine var_spec_grad(p, ux,uz)
     type (spec), intent(in)  :: p
     type (spec), intent(out) :: ux,uz
     _loop_mn_vars
     
     _loop_mn_begin
     ux%Re(:,m,n) = -p%Im(:,m,n)*ad_m1(m)
     ux%Im(:,m,n) =  p%Re(:,m,n)*ad_m1(m)
     uz%Re(:,m,n) = -p%Im(:,m,n)*ad_n1(nn)
     uz%Im(:,m,n) =  p%Re(:,m,n)*ad_n1(nn)
     _loop_mn_end
     
   end subroutine var_spec_grad
   
 end module variables
