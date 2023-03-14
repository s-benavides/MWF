!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
module modes
  !*************************************************************************
  use parameters
  implicit none
  save
  
contains

  function udotgradu(v,vx,vz) result(ans)
    double precision :: ans(i_KK)
    double precision,intent(in) :: v(i_KK),vx(i_KK),vz(i_KK)
    ans(1) =v(1)*vx(1) + (v(2)*vx(2))/2 + (v(3)*vx(3))/2 + (v(4)*vx(4))/2 + v(8)*vz(1) + (v(9)*vz(2))/2 + (v(10)*vz(3))/2 + (v(11)*vz(4))/2 + (d_beta*v(2)*v(5))/2 - d_beta*v(3)*v(6) + (3*d_beta*v(4)*v(7))/2
    ans(2) =v(1)*vx(2) + v(2)*vx(1) - (v(2)*vx(3))/2 - (v(3)*vx(2))/2 + (v(3)*vx(4))/2 + (v(4)*vx(3))/2 + v(8)*vz(2) + v(9)*vz(1) - (v(9)*vz(3))/2 - (v(10)*vz(2))/2 + (v(10)*vz(4))/2 + (v(11)*vz(3))/2 + (d_beta*v(2)*v(6))/2 - d_beta*v(3)*v(5) + d_beta*v(3)*v(7) - (3*d_beta*v(4)*v(6))/2
    ans(3) =v(1)*vx(3) - (v(2)*vx(2))/2 + v(3)*vx(1) + (v(2)*vx(4))/2 + (v(4)*vx(2))/2 + v(8)*vz(3) - (v(9)*vz(2))/2 + v(10)*vz(1) + (v(9)*vz(4))/2 + (v(11)*vz(2))/2 + (d_beta*v(2)*v(5))/2 + (d_beta*v(2)*v(7))/2 + (3*d_beta*v(4)*v(5))/2
    ans(4) =v(1)*vx(4) + (v(2)*vx(3))/2 + (v(3)*vx(2))/2 + v(4)*vx(1) + v(8)*vz(4) + (v(9)*vz(3))/2 + (v(10)*vz(2))/2 + v(11)*vz(1) + (d_beta*v(2)*v(6))/2 - d_beta*v(3)*v(5)
    ans(5) =v(1)*vx(5) + (v(2)*vx(6))/2 + (v(3)*vx(5))/2 + (v(3)*vx(7))/2 + (v(4)*vx(6))/2 + v(8)*vz(5) + (v(9)*vz(6))/2 + (v(10)*vz(5))/2 + (v(10)*vz(7))/2 + (v(11)*vz(6))/2 + (d_beta*v(5)*v(6))/2 - (d_beta*v(6)*v(7))/2
    ans(6) =v(1)*vx(6) + (v(2)*vx(5))/2 - (v(2)*vx(7))/2 + (v(4)*vx(5))/2 + v(8)*vz(6) + (v(9)*vz(5))/2 - (v(9)*vz(7))/2 + (v(11)*vz(5))/2 - (d_beta*v(5)*v(5))/2 - d_beta*v(5)*v(7)
    ans(7) =v(1)*vx(7) - (v(2)*vx(6))/2 + (v(3)*vx(5))/2 + v(8)*vz(7) - (v(9)*vz(6))/2 + (v(10)*vz(5))/2 + (3*d_beta*v(5)*v(6))/2
    ans(8) =v(1)*vx(8) + (v(2)*vx(9))/2 + (v(3)*vx(10))/2 + (v(4)*vx(11))/2 + v(8)*vz(8) + (v(9)*vz(9))/2 + (v(10)*vz(10))/2 + (v(11)*vz(11))/2 + (d_beta*v(5)*v(9))/2 - d_beta*v(6)*v(10) + (3*d_beta*v(7)*v(11))/2
    ans(9) =v(1)*vx(9) + v(2)*vx(8) - (v(2)*vx(10))/2 - (v(3)*vx(9))/2 + (v(3)*vx(11))/2 + (v(4)*vx(10))/2 + v(8)*vz(9) + v(9)*vz(8) - (v(9)*vz(10))/2 - (v(10)*vz(9))/2 + (v(10)*vz(11))/2 + (v(11)*vz(10))/2 - d_beta*v(5)*v(10) + (d_beta*v(6)*v(9))/2 - (3*d_beta*v(6)*v(11))/2 + d_beta*v(7)*v(10)
    ans(10) =v(1)*vx(10) - (v(2)*vx(9))/2 + v(3)*vx(8) + (v(2)*vx(11))/2 + (v(4)*vx(9))/2 + v(8)*vz(10) - (v(9)*vz(9))/2 + v(10)*vz(8) + (v(9)*vz(11))/2 + (v(11)*vz(9))/2 + (d_beta*v(5)*v(9))/2 + (3*d_beta*v(5)*v(11))/2 + (d_beta*v(7)*v(9))/2
    ans(11) =v(1)*vx(11) + (v(2)*vx(10))/2 + (v(3)*vx(9))/2 + v(4)*vx(8) + v(8)*vz(11) + (v(9)*vz(10))/2 + (v(10)*vz(9))/2 + v(11)*vz(8) - d_beta*v(5)*v(10) + (d_beta*v(6)*v(9))/2
  end function udotgradu

  function udotgradu2(u,ux,uz) result(ans)
    double precision :: ans(i_KK)
    double precision,intent(in) :: u(i_KK),ux(i_KK),uz(i_KK)
    ans(1:i_K0)= nluw(u(1:i_K0),ux(1:i_K0)) + nluyv(u(1:i_K0),u(i_K0+1:2*i_K0-1)) + nluw(u(2*i_K0:i_KK),uz(1:i_K0))
    ans(i_K0+1:2*i_K0-1)= nluv(u(1:i_K0),ux(i_K0+1:2*i_K0-1)) + nlvvy(u(i_K0+1:2*i_K0-1))        + nluv(u(2*i_K0:i_KK),uz(i_K0+1:2*i_K0-1))
    ans(2*i_K0:i_KK)=nluw(u(1:i_K0),ux(2*i_K0:i_KK))+ nluyv(u(2*i_K0:i_KK),u(i_K0+1:2*i_K0-1))+ nluw(u(2*i_K0:i_KK),uz(2*i_K0:i_KK))
  end function udotgradu2
  
  function nluw(u,w) result(ans)
    double precision :: ans(0:i_K1)
    double precision,intent(in) :: u(0:i_K1),w(0:i_K1)
    ans(0)=u(0)*w(0) + (u(1)*w(1))/2 + (u(2)*w(2))/2 + (u(3)*w(3))/2
    ans(1)=u(0)*w(1) + u(1)*w(0) - (u(1)*w(2))/2 - (u(2)*w(1))/2 + (u(2)*w(3))/2 + (u(3)*w(2))/2	
    ans(2)=u(0)*w(2) - (u(1)*w(1))/2 + u(2)*w(0) + (u(1)*w(3))/2 + (u(3)*w(1))/2
    ans(3)=u(0)*w(3) + (u(1)*w(2))/2 + (u(2)*w(1))/2 + u(3)*w(0)
    RETURN
  end function nluw

  function nluyv(u,v) result(ans)
    double precision :: ans(0:i_K1)
    double precision,intent(in) :: u(0:i_K1),v(1:i_K1)
    ans(0)=(d_beta*u(1)*v(1))/2 - d_beta*u(2)*v(2) + (3*d_beta*u(3)*v(3))/2
    ans(1)=(d_beta*u(1)*v(2))/2 - d_beta*u(2)*v(1) + d_beta*u(2)*v(3) - (3*d_beta*u(3)*v(2))/2
    ans(2)=(d_beta*u(1)*v(1))/2 + (d_beta*u(1)*v(3))/2 + (3*d_beta*u(3)*v(1))/2
    ans(3)=(d_beta*u(1)*v(2))/2 - d_beta*u(2)*v(1)
    RETURN
  end function nluyv
  
  function nluv(u,v) result(ans)
    double precision :: ans(1:i_K1)
    double precision,intent(in) :: u(0:i_K1),v(1:i_K1)
    ans(1)=u(0)*v(1) + (u(1)*v(2))/2 + (u(2)*v(1))/2 + (u(2)*v(3))/2 + (u(3)*v(2))/2
    ans(2)=u(0)*v(2) + (u(1)*v(1))/2 - (u(1)*v(3))/2 + (u(3)*v(1))/2
    ans(3)=u(0)*v(3) - (u(1)*v(2))/2 + (u(2)*v(1))/2
  end function nluv
  
  function nlvvy(v) result(ans)
    double precision :: ans(1:i_K1)
    double precision,intent(in) :: v(1:i_K1)
    ans(1)=(d_beta*v(1)*v(2))/2 - (d_beta*v(2)*v(3))/2
    ans(2)=- (d_beta*v(1)**2)/2 - d_beta*v(1)*v(3)
    ans(3)=(3*d_beta*v(1)*v(2))/2
  end function nlvvy
  
  function nlvv(v) result(ans)
    double precision :: ans(0:i_K1)
    double precision,intent(in) :: v(1:i_K1)
    ans(0)=v(1)**2/2 + v(2)**2/2 + v(3)**2/2
    ans(1)=v(1)*v(2) - v(2)*v(3)
    ans(2)=v(1)**2/2 + v(3)*v(1)
    ans(3)=v(1)*v(2)
  end function nlvv
  
  function nlvv2(v1,v2) result(ans)
    double precision :: ans(0:i_K1)
    double precision,intent(in) :: v1(1:i_K1),v2(1:i_K1)
    ans(0)=v1(1)*v2(1)/2 + v1(2)*v2(2)/2 + v1(3)*v2(3)/2
    ans(1)=v2(2)*v1(1)/2 - v2(2)*v1(3)/2 + v1(2)*v2(1)/2 - v1(2)*v2(3)/2
    ans(2)=v1(1)*v2(1)/2 + v1(1)*v2(3)/2 + v1(3)*v2(1)/2
    ans(3)=v1(1)*v2(2)/2 + v2(1)*v1(2)/2
  end function nlvv2
  
  double precision function velU(V,y) 
    double precision, intent(in) :: V(i_KK),y
    velU= V(1) &
         +dsin(d_beta*y) * V(2) &
         +dcos(2d0*d_beta*y) * V(3)&
         +dsin(3d0*d_beta*y) * V(4)
    return
  end function velU
  
  double precision function velW(V,y)
    double precision, intent(in) :: V(i_KK),y
    velW = V(8)&
         +dsin(d_beta*y) * V(9)&
         +dcos(2d0*d_beta*y) * V(10)&
         +dsin(3d0*d_beta*y) * V(11)
    return
  end function velW

  double precision function velV(V,y)
    double precision, intent(in) :: V(i_KK),y
    velV= dcos(d_beta*y) * V(5)&
         +dsin(2d0*d_beta*y) * V(6)&
         +dcos(3d0*d_beta*y) * V(7)
    return
  end function velV

end module modes

