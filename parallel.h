#ifndef parallel_h
#define parallel_h
!//***************************************************************/
#define _Np 128
#define _CPUTIME
!/***************************************************************/

!/*-------------------------------------------------------------*/
#if _Np != 1
#define _MPI
#endif

!/*-------------------------------------------------------------*/
#ifndef _MPI
!/*-------------------------------------------------------------*/

#define _loop_mn_vars 	\
   integer :: m,n,mm,nn

#define _loop_kmn_vars   \
   integer :: m,mm,n,nn,k,k1,k2,k3

#define _loop_mn_begin	\
   do n = 0, i_NN1;	\
      nn=n; \
      do m = 0, i_M1		

#define _loop_phy_begin  \
   do m = 0, i_3M-1;  \
      mm=m; \
      do n = 0, i_3N-1    

#define _loop_kmn_begin  \
   do n = 0, i_NN1 ;   \
      nn = n;    \
      do m = 0,i_M1;     \
         do k = 1,i_K

#define _loop_k0mn_begin  \
  do n = 0, i_NN1;\
   nn = n;\
      do m = 0,i_M1;\
         do k = 0,i_K0-1;\
            k1 = k + 1;\
            k2 = k + i_K0 

#define _loop_kmn_end \
         end do;     \
      end do;     \
   end do

#define _loop_mn_end	\
      end do;		\
   end do

!/*-------------------------------------------------------------*/
#else  /* def _MPI */
!/*-------------------------------------------------------------*/

#define _loop_mn_vars	\
   integer :: pH0,pH1,m,n,mm,nn

#define _loop_mn_begin	\
   pH0=var_N%pH0 ; 	\
   pH1=var_N%pH1 ;	\
   do n = 0, pH1 ;	\
      nn = n + pH0;    \
      do m = 0,i_M1

#define _loop_phy_begin  \
   pH0=var_M%pH0 ;   \
   pH1=var_M%pH1 ;   \
   do m = 0, pH1 ;   \
      mm = m + pH0;    \
      do n = 0,i_3N-1

#define _loop_mn_end	\
      end do;		\
   end do

#define _loop_kmn_vars   \
   integer :: pH0,pH1,m,mm,n,nn,k,k1,k2,k3

#define _loop_kmn_begin  \
   pH0=var_N%pH0 ;   \
   pH1=var_N%pH1 ;   \
   do n = 0, pH1 ;   \
      nn = n + pH0;    \
      do m = 0,i_M1;     \
         do k = 1,i_K

#define _loop_k0mn_begin  \
   pH0=var_N%pH0 ;   \
   pH1=var_N%pH1 ;   \
   do n = 0, pH1 ;   \
      nn = n + pH0;    \
      do m = 0,i_M1;     \
         do k = 0,i_K0-1;     \
            k1 = k + 1;     \
            k2 = k + i_K0 


#define _loop_kmn_end \
         end do;     \
      end do;     \
   end do

!/*-------------------------------------------------------------*/
#endif /* def _MPI */
!/*-------------------------------------------------------------*/

#endif /* ndef parallel_h */
