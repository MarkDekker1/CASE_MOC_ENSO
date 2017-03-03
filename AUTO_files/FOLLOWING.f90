!----------------------------------------------------------------------
!   Model subroutines
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----


      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION Tr0,alpha,H,hstar,eps,L,r,z0,mu,zeta, kwant,tau_ext,Tr,T1,T2,h1,Tsub,tau, Hm, b, beta, uv, wv,h2,k,Q

	tau_ext=PAR(1)
	eps = PAR(2)
	zeta = PAR(3)
	mu = PAR(4)

	Tr0 = 16.
	alpha = 1./180.
	H = 100.
	hstar = 62.
	L = 15e6
	r = 1./400.
	z0 = 75.
	!mu = 0.0026
	Q = 22.
	beta=mu/0.02
	b = 22.*beta/(mu*L)
	Tr = 29.5
	Hm = 50.

	T1 = U(1)
	T2 = U(2)
	h1 = U(3)

	tau = tau_ext+mu/beta*(T2-T1)
	uv = eps*beta*tau*L/2.
	wv = -zeta*beta*tau*Hm
	h2 = h1+b*L*tau

	Tsub = Tr-(Tr-Tr0)/2.*(1.-tanh((H-z0+h2)/hstar))

	F(1) = -alpha*(T1-Tr)-(T2-T1)*uv/(L/2.)
	F(2) = -alpha*(T2-Tr)-(T2-Tsub)*wv/Hm
	F(3) = r*(-h1-b*L*tau/2.)

      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,Z)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: Z

!      Set the parameters
       PAR(1)=0.	!tau_ext
       PAR(2)=0.05 	!eps
       PAR(3)=1.3	!zeta
       PAR(4)=0.0026	!mu

!      Set the variables equilibria
       U(1)= 27.807842171280281
       U(2)= 19.304052596936508
       U(3)= 93.541685317706509

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
