!----------------------------------------------------------------------
!   Model subroutines
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----


      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION tr, H, td, ta, qtje, V, alphat, alphas, S0, Th, rho0, Fs, DT, DS, Drho, Q

	tr = 250
	H = 4500
	td = 180*365
	ta = 29*365
	qtje = 192e10*3600*24
	V = 300*4.5*8200*1e9
	alphat = 1e-4
	alphas = 76e-5
	S0 = 35
	!Th = 25
	rho0 = 1

	Fs = PAR(1)
	Th = PAR(2)

	DT = U(1)
	DS = U(2)

	Drho = rho0*(alphat*DT-alphas*DS)
	Q = 1/td+qtje*Drho**2/(V)

	F(1) = -1/tr*(DT-Th)-Q*DT
	F(2) = Fs*S0/H-Q*DS

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
       PAR(1)=0.	!Fs
       PAR(2)=0.	!Th


!      Set the variables equilibria
       U(1)=0.
       U(2)=0.

      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
