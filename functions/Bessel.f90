module bessel
! Bessel functions and legendre polynomial implementation


contains

	SUBROUTINE STVH0(X,SH0)
	!       =============================================
	!       Purpose: Compute Struve function H0(x)
	!       Input :  x   --- Argument of H0(x) ( x ò 0 )
	!       Output:  SH0 --- H0(x)
	!       =============================================
			IMPLICIT DOUBLE PRECISION (A-H,O-Z)
			integer K,KM
			!real*8 SH0
			
			PIi=3.141592653589793D0
			S=1.0D0
			R=1.0D0
			IF (X.LE.20.0D0) THEN
			   A0=2.0*X/PIi
			   DO 10 K=1,60
				  R=-R*X/(2.0D0*K+1.0D0)*X/(2.0D0*K+1.0D0)
				  S=S+R
				  IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
	10         CONTINUE
	15         SH0=A0*S
			ELSE
			   KM=INT(.5*(X+1.0))
			   IF (X.GE.50.0) KM=25
			   DO 20 K=1,KM
				  R=-R*((2.0D0*K-1.0D0)/X)**2
				  S=S+R
				  IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 25
	20         CONTINUE
	25         T=4.0D0/X
			   T2=T*T
			   P0=((((-.37043D-5*T2+.173565D-4)*T2-.487613D-4)     &
				  *T2+.17343D-3)*T2-.1753062D-2)*T2+.3989422793D0
			   Q0=T*(((((.32312D-5*T2-.142078D-4)*T2+.342468D-4)*  &
				  T2-.869791D-4)*T2+.4564324D-3)*T2-.0124669441D0)
			   TA0=X-.25D0*PIi
			   BY0=2.0D0/DSQRT(X)*(P0*DSIN(TA0)+Q0*DCOS(TA0))
			   SH0=2.0D0/(PIi*X)*S+BY0
			ENDIF
			RETURN
	END

	FUNCTION BESSY (N,X)
	! ------------------------------------------------------------------
	!     This subroutine calculates the second kind Bessel Function of
	!     integer order N, for any real X. We use here the classical
	!     recursive formula. 
	! ------------------------------------------------------------------
		  IMPLICIT NONE
		  INTEGER N, J
		  REAL *8 X,TOX,BY,BYM,BYP,BESSY!,BESSY0,BESSY1,
		  IF (N.EQ.0) THEN
		  BESSY = BESSY0(X)
		  RETURN
		  ENDIF
		  IF (N.EQ.1) THEN
		  BESSY = BESSY1(X)
		  RETURN
		  ENDIF
		  IF (X.EQ.0.) THEN
		  BESSY = -1.E30
		  RETURN
		  ENDIF
		  TOX = 2./X
		  BY  = BESSY1(X)
		  BYM = BESSY0(X)
		  DO 11 J = 1,N-1
		  BYP = J*TOX*BY-BYM
		  BYM = BY
		  BY  = BYP
	   11 CONTINUE
		  BESSY = BY
		  RETURN
		  END

	FUNCTION BESSY0 (X)
		  IMPLICIT NONE
		  REAL *8 X,FS,FR,Z,FP,FQ,XX,BESSY0!,BESSJ0,
	! ---------------------------------------------------------------------
	!     This subroutine calculates the Second Kind Bessel Function of
	!     order 0, for any real number X. The polynomial approximation by
	!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	!     REFERENCES:
	!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	!     VOL.5, 1962.
	! ---------------------------------------------------------------------
		  REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
				   ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
		  DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
		  -.2073370639D-5,.2093887211D-6 /
		  DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3,  &
		  -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
		  DATA R1,R2,R3,R4,R5,R6 /-2957821389.D0,7062834065.D0, &
		  -512359803.6D0,10879881.29D0,-86327.92757D0,228.4622733D0 /
		  DATA S1,S2,S3,S4,S5,S6 /40076544269.D0,745249964.8D0, &
		  7189466.438D0,47447.26470D0,226.1030244D0,1.D0 /
		  IF (X.EQ.0.D0) THEN
		  BESSY0 = -1.E30
		  RETURN
		  ENDIF
		  IF (X.LT.8.D0) THEN
		  Y = X*X
		  FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
		  FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
		  BESSY0 = FR/FS+.636619772D0*BESSJ0(X)*LOG(X)
		  ELSE
		  Z = 8.D0/X
		  Y = Z*Z
		  XX = X-.785398164D0
		  FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
		  FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
		  BESSY0 = SQRT(.636619772D0/X)*(FP*SIN(XX)+Z*FQ*COS(XX))
		  ENDIF
		  RETURN
		  END
		  
	FUNCTION BESSY1 (X)
		  IMPLICIT NONE
		  REAL *8 X,FR,FS,Z,FP,FQ,XX,BESSY1!,BESSJ1
	! ----------------------------------------------------------------------
	!     This subroutine calculates the Second Kind Bessel Function of
	!     order 1, for any real number X. The polynomial approximation by
	!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	!     REFERENCES:
	!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	!     VOL.5, 1962.
	! ----------------------------------------------------------------------
		  REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
				   ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6,S7
		  DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4, &
		  .2457520174D-5,-.240337019D-6 /
		  DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,  &
		  .8449199096D-5,-.88228987D-6,.105787412D-6 /
		  DATA R1,R2,R3,R4,R5,R6 /-.4900604943D13,.1275274390D13,   &
		  -.5153438139D11,.7349264551D9,-.4237922726D7,.8511937935D4 /
		  DATA S1,S2,S3,S4,S5,S6,S7 /.2499580570D14,.4244419664D12, &
		  .3733650367D10,.2245904002D8,.1020426050D6,.3549632885D3,1.D0 /
		  IF (X.EQ.0.) THEN
		  BESSY1 = -1.E30
		  RETURN
		  ENDIF
		  IF (X.LT.8.) THEN
		  Y = X*X
		  FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
		  FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*(S6+Y*S7)))))
		  BESSY1 = X*(FR/FS)+.636619772*(BESSJ1(X)*LOG(X)-1./X)
		  ELSE
		  Z = 8./X
		  Y = Z*Z
		  XX = X-2.356194491
		  FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
		  FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
		  BESSY1 = SQRT(.636619772/X)*(SIN(XX)*FP+Z*COS(XX)*FQ)
		  ENDIF
		  RETURN
		  END
		  
	FUNCTION BESSJ0 (X)
		  IMPLICIT NONE
		  REAL *8 X,AX,FR,FS,Z,FP,FQ,XX,BESSJ0
	! --------------------------------------------------------------------
	!     This subroutine calculates the First Kind Bessel Function of
	!     order 0, for any real number X. The polynomial approximation by
	!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	!     REFERENCES:
	!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	!     VOL.5, 1962.
	! ---------------------------------------------------------------------
		  REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
				   ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
		  DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
		  -.2073370639D-5,.2093887211D-6 /
		  DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
		  -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
		  DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
		  651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
		  DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
		  9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
		  IF(X.EQ.0.D0) GO TO 1
		  AX = ABS (X)
		  IF (AX.LT.8.) THEN
		  Y = X*X
		  FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
		  FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
		  BESSJ0 = FR/FS
		  ELSE
		  Z = 8./AX
		  Y = Z*Z
		  XX = AX-.785398164
		  FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
		  FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
		  BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
		  ENDIF
		  RETURN
		1 BESSJ0 = 1.D0
		  RETURN
		  END

	FUNCTION BESSJ1 (X)
		  IMPLICIT NONE
		  REAL *8 X,AX,FR,FS,Z,FP,FQ,XX,BESSJ1
	! ---------------------------------------------------------------------
	!     This subroutine calculates the First Kind Bessel Function of
	!     order 1, for any real number X. The polynomial approximation by
	!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	!     REFERENCES:
	!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	!     VOL.5, 1962.
	! ---------------------------------------------------------------------
		  REAL *8 Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
				   ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
		  DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
		  .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
		  DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
		  .8449199096D-5,-.88228987D-6,.105787412D-6 /
		  DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, & 
		  242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
		  DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
		  18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

		  AX = ABS(X)
		  IF (AX.LT.8.) THEN
		  Y = X*X
		  FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
		  FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
		  BESSJ1 = X*(FR/FS)
		  ELSE
		  Z = 8./AX
		  Y = Z*Z
		  XX = AX-2.35619491
		  FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
		  FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
		  BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
		  ENDIF
		  RETURN
		  END
			  
	subroutine pmns_polynomial_value ( mm, n, m, x, cx )

	!*****************************************************************************80
	!
	!! PMNS_POLYNOMIAL_VALUE: sphere-normalized Legendre polynomial Pmns(n,m,x).
	!
	!  Discussion:
	!
	!    The unnormalized associated Legendre functions P_N^M(X) have
	!    the property that
	!
	!      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX 
	!      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
	!
	!    By dividing the function by the square root of this term,
	!    the normalized associated Legendre functions have norm 1.
	!
	!    However, we plan to use these functions to build spherical
	!    harmonics, so we use a slightly different normalization factor of
	!
	!      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) ) 
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license. 
	!
	!  Modified:
	!
	!    03 May 2013
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Reference:
	!
	!    Milton Abramowitz, Irene Stegun,
	!    Handbook of Mathematical Functions,
	!    National Bureau of Standards, 1964,
	!    ISBN: 0-486-61272-4,
	!    LC: QA47.A34.
	!
	!  Parameters:
	!
	!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
	!
	!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
	!    function, which must be at least 0.
	!
	!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
	!    which must be at least 0, and no greater than N.
	!
	!    Input, real ( kind = 8 ) X(MM), the evaluation points.
	!
	!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
	!
	  implicit none

	  integer ( kind = 4 ) mm
	  integer ( kind = 4 ) n

	  real ( kind = 8 ) cx(mm,0:n)
	  real ( kind = 8 ) factor
	  integer ( kind = 4 ) j
	  integer ( kind = 4 ) m
	  !real ( kind = 8 ) r8_factorial
	  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
	  real ( kind = 8 ) x(mm)

	  cx(1:mm,0:n) = 0.0D+00

	  if ( m <= n ) then
		cx(1:mm,m) = 1.0D+00
		factor = 1.0D+00
		do j = 1, m
		  cx(1:mm,m) = - cx(1:mm,m) * factor * sqrt ( 1.0D+00 - x(1:mm) ** 2 )
		  factor = factor + 2.0D+00
		end do
	  end if

	  if ( m + 1 <= n ) then
		cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
	  end if

	  do j = m + 2, n
		cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
					 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &
					 / real (     j - m,     kind = 8 )
	  end do
	!
	!  Normalization.
	!
	  do j = m, n
		factor = sqrt ( ( real ( 2 * j + 1, kind = 8 ) * r8_factorial ( j - m ) ) &
		  / ( 4.0D+00 * r8_pi * r8_factorial ( j + m ) ) )
		cx(1:mm,j) = cx(1:mm,j) * factor
	  end do

	  return
	end 
	   
	function r8_factorial ( n ) result(r8)
	implicit none

	  real ( kind = 8 ) r8
	  integer ( kind = 4 ) i
	  integer ( kind = 4 ),intent(in):: n

	  r8= 1.0D+00

	  do i = 1, n
		r8 = r8 * real ( i, kind = 8 )
	  end do

	  return
	end  
		

end module bessel