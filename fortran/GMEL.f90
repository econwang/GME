! python c:\python27\scripts\f2py.py -c -m fgmel gmel.f90 --fcompiler=gnu95 -lnlopt -lm
MODULE DATA
  IMPLICIT NONE
  REAL*8, ALLOCATABLE :: Y(:), X(:,:), Z(:,:), V(:)
END MODULE DATA

SUBROUTINE ESTIMATE(BETA,VCOV,H,Y1,X1,Z1,V1,T,K,M,J)
  USE DATA
  IMPLICIT NONE
  integer T, K, M, J
  REAL*8, INTENT(OUT), DIMENSION(K) :: BETA
  REAL*8, INTENT(OUT), DIMENSION(K,K) :: VCOV
  REAL*8, INTENT(OUT) :: H
  REAL*8, INTENT(IN), DIMENSION(T) :: Y1
  REAL*8, INTENT(IN), DIMENSION(T,K) :: X1
  REAL*8, INTENT(IN), DIMENSION(K,M) :: Z1
  REAL*8, INTENT(IN), DIMENSION(J) :: V1
  INTEGER*8 opt
  REAL*8 LB(T), UB(T), MAXH, lamest(T), DELTAS, OMES
  REAL*8 P(K,M), W(T,J), OMEGA(K,M), OMEGAS(K), LAM0(T), XOUT(K,K), XOUTINV(K,K)
  INTEGER IRES
  EXTERNAL DUAL_OBJ
  INCLUDE 'nlopt.f'
  Y = Y1
  X = X1
  Z = Z1
  V = V1
  lb = -1.0D5
  ub = 1.0D5
  LAM0 = 0.0D0
  CALL NLO_CREATE(OPT,NLOPT_LD_TNEWTON,T)
  CALL nlo_set_lower_bounds(ires, opt, lb)
  call nlo_set_upper_bounds(ires, opt, ub)
  call nlo_set_vector_storage(ires, opt, 200)
  call nlo_set_min_objective(ires, opt, DUAL_OBJ, 0)
  call nlo_set_xtol_rel(ires, opt, 1.0D-10)
  lamest = lam0
  call nlo_optimize(ires, opt, LAMEST, MAXH)
  call nlo_destroy(opt)
  OMEGA=EXP(-Z*SPREAD(MATMUL(TRANSPOSE(X),LAMEST),2,M))
  OMEGAS = SUM(OMEGA,2)
  P = OMEGA/SPREAD(OMEGAS,2,M)
  BETA = SUM(Z*P,2)
  H = MAXH
  W = EXP(-1*MATMUL(RESHAPE(LAMEST,(/T,1/)),RESHAPE(V,(/1,J/))))/ &
      spread(sum(EXP(-1*MATMUL(RESHAPE(LAMEST,(/T,1/)),RESHAPE(V,(/1,J/)))),2),2,J)
  DELTAS = SUM(LAMEST**2)/T
  OMES = (SUM(1/(sum((spread(v,1,T)**2)*W,2)-(sum(spread(v,1,T)*W,2))**2))/T)**2
  XOUT = MATMUL(TRANSPOSE(X),X)
  CALL inverse(XOUT,XOUTINV,K)
  VCOV = (DELTAS/OMES)*XOUTINV
  RETURN
end subroutine estimate

SUBROUTINE DUAL_OBJ(OBJ, T, LAM, GRAD, NEED_GRADIENT)
  USE DATA
  IMPLICIT NONE
  REAL*8 OBJ, LAM(T), GRAD(T)
  INTEGER T, NEED_GRADIENT, M, J
  M = SIZE(Z,2)
  J = SIZE(V)

  if (need_gradient.ne.0) then
  grad = y - matmul(x,SUM(Z*EXP(-Z*SPREAD(MATMUL(TRANSPOSE(X),LAM),2,M)) &
         /SPREAD(SUM(EXP(-Z*SPREAD(MATMUL(TRANSPOSE(X),LAM),2,M)),2),2,M),2)) - &
         sum(spread(v,1,T)*EXP(-1*MATMUL(RESHAPE(LAM,(/T,1/)),RESHAPE(V,(/1,J/))))/ &
         spread(sum(EXP(-1*MATMUL(RESHAPE(LAM,(/T,1/)),RESHAPE(V,(/1,J/)))),2),2,J),2)
  endif

  OBJ = DOT_PRODUCT(Y,LAM) + &
        SUM(LOG(SUM(EXP(-Z*SPREAD(MATMUL(TRANSPOSE(X),LAM),2,M)),2))) + &
        SUM(LOG(SUM(EXP(-1*MATMUL(RESHAPE(LAM,(/T,1/)),RESHAPE(V,(/1,J/)))),2)))
  RETURN
END SUBROUTINE DUAL_OBJ


 subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      enddo
   enddo
enddo

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
enddo
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  enddo
enddo

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    enddo
  enddo
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    enddo
    x(i) = x(i)/u(i,i)
  enddo
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
enddo
end subroutine inverse