SUBROUTINE tridiag(a,b,c,r,u,N)
!  solves a tridiagonal linear system  Au = r

!  algorithm comes from Numerical Recipes book

!  a,b,c are vectors of length N giving the non-zero elements of matrix A
!  a(2:N) are the entries to the left of the main diagonal (a(1) is ignored)
!  b(1:N) are main diagonal elements
!  c(1:N-1) are to the right of main diagonal (c(N) is ignored)
!  r is a vector of length N giving RHS of matrix equation
!  u is the unknown vector we are solving for
	INTEGER N,j
	REAL*8, DIMENSION(N) :: a,b,c,r,u
	REAL*8, DIMENSION(N) :: gam
	REAL*8 bet
	
	bet = b(1)
	IF (bet == 0.) THEN
		PRINT *, 'First diagonal element is zero'
		STOP
	END IF
	u(1) = r(1)/bet
	DO j=2,N
		gam(j) = c(j-1)/bet
		bet = b(j)-a(j)*gam(j)
		IF (bet == 0.) THEN
			PRINT *, 'tridiag trying to divide by zero'
			STOP
		END IF
		u(j) = (r(j)-a(j)*u(j-1))/bet
	END DO
	DO j=N-1,1,-1
		u(j) = u(j) - gam(j+1)*u(j+1)
	END DO
	RETURN
END SUBROUTINE tridiag	

subroutine polyfit(n,X,Y,m,coeffs)

integer n, m  ! number of points and degree of polynomial

integer i,ij,j,k,n1,m1,m2
real*8  C(m+2,m+2), A(m+2), B(m+2), Xc(m+2), Yx(m+2)
real*8  p,xx,s,yc
real*8  X(n), Y(n), coeffs(m+1)

coeffs = 0.d0
  n1=n; m1=m+1; m2=m+2
  do k=1, m2
    Xc(k)=0.d0
    do i=1, n1
	  Xc(k) = Xc(k) + X(i)**k
    end do
  end do
  yc=0.d0
  do i=1, n1  
    yc = yc + Y(i)
  end do
  do k=1, m
    Yx(k)=0.d0
    do i=1, n1
	  Yx(k) =  Yx(k) + Y(i)*X(i)**k
    end do 
  end do
  do i=1, m1
	do j=1, m1
      ij=i+j-2
      if (i==1.and.j==1)  then
	    C(1,1) = n1
      else 
	    C(i,j)=Xc(ij)
      end if 
    end do
  end do
  B(1)=yc; 
  do i=2,m1
    B(i)=Yx(i-1)
  end do 
  do k=1, m
    do i=k+1, m1
      B(i) = B(i) - C(i,k)/C(k,k)*B(k)
      do j=k+1, m1
        C(i,j) = C(i,j) - C(i,k)/C(k,k)*C(k,j)
      end do
    end do
  end do
  coeffs(m1)=B(m1)/C(m1,m1)
  do i=m, 1, -1
    s=0.d0
    do k=i+1, m1  
	  s = s + C(i,k)*coeffs(k)
    end do 
    coeffs(i) = (B(i)-s)/C(i,i)
  end do
end subroutine polyfit