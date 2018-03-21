!  compile using "gfortran -xf95-cpp-input"
!    to invoke the preprocessor

PROGRAM main

use emombm
IMPLICIT NONE
	INTEGER i, ice
	REAL*8, DIMENSION(J) :: Quold, Qlold, Tsold, Qunew, Qlnew, Tsnew
	REAL*8 :: Taveold, Tavenew
	REAL*8, DIMENSION(J+1) :: Uu,Ul,tau
	
	call emombm_setup()

	print *, ' '
	print *, 'Running the EMomBM once'
	print *, ' ' 

	if (.NOT. readInputFile ) then
		!  fixed initial condition
		Quold = dble(2.0)*f0/sqrt(2.)*sin(phiq)
		Qlold = Quold
		Taveold = 11.
		Tsold = 40.-40.*2./pi*phiq
	else
		!  get initial condition from a file
		print *, 'Loading initial conditions from file ', inputfilename
		open(unit=10,file=inputfilename)
		read(unit=10,fmt=*) Quold,Qlold
		read(unit=10,fmt=*) Taveold
		read(unit=10,fmt=*) Tsold
		close(unit=10)
	end if
	
	!  main time loop
time:	DO i=1,numsteps
			call timestep_emombm(Quold,Qlold,Taveold,Tsold,Qunew,Qlnew,Tavenew,Tsnew,ice)
			Quold = Qunew
			Qlold = Qlnew
			Taveold = Tavenew
			Tsold = Tsnew
			!PRINT *, ice
		END DO time
	
	! final winds
	CALL solveU(Quold,Qlold,Uu,Ul,J)
	tau = eps*(3*Ul - Uu)/2
	
	100 format(F14.8)
	101 format(E20.10)
	open(unit=10,file='output.dat')
	write(unit=10,fmt=101) Qunew,Qlnew
	write(unit=10,fmt=100) Tavenew
	write(unit=10,fmt=100) Tsnew
	write(unit=10,fmt=101) Uu,Ul,tau
	close(unit=10)
	print *, 'Results saved to file output.dat'
END PROGRAM main

