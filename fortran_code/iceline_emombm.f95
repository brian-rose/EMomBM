!  compile using "gfortran -xf95-cpp-input"
!    to invoke the preprocessor

PROGRAM iceline_emombm

use emombm
IMPLICIT NONE

!  Declarations
	INTEGER, PARAMETER :: maxSizeS0array = 5000
	INTEGER i, ice, ess
	REAL*8, DIMENSION(J) :: Quold, Qlold, Tsold, Qunew, Qlnew, Tsnew
	REAL*8 ::  Taveold, Tavenew
	!REAL*8, PARAMETER :: S0max = 1450., S0min = 1150., deltaS0 = 1.
	!INTEGER, PARAMETER :: sizeS0array = NINT((S0max-S0min)/deltaS0*2) + 1
	REAL*8, DIMENSION(maxSizeS0array) :: S0array
	INTEGER :: sizeS0array, ios
	CHARACTER(LEN=100) :: S0file,outputfile,iceoutfile
	REAL*8 value, joe
	REAL*8 :: maxS0,minS0,deltaS0
	LOGICAL done
	
	namelist /icelineparam/ S0file,outputfile,iceoutfile,maxS0,minS0,deltaS0

	open(unit=8,file='paramfile')
	read(unit=8,nml=icelineparam)
	close(unit=8)
	
	call emombm_setup()

	!  read in the array of solar constants
	S0array = -1. !pad the array with extra -1 values... we declare a larger array than we need so that we can run with different S0 files without recompiling.
	!open(unit=9,file=S0file)
	!ios = 0
	!ess = 0
	!DO WHILE(ios==0)
	!	read(unit=9,fmt=103,iostat=ios) value
	!	if (ios==0) then
	!		ess = ess+1
	!		S0array(ess) = value
	!	end if
	!END DO
	!close(unit=9)
	!sizeS0array = ess
	done = .FALSE.
	S0array(1) = maxS0
	ess = 2
	DO WHILE(.NOT. done)
		joe = S0array(ess-1) - deltaS0
		IF ( joe .GE. minS0 ) THEN
			S0array(ess) = joe
		ELSE
			S0array(ess) = S0array(ess-1)+deltaS0
			done = .TRUE.
		END IF
		ess = ess+1
	END DO
	done = .FALSE.
	DO WHILE(.NOT. done)
		joe = S0array(ess-1)+deltaS0
		IF (joe .LE. maxS0 ) THEN
			S0array(ess) = joe
		ELSE
			sizeS0array = ess-1
			done = .true.
		END IF
		ess = ess+1
	END DO
	
	print *, ' '
	print *, 'Running an iceline experiment with sizeS0array = '
	print *, sizeS0array
	print *, ' '
	
	!  initial condition
	Quold = dble(2.0)*f0/sqrt(2.)*sin(phiq)
	Qlold = Quold
	Taveold = 11.
	Tsold = 40.-40.*2./pi*phiq	
	
	! open output files
	open(unit=9, file=S0file)
	open(unit=10, file=outputfile)
	open(unit=11, file=iceoutfile)
	
	DO ess=1,sizeS0array
		S0 = S0array(ess)
		S = S0/4*(1. + s2/2*(3*sin(phiq)**2 -1) )
		print *, S0
time:	DO i=1,numsteps
			call timestep_emombm(Quold,Qlold,Taveold,Tsold,Qunew,Qlnew,Tavenew,Tsnew,ice)
			Quold  = Qunew
			Qlold = Qlnew
			Taveold = Tavenew
			Tsold = Tsnew
		END DO time
		PRINT *, ice
		write(unit=9,fmt=100) S0
		write(unit=10, fmt=99) Qunew,Qlnew
		write(unit=10,fmt=100) Tavenew
		write(unit=10,fmt=100) Tsnew
		write(unit=11,fmt=102) ice
	END DO
	
	close(unit=9)
	close(unit=10)
	close(unit=11)
	
	print *, 'Done. Go home.'
	
	99 format(E20.10)
	101 format(F14.2)
	100 format(F14.8)
	102 format(I10)
	103 format(E20.7)

END PROGRAM iceline_emombm