!  compile using "gfortran -xf95-cpp-input"
!    to invoke the preprocessor

PROGRAM main

IMPLICIT NONE
#include "SIZE.h" 

	INTEGER, PARAMETER :: maxSizeS0array = 5000
	REAL*8, DIMENSION(J) :: Qu, Ql, Ts	
	REAL*8 :: Tave, targetS0, value
	INTEGER :: ios, ess, numHits
	CHARACTER(LEN=30) :: S0file,inputfile,S0string
	CHARACTER(LEN=7) :: outputfileroot
	CHARACTER :: numString
	REAL*8, DIMENSION(maxSizeS0array) :: S0array
	INTEGER, DIMENSION(10) :: indices
	
	S0file = 'S0array.dat'
	inputfile = 'iceline_output.dat'
	outputfileroot = 'output_'
	print *, 'Enter the target value for S0'
	read *, targetS0
	print *, 'Searching file ' // inputfile // ' for target S0 = ', targetS0
	
	S0array = -1. !pad the array with extra -1 values... we declare a larger array than we need so that we can run with different S0 files without recompiling.
	open(unit=9,file=S0file)
	open(unit=10,file=inputfile)
	ios = 0
	ess = 0
	numHits = 0
	DO WHILE(ios==0)
		read(unit=9,fmt=*,iostat=ios) value
		if (ios==0) then
			ess = ess+1
			S0array(ess) = value
			read(unit=10,fmt=*) Qu,Ql
			read(unit=10,fmt=*) Tave
			read(unit=10,fmt=*) Ts
			if (value==targetS0) then
				numHits = numHits + 1
				write(numString,'(i1)') numHits
				open(unit=11,file=(outputfileroot // numString // '.dat'))
				write(unit=11,fmt=101) Qu,Ql
				write(unit=11,fmt=100) Tave
				write(unit=11,fmt=100) Ts
				close(unit=11)
			end if
		end if
	END DO
	close(unit=9)
	close(unit=10)
	
	100 format(F14.8)
	101 format(E20.10)
	103 format(E20.7)

	print *, 'Found ', numHits, ' matches.'

END PROGRAM main

