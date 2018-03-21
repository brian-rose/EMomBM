module emombm
implicit none
! preprocessor directive that defines grid size parameter J
#include "SIZE.h" 

!  declare parameter variables that will be read in from input file
	INTEGER :: numsteps
	LOGICAL :: readInputFile
	CHARACTER(LEN=50) :: inputfilename
	REAL*8 :: pi,delt,a,Tf,S0,s2,ai,a0,a2,Aup,Bup,Aout,Bout,Ca,Co	
	REAL*8 :: p0,f0,R,g,eps,Ld,smallkl,mu	
	REAL*8 :: dphi,Tavefactor!,Qmodfactor
	REAL*8, DIMENSION(J+1) :: phiu
	REAL*8, DIMENSION(J) :: phiq, S, coalbedo_o, coalbedo_i, coalbedo,taufactor,Qmodfactor
	REAL*8, DIMENSION(J+1) :: Ku, Kl, Ko, Ai1, Ai2, Ai3, Ad2,Kofactor
	REAL*8, DIMENSION(J) :: Al1,Al2,Al3,Au1,Au2,Au3,As1,As2,As3

contains

	subroutine emombm_setup()
		INTEGER i
		CHARACTER :: text
		INTEGER :: ios = 0
		
		namelist /commonparam/ numsteps,pi,delt,a,Tf,S0,s2,ai,a0,a2,Aup,Bup,Aout,Bout,Ca,Co,readInputFile,inputfilename
		namelist /emombmparam/ f0,R,g,p0,eps,Ld,smallkl,mu

		open(unit=8,file='paramfile')
			!  Echo the input files
		print *, ' '
		print *, 'Initializing the EMomBM with number of grid point:'
		print *, J
		print *, 'and with the following parameters:'
		print *, ' ' 
		print *, '*************  start of file PARAMFILE  **************'
		do while (ios/=-1)
			read(unit=8,fmt='(A)',iostat=ios,advance='NO') text
			if (ios==0) then
				write(*,'(A)',advance='NO') text
			elseif(ios==-2) then
				write(*,*) ' '
			end if
		end do
		print *, ' '
		print *, ' '
		print *, '*************  end of file PARAMFILE  **************'
		close(unit=8)
		
		open(unit=8,file='paramfile')
		read(unit=8,nml=commonparam)
		read(unit=8,nml=emombmparam)
		close(unit=8)
		
		dphi = pi/2/J
		!  here we use an implicit loop to initialize these arrays
		phiu = (/ ((i-1)*dphi, i=1,J+1) /)
		phiq = (/ ((i-0.5)*dphi, i=1,J) /)
		
		S = S0/4*(1. + s2/2*(3*sin(phiq)**2 -1) )  ! generates an array S
		coalbedo_o = a0 + a2/2*(3*sin(phiq)**2 -1)
		coalbedo_i = ai
		
		Tavefactor = delt/Ca/sum(cos(phiq))
		!Qmodfactor = delt*R/Ca/Ld**2/f0
		Qmodfactor = delt*R/Ca/Ld**2/f0/(1.-0.5*cos(phiq))
		taufactor = 2*delt*g/p0/a/cos(phiq)/dphi
		Kofactor = a*mu/f0/Co*cos(phiu)
		
		!  sets up the terms of the tridiagonal QGPV inversion operators
		Ai1 = (/ dble(0.), cos(phiu(1:J-1))/cos(phiq(1:J-1)), dble(0.) /)
		Ai2 = (/ dble(1.), -(1./cos(phiq(2:J))+1./cos(phiq(1:J-1)))*cos(phiu(2:J)), dble(1.) /)
		Ai3 = (/ dble(0.), cos(phiu(3:J+1))/cos(phiq(2:J)), dble(0.) /)		
		!Ad2 = Ai2 - 2./(Ld**2)*(a**2)*(dphi**2)
		Ad2 = Ai2 - 2./(Ld**2)/(1.-0.5*cos(phiu))*(a**2)*(dphi**2)

	end subroutine emombm_setup

	SUBROUTINE timestep_emombm(Quold,Qlold,Taveold,Tsold,Qunew,Qlnew,Tavenew,Tsnew,ice)
		REAL*8, DIMENSION(J) :: Tsold, Tsnew, Quold, Qlold, Qunew, Qlnew
		REAL*8 Taveold, Tavenew, maxUd, smallku
		INTEGER ice
		REAL*8, DIMENSION(J) :: Quhalf, Qlhalf, Tshalf, Fs, Qdot, Ta
		REAL*8, DIMENSION(J+1) :: Uu, Ul, Ud, tau, Y, Ko, Ku, Kl, curltau
		
		CALL solveU(Quold,Qlold,Uu,Ul,J)
		tau = eps*(3*Ul - Uu)/2
		!  Ko goes like curl(tau) squared
		curltau = (/ dble(0.), ( tau(3:J+1)*cos(phiu(3:J+1)) - tau(1:J-1)*cos(phiu(1:J-1)) ) /cos(phiu(2:J))/2/dphi, dble(0.) /)
		Ko = Kofactor*curltau**2
		
		Ud = Uu-Ul
		CALL temperature(Ud,Taveold,Ta)
		CALL heating(Ta,Tsold,S,coalbedo_o,coalbedo_i,ice,Fs,Qdot)
		Tavenew = Taveold + Tavefactor*sum(cos(phiq)*Qdot)
		Tshalf = Tsold + delt/Co*Fs
		Quhalf = Quold - Qmodfactor*Qdot
		Qlhalf = Qlold + Qmodfactor*Qdot + taufactor*( cos(phiu(2:J+1))*tau(2:J+1) - cos(phiu(1:J))*tau(1:J) )
		maxUd = maxval(abs(Ud))
		IF ( maxUd > 0 ) THEN
			Y = Ud / maxUd
		ELSE
			Y=0
		END IF
		CALL computeku(Quold,Qlold,Y,smallkl,smallku)
		IF ( smallku>0 ) THEN
			Ku = smallku*Y; Kl = smallkl*Y
		ELSE
			Ku = 0; Kl = 0;
		END IF
		Au1 = -cos(phiu(1:J))*Ku(1:J)/cos(phiq)*delt/a**2/dphi**2
		Au3 = -cos(phiu(2:J+1))*Ku(2:J+1)/cos(phiq)*delt/a**2/dphi**2
		Au2 = 1 + ( (/ dble(0.), -Au1(2:J) /) + (/ -Au3(1:J-1),dble(0.) /) )  ! atmosphere, main diagonal
		Al1 = -cos(phiu(1:J))*Kl(1:J)/cos(phiq)*delt/a**2/dphi**2
		Al3 = -cos(phiu(2:J+1))*Kl(2:J+1)/cos(phiq)*delt/a**2/dphi**2
		Al2 = 1 + ( (/ dble(0.), -Al1(2:J) /) + (/ -Al3(1:J-1),dble(0.) /) )
		As1 = -cos(phiu(1:J))*Ko(1:J)/cos(phiq)*delt/a**2/dphi**2
		As3 = -cos(phiu(2:J+1))*Ko(2:J+1)/cos(phiq)*delt/a**2/dphi**2
		As2 = 1 + ( (/ dble(0.), -As1(2:J) /) + (/ -As3(1:J-1),dble(0.) /) )  ! ocean, main diagonal
		SELECT CASE (ice)
		CASE (:0)  ! zero and below = no ice... no changes
		CASE (1)  ! ice at first gridpoint = snowball, no diffusion
			As2 = 1;  As1 =  0;  As3 = 0
		CASE(J:)  !  ice at last gridpoint only
			Tshalf(ice-1) = Tshalf(ice-1) - Tf*As3(ice-1)
			Tshalf(ice) = Tshalf(ice) + Tf*As1(ice)
			As2(ice) = 1;  As3(ice-1) = 0;
		CASE DEFAULT  ! ice edge in interior
			Tshalf(ice-1) = Tshalf(ice-1) - Tf*As3(ice-1)
			Tshalf(ice) = Tshalf(ice) + Tf*As1(ice)
			As2(ice:J) = 1; As3(ice-1:J) = 0;  As1(ice+1:J) = 0
		END SELECT
		CALL tridiag(Au1,Au2,Au3,Quhalf,Qunew,J)
		CALL tridiag(Al1,Al2,Al3,Qlhalf,Qlnew,J)  ! diffusion operator, atmosphere
		CALL tridiag(As1,As2,As3,Tshalf,Tsnew,J)  ! same for ocean	
	END SUBROUTINE timestep_emombm

	SUBROUTINE solveU(Qu,Ql,Uu,Ul,J)
		REAL*8, DIMENSION(J) :: Qu,Ql,Qi,Qd
		REAL*8, DIMENSION(J+1) :: Uu,Ul,Ui,Ud,Bi,Bd
		INTEGER J
		
		Qd = Qu-Ql
		Qi = Qu+Ql
		
		Bi = (/dble(0.), a*dphi*(4.*f0/sqrt(2.)*(sin(phiq(2:J)) - sin(phiq(1:J-1))) -(Qi(2:J)-Qi(1:J-1)) ), dble(0.)/)
		Bd = -a*dphi*(/ dble(0.), (Qd(2:J)-Qd(1:J-1)  ), dble(0.) /)
		
		CALL tridiag(Ai1,Ai2,Ai3,Bi,Ui,J)
		CALL tridiag(Ai1,Ad2,Ai3,Bd,Ud,J)
		Uu = (/ dble(0.), (Ui(2:J)+Ud(2:J))/2., dble(0.) /)
		Ul = (/ dble(0.), (Ui(2:J)-Ud(2:J))/2., dble(0.) /)
	END SUBROUTINE solveU
		
	SUBROUTINE temperature(Ud,Tave,Ta)
		REAL*8, DIMENSION(J+1) :: Ud
		REAL*8 :: Tave
		REAL*8, DIMENSION(J) :: Ta, temp
		INTEGER i
		
		temp(1) = (Ud(2)+Ud(1))/2
		DO i=2,J
			temp(i) = temp(i-1) + (Ud(i+1)+Ud(i))/2
		END DO
		Ta = -f0*a/R*dphi*(temp - sum(cos(phiq)*temp)/sum(cos(phiq))) + Tave		
	END SUBROUTINE temperature
	
	SUBROUTINE heating(Ta,Ts,S,coalbedo_o,coalbedo_i,ice,Fs,Qdot)
		INTEGER ice,i
		REAL*8, DIMENSION(J) :: Ta,Ts,S,coalbedo_o,coalbedo_i,coalbedo,Fs,Qdot,Fup
		!REAL*8 Aout,Bout,Aup,Bup,Tf
		LOGICAL found
		
		i = 0
		found = .FALSE.
		DO WHILE (.NOT. found)
			i = i + 1
			IF ( i>J ) THEN
				ice = -1 ! flag for no ice
				found = .TRUE.
			ELSEIF (Ts(i) < Tf) THEN
				ice = i
				found = .TRUE.
			END IF
		END DO
		coalbedo = coalbedo_o
		IF ( ice > 0 )  THEN
			coalbedo(ice:J) = coalbedo_i(ice:J)
		END IF
		Fup = Aup + Bup*(Ts-Ta)
		Fs = coalbedo*S - Fup
		Qdot = Fup - Aout - Bout*Ta
		RETURN
	END SUBROUTINE heating
	
	SUBROUTINE computeku(Qu,Ql,Y,smallkl,smallku)
		REAL*8, DIMENSION(J) :: Qu,Ql
		REAL*8, DIMENSION(J+1) :: Y
		REAL*8 :: smallkl, smallku, gradQubar, gradQlbar
		
		gradQubar = sum( cos(phiu)*Y*(/ dble(0.), Qu(2:J)-Qu(1:J-1), dble(0.) /) )
		gradQlbar = sum( cos(phiu)*Y*(/ dble(0.), Ql(2:J)-Ql(1:J-1), dble(0.) /) )
		IF ( gradQubar == 0 ) THEN
			smallku = 0
		ELSE
			smallku = -smallkl*gradQlbar/gradQubar
		END IF
	END SUBROUTINE computeku
end module emombm