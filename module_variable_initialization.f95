! ############################################################################################### !
!              Canonical Monte Carlo algorithm for ellipsoid-of-revolution molecules              !
!      This module initialize common variables (number of particles, reduced number density,      !
!       reduced temperature etc.), Monte Carlo parameters (total number of cycles, number of      ! 
! equilibration cycles etc.), and Block Averaging parameters (maximum/minimum number of blocks).  !
!  This module also initialize some inquiry (character) variables, allowing the user to control   !
! which results will be written out in external files and to enable post-processing subroutines.  !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Original Developer: Joyce Tavares Lopes                             !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                            June 20th                                            !
! ############################################################################################### !
! Main Reference:                   M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !

module initvar

	! Uses one module: global variables
	use globalvar

	implicit none

	! *************************************************************************************** !
	!                                 VARIABLE INITIALIZATION                                 !
	!   Most variables should be first specified in an input file ('mchgo_input.ini' file).   !
	!          Some variables are on-the-fly (OTF) and must be specified by the user          !
	!                           at the beginning of the simulation.                           !
	!                   We provided an example .ini file to guide the user.                   !
	!                 Please also check our README file for more information.                 !
	! *************************************************************************************** !

	contains

	! *************************************************************************************** !
	!                            Initialization of common variables                           !
	! *************************************************************************************** !
	subroutine common_var()

		implicit none

		! *******************************************************************************
		! Reduced number density (OTF variable)
		! *******************************************************************************
		rdensity_loop: do
			print *, "Enter a reduced number density: "
			read ( *, * ) rho
			! Cannot be 0 or negative
			if ( rho > 0.d0 ) then
				exit rdensity_loop
			else
				print *, "Invalid answer. Enter a non-zero positive value."
			end if
		end do rdensity_loop

		! *******************************************************************************
		! Ellipsoid elongation parameter (OTF variable)
		! *******************************************************************************
		elongation_loop: do
			print *, "Enter an elongation for ellipsoidal particles: "
			read ( *, * ) elongation
			! Cannot be 0 or negative
			if ( elongation > 0.d0 ) then
				exit elongation_loop
			else
				print *, "Invalid answer. Enter a non-zero positive value."
			end if
		end do elongation_loop

		! *******************************************************************************
		! Output file descriptors (elongation [1] and reduced number density [2])
		! *******************************************************************************
		!  Might be necessary to change formats if the number of decimal places is
		!  higher than 2 or the number of characters is higher than 4.
		! *******************************************************************************
		format_file1 = "(F4.2)"
		write ( descriptor_file1, format_file1 ) elongation
		format_file2 = "(F4.2)"
		write ( descriptor_file2, format_file2 ) rho

		! *******************************************************************************
		! External variables
		! *******************************************************************************
		open ( UNIT = 10, FILE = "mchgo_input.ini" )

		! *******************************************************************************
		! Number of particles
		! *******************************************************************************
		read ( 10, * ) get, n_particles

		! *******************************************************************************
		! Cube root check
		! *******************************************************************************
		!  Checks whether the number of particles entered by the user is valid based on 
		!  each molecular structure.
		! *******************************************************************************
		particle_loop: do
			if ( config_selec(1) ) then
				check_root = dble( n_particles ) ** ( 1.d0 / 3.d0 )
				if ( dabs( check_root - dnint( check_root ) ) <= 1.d-10 ) then
					exit particle_loop
				else
					print *, "Invalid number of particles. Exiting..."
					call exit()
				end if
			else if ( config_selec(2) ) then
				check_root = ( 0.5d0 * dble( n_particles ) ) ** ( 1.d0 / 3.d0 )
				if ( dabs( check_root - dnint( check_root ) ) <= 1.d-10 ) then
					exit particle_loop
				else
					print *, "Invalid number of particles. Exiting..."
					call exit()
				end if
			else if ( config_selec(3) ) then
				check_root = ( 0.25d0 * dble( n_particles ) ) ** ( 1.d0 / 3.d0 )
				if ( dabs( check_root - dnint( check_root ) ) <= 1.d-10 ) then
					exit particle_loop
				else
					print *, "Invalid number of particles. Exiting..."
					call exit()
				end if
			end if
		end do particle_loop

		! *******************************************************************************
		! Reduced temperature
		! *******************************************************************************
		read ( 10, * ) get, temp

		! *******************************************************************************
		! Number of attractive range parameters
		! *******************************************************************************
		read ( 10, * ) get, n_lambda

		allocate ( lambda(n_lambda) )

		! *******************************************************************************
		! Attractive range parameter (λ)
		! *******************************************************************************
		!  Should the user choose multiple attractive range parameters, please be aware
		!  that all parameters must be entered in the same line at the .ini input file.
		!  See our example .ini file for more details.
		! *******************************************************************************
		read ( 10, * ) get, lambda

		! *******************************************************************************
		! Quarternion angle in degrees (for initial configuration only)
		! *******************************************************************************
		read ( 10, * ) get, quaternion_angle

		close ( 10 )

	end subroutine common_var

	! *************************************************************************************** !
	!                        Initialization of Monte Carlo parameters                         !
	! *************************************************************************************** !
	subroutine montecarlo_var()

		implicit none

		! *******************************************************************************
		! Simulation parameters
		! *******************************************************************************
		open ( UNIT = 10, FILE = "mchgo_input.ini" )

		! Common variables will not be read again
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber

		! *******************************************************************************
		! Total number of cycles
		! *******************************************************************************
		read ( 10, * ) get, max_cycles

		! *******************************************************************************
		! Number of equilibration cycles
		! *******************************************************************************
		read ( 10, * ) get, n_equil

		! *******************************************************************************
		! Result saving frequency
		! *******************************************************************************
		!  [1    = highest frequency, all results will be written out for each cycle ]
		!  [1000 = lower frequency, all results will be written out every 1000 cycles]
		! *******************************************************************************
		read ( 10, * ) get, n_save

		! *******************************************************************************
		! Maximum displacement adjustment frequency
		! *******************************************************************************
		!  [1    = highest frequency, maximum displacement will be adjusted for each cycle   ]
		!  [200  = moderate frequency, maximum displacement will be adjusted every 200 cycles]
		! *******************************************************************************
		read ( 10, * ) get, n_adjust

		! *******************************************************************************
		! Maximum translational displacement
		! *******************************************************************************
		read ( 10, * ) get, max_trans

		! *******************************************************************************
		! Maximum rotational displacement
		! *******************************************************************************
		read ( 10, * ) get, max_rot

		close ( 10 )

	end subroutine montecarlo_var

	! *************************************************************************************** !
	!                      Initialization of Block Averaging parameters                       !
	! *************************************************************************************** !
	subroutine blockav_var()

		implicit none

		! *******************************************************************************
		! Block Averaging parameters
		! *******************************************************************************
		open ( UNIT = 10, FILE = "mchgo_input.ini" )

		! Common variables and simulation parameters will not be read again
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber

		! *******************************************************************************
		! Minimum number of blocks (see Allen and Tildesley, 2nd Edition, page 282)
		! *******************************************************************************
		read ( 10, * ) get, min_blocks

		! *******************************************************************************
		! Maximum number of blocks (see Allen and Tildesley, 2nd Edition, page 282)
		! *******************************************************************************
		read ( 10, * ) get, max_blocks

		close ( 10 )

	end subroutine blockav_var

	! *************************************************************************************** !
	!                       Initialization of Inquiry/Control variables                       !
	! *************************************************************************************** !
	subroutine inquery_var()

		implicit none

		integer*8	:: i	! Counter

		! *******************************************************************************
		! Inquiry variables
		! *******************************************************************************
		open ( UNIT = 10, FILE = "mchgo_input.ini" )

		! Previous variables and parameters will not be read again
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber
		read ( 10, * ) get, anynumber

		! *******************************************************************************
		! Trajectory inquiry
		! *******************************************************************************
		!  Inquires whether trajectory of particles shall be written out.
		!  Answer Y (Yes) to write out trajectories or
		!         N (No) to ignore it
		! *******************************************************************************
		read ( 10, * ) get, traj_inq

		! Transforms characters into logical variables
		if ( traj_inq == "y" .or. traj_inq == "Y" ) then
			traj_check = .true.
		else if ( traj_inq == "n" .or. traj_inq == "N" ) then
			traj_check = .false.
		end if

		! *******************************************************************************
		! TPT coefficients inquiry
		! *******************************************************************************
		!  Inquires whether the post-processing block averaging subroutine will be called.
		!  This subroutine calculates the first- and second-order TPT coefficients on the fly.
		!  Answer Y (Yes) to calculate TPT coefficients or
		!         N (No) to ignore it
		! *******************************************************************************
		read ( 10, * ) get, coef_inq

		! Transforms characters into logical variables
		if ( coef_inq == "y" .or. coef_inq == "Y" ) then
			coef_check = .true.
		else if ( coef_inq == "n" .or. coef_inq == "N" ) then
			coef_check = .false.
		end if

		! *******************************************************************************
		! Potential inquiry
		! *******************************************************************************
		!  Inquires whether all potential results will be written out or only production-related ones.
		!  Answer Y (Yes) to write out both equilibration- and production-related results or
		!         N (No) to write out only production-related results
		! *******************************************************************************
		read ( 10, * ) get, pot_inq

		! Transforms characters into logical variables
		if ( pot_inq == "y" .or. pot_inq == "Y" ) then
			potential_check = .true.
		else if ( pot_inq == "n" .or. pot_inq == "N" ) then
			potential_check = .false.
		end if

		close ( 10 )

	end subroutine inquery_var

end module initvar
