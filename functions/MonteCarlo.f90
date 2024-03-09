module monteCarlo
	use fgsl
	use, intrinsic :: iso_c_binding
	use bessel
	use tightBinding
	implicit none

	! parameters of the integral
	type parameters
		real(fgsl_double):: Ion1(3),Ion2(3),Ion3(3),Ion4(3) ! atom positions
		integer :: orbi,orbj,orbk,orbl ! orbital angular momentum
		real*8 limQ0,alpha ! parameters for Keldysh screening
		character*2 potential ! choice of screening
	end type parameters

	!integral parameters to display
	type dispParameters
		integer:: at1,at2,at3,at4 ! atom numbers from the grid
		integer :: orbi,orbj,orbk,orbl ! orbital angular momentum
	end type dispParameters

	! number of evaluations for integration
	integer(fgsl_size_t), parameter :: fac = 1000 !keep constant
	integer(fgsl_size_t), parameter :: init = 100
	integer(fgsl_size_t), parameter :: initNum = init*fac
	! parameters of the Slater orbitals used throughout for MoS2
	real(c_double), parameter :: ksiM=3.111_c_double
	real(c_double), parameter :: ksiS=1.8273_c_double
	

  
contains

	function g(rvec, n, par) bind(c)
		! Integrand
		! INPUTS: - position vector in 3D
		!		  - size of the problem: 6 dims
		!		  - parameter  pointer
		
		integer(c_size_t), value :: n
		type(c_ptr), value :: rvec, par
		real(c_double) :: g

		real(c_double) :: rR1,rR2,pot,dr
		real(c_double), dimension(:), pointer :: v
		real*8 alpha,limQ0
		complex*16 orbi,orbj,orbk,orbl
		type(parameters), pointer :: params
		
		! assign targets from C to pointers in Fortran
		call c_f_pointer(rvec, v, [n])
		call c_f_pointer(par, params)
		
		! select orbital in the integral and build Slater orbitals
		orbi=selectOrbital(v,params%Ion1,params%orbi,ksiS,ksiM)
		orbj=selectOrbital(v,params%Ion2,params%orbj,ksiS,ksiM)
		orbk=selectOrbital(v,params%Ion3,params%orbk,ksiS,ksiM)
		orbl=selectOrbital(v,params%Ion4,params%orbl,ksiS,ksiM)

		! parameter for Keldysh interaction
		alpha=params%alpha

		! select Coulomb or Keldysh (default is Coulomb)
		select case(params%potential)
			case('Co') ! Coulomb
				! Coulomb interaction: 2/r (in effective Rydberg units)
				dr=sqrt((v(4)-v(1))**2+(v(5)-v(2))**2+(v(6)-v(3))**2)
				pot=1.0_c_double/dr
			case('Ke') ! Keldysh
				dr=sqrt((v(4)-v(1))**2+(v(5)-v(2))**2)				
				pot=keldyshFunR(dr,alpha)
				
			case default
				dr=sqrt((v(4)-v(1))**2+(v(5)-v(2))**2+(v(6)-v(3))**2)
				pot=1.0_c_double/dr		
		end select
		
		! integrand is psi1*psi2*potential*psi3*psi4
		g=2.0_c_double*conjg(orbi)*conjg(orbj)*pot*orbk*orbl*Ryd

	end function g
	
	
	
	function selectOrbital(v,Icoord,m,ksiS,ksiM) result(orbi)
		! selects one of 4 orbitals in the integral
		! INPUTS: - 3D position vector
		!		  - 3D position of one of the 4 atoms in the integral
		!		  - m_p/m_d angular momentum of the orbital, 
		!			9 stands for m_d=0, 0 is m_p=0
		!		  - ksi Slater orbital parameters for both atom types
		! OUTPUT: orbital function value at point v
		
		implicit none
		complex*16 orbi	
		real(c_double):: Ix,Iy,Iz
	
		Ix=Icoord(1)
		Iy=Icoord(2)
		Iz=Icoord(3)
	
		select case (m)
			case(-1) ! m_p=-1
				orbi=slaterOrb(v(1),v(2),v(3),Ix,Iy,Iz,3,1,-1,ksiS)/sqrt(2.0_c_double)
				orbi=orbi+slaterOrb(v(1),v(2),v(3),Ix4,Iy4,-Iz4,3,1,-1,ksiS)/sqrt(2.0_c_double)
			case(0) ! m_p=0
				orbi=slaterOrb(v(1),v(2),v(3),Ix,Iy,Iz,3,1,0,ksiS)/sqrt(2.0_c_double)
				orbi=orbi-slaterOrb(v(1),v(2),v(3),Ix4,Iy4,-Iz4,3,1,0,ksiS)/sqrt(2.0_c_double)
			case(1) ! m_p=1
				orbi=(slaterOrb(v(1),v(2),v(3),Ix,Iy,Iz,3,1,1,ksiS)/sqrt(2.0_c_double)
				orbi=orbi+slaterOrb(v(1),v(2),v(3),Ix4,Iy4,-Iz4,3,1,1,ksiS)/sqrt(2.0_c_double)
			case(-2) ! m_d=-2
				orbi=slaterOrb(v(1),v(2),v(3),Ix,Iy,Iz,4,2,-2,ksiM)
			case(9) ! m_d=0
				orbi=slaterOrb(v(1),v(2),v(3),Ix,Iy,Iz,4,2,0,ksiM)
			case(2) ! m_d=2
				orbi=slaterOrb(v(1),v(2),v(3),Ix,Iy,Iz,4,2,2,ksiM)		
		end select
	
	end function
    
	function keldyshFunR(r,alpha) result(valu)
		! calculates Keldysh screening to Coulomb interaction
		! INPUTS: 3D position vector, alpha parameter
		! OUTPUT: Keldysh screened interaction
		
		implicit none
		real*8,intent(in):: r,alpha
		real*8 valu
		
		real*8 struv0,argu,bess
		
		argu=r/2.0_c_double/pi/alpha
		call STVH0(argu,struv0)
		bess=BESSY(0,argu)
		valu=1.0_c_double/4.0_c_double/alpha*(struv0-bess)	
		
	end function keldyshFunR

	function slaterOrb(x,y,z,xc,yc,zc,n,l,m,ksi) result(valu)
		! builds a Slater orbital using spherical harmonics
		! INPUTS: - 3D position
		!		  - 3D position of the center
		!		  - spherical harmonics quantum numbers n,l,m
		!		  - ksi parameter
		
		implicit none
		real*8,intent(in):: x,y,z,xc,yc,zc,ksi
		integer,intent(in):: n,l,m
		complex*16 valu
		
		real*8 rad,coef
		integer i,nfac
		complex*16 spherical
		
		nfac = PRODUCT([(i,i=1,(2*n))])
		coef=(2.0_c_double*ksi)**n*sqrt(2.0_c_double*ksi/nfac)
		rad=sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)
		spherical=spherH(x-xc,y-yc,z-zc,l,m)
		valu=coef*rad**(n-1)*exp(-ksi*rad*1.0_c_double)*spherical*1.0_c_double

	end function slaterOrb
		 
	function spherH(x,y,z,l,m) result(valu)
		! spherical harmonics
		! INPUTS: - 3D position
		!		  - quantum numbers l,m
		! OUTPUT: spherical harmonic
		
		implicit none
		real*8,intent(in):: x,y,z
		integer,intent(in):: l,m
		complex*16 valu
		
		integer mab,i,lmmfac,lmpfac
		real*8 small,r,phi,zr,fac
		real*8,allocatable:: vals(:,:)
		
		allocate(vals(1,l+1))
		valu=(0.0_c_double,0.0_c_double); r=0.0_c_double
		small=0.00001_c_double
		
		mab=abs(m)	
		r=sqrt(x**2+y**2+z**2)
		phi=atan2(y,x)
		if (abs(r).lt.small) then
			zr=0.0_c_double
		else
			zr=z/r
		end if
		
		! call legendre polynomial
		call pmns_polynomial_value ( 1, l, mab, [zr], vals )		
		
		valu=vals(1,l+1)*exp(im*m*phi)


	end function spherH



	subroutine display_results(title,calls,resul,error,pars)
		! displays results of the Vegas algorithm convergence
		! INPUTS: - step title
		!		  - number of evaluations
		!		  - step result
		!		  - step convergence error
		!		  - parameters: atom numbers and orbitals etc
		
		implicit none
		character(kind=fgsl_char,len=*), intent(in) :: title
		real(fgsl_double), intent(in) :: result, error
		type(dispParameters), target :: pars
		integer(fgsl_size_t):: calls 

		write(11, '(A)') trim(title)
		write(11,*) 'calls:', calls
		write(11, *) 'ijkl:',pars%at1,pars%at2,pars%at3,pars%at4,'   orbitals:    ',pars%orbi,pars%orbj,pars%orbk,pars%orbl
		write(11, '(''result = '',F10.6)') resul
		write(11, '(''error  = '',F10.6)') error
		flush(11)

	end subroutine display_results

	subroutine Vegas(xl,xu,calls,params,dispParam)
		! computes the integral using Vegas algorithm
		! INPUTS: - boundaries of the integration
		!		  - number of evaluations
		!		  - parameters
		
		
		implicit none
		type(parameters), target :: params
		type(dispParameters), target :: dispParam
		type(fgsl_file) :: file
		type(fgsl_rng) :: r
		type(fgsl_rng_type) :: t
		type(fgsl_monte_vegas_state) :: v
		type(fgsl_monte_function) :: gfun
		type(c_ptr) :: ptr
		
		integer(fgsl_size_t):: calls 
		integer(fgsl_size_t) :: its
		integer(fgsl_int) :: status, stage, mode, verbose
		real(fgsl_double) :: chisq, res, err, y, yy, yyy
		real(fgsl_double) :: xl(6), xu(6)
		
		integer i
		
		ptr = c_loc(params)
		t = fgsl_rng_env_setup()
		r = fgsl_rng_alloc(t)

		!!! PREP
		! 6_fgsl_size_t denotes 6-dimensional integral
		! g is the integrated function
		! xl, xu are boundaries of integration
		! initNum parameter sets trial evaluation number, 
		! later bigger number "calls" is used
		gfun = fgsl_monte_function_init(g, 6_fgsl_size_t, ptr)
		v = fgsl_monte_vegas_alloc(6_fgsl_size_t)
		status = fgsl_monte_vegas_integrate(gfun, xl, xu, 6_fgsl_size_t, &
		   initNum, r, v, res, err)

		!!! MONTE CARLO
		do
			status = fgsl_monte_vegas_integrate(gfun, xl, xu, 6_fgsl_size_t, &
				  calls, r, v, res, err)
			call fgsl_monte_vegas_getparams(v, y, yy, chisq, yyy, &
			   its, stage, mode, verbose, file)			
			
			flush(11)
			if (abs(chisq - 1.0_fgsl_double) <= 0.5_fgsl_double) exit
		end do


		write(11,*) char(10)
		write(11,*) '==============================='
		call display_results('VEGAS CONVERGED',calls,res,err,dispParam)	  
		call fgsl_monte_vegas_free(v)
		call fgsl_monte_function_free(gfun)
		call fgsl_rng_free(r)

	end subroutine


  
end module monteCarlo

