module tightBinding
	use fgsl
	implicit none

	! d perpendicular parameter of the MoS2 hexagonal lattice
	real*8, parameter :: dperp=3.135_c_double/2.0_c_double  
	real*8, parameter :: alat=1.836_c_double ! hexagonal lattice constant
	real(fgsl_double), parameter :: bohr = 0.529177_fgsl_double ! effective bohr radius
	real(fgsl_double), parameter :: Ryd = 13.605693_fgsl_double ! effective Rydberg constant
	complex*16,parameter:: im=(0.0_c_double,1.0_c_double) ! imaginary unit
	real*8,parameter:: pi=3.1415926

	

contains

	subroutine Set_TBmodel(a1,a2,be)
		! sets the details of the tight binding model to make a lattice
		! OUTPUTS: hexagonal lattice vectors
		
		implicit none
		real*8,intent(out):: a1(2),a2(2),be(2)
		logical orient !!! true is be along y
		
		orient=.false. ! orientation of the lattice

		! define hexagonal lattice vectors
		if (orient) then
			a1=alat*sqrt(3.0_c_double)*[-1.0_c_double/2.0_c_double,sqrt(3.0_c_double)/2.0_c_double]
			a2=alat*sqrt(3.0_c_double)*[1.0_c_double/2.0_c_double,sqrt(3.0_c_double)/2.0_c_double]
			be=alat*[0.0_c_double,1.0_c_double]
		else
			a1=alat*sqrt(3.0_c_double)*[0.0_c_double,1.0_c_double]
			a2=alat*sqrt(3.0_c_double)*[sqrt(3.0_c_double)/2.0_c_double,1.0_c_double/2.0_c_double]
			be=alat*[1.0_c_double,0.0_c_double]	
		end if


	end subroutine Set_TBmodel
  
	subroutine rectQD_coord(m,n,coordA,coordB)
		! builds a rectangular flake of MoS2 for selecting atoms for the integral
		! INPUTS: size of the flake in 2D
		! OUTPUTS: coordinates of atoms A & B
		
		implicit none
		integer,intent(in)::n,m
		real*8,allocatable,intent(out):: coordA(:,:),coordB(:,:)

		integer ile,i,j,ind(1)
		real*8 uniA(2,2),uniB(2,2),initpos(2),addA(2,2),addB(2,2),R1(2),R2(2)
		logical orient !!! true is be along y
		real*8 a1(2),a2(2),be(2)


		ile=n*m*2
		allocate(coordA(ile,2),coordB(ile,2))
		initpos=[0.0_c_double,0.0_c_double]
		orient=.false.

		a1=0.0_c_double; a2=0.0_c_double; be=0.0_c_double;
		call Set_TBmodel(a1,a2,be)

		if (orient) then
			R1=a2-a1
			R2=a2+a1
			uniA(1,:)=initpos
			uniB(1,:)=initpos+be
			uniB(2,:)=initpos+be-a1
			uniA(2,:)=initpos+a2
		else
			R1=-a1
			R2=-a1+2*a2
			uniA(1,:)=initpos
			uniB(1,:)=initpos+be
			uniB(2,:)=initpos+be-a2
			uniA(2,:)=initpos-a1+a2
		end if

		coordA=0.0_c_double
		coordB=0.0_c_double
		coordA(1:2,:)=uniA
		coordB(1:2,:)=uniB

		do i=0,n-1
			do j=0,m-1
				if (i.gt.0 .or. j.gt.0) then
					addA=uniA+spread(R1*j+R2*i,1,2)
					addB=uniB+spread(R1*j+R2*i,1,2)
					ind=sub2ind([m, n],[j+1],[i+1])
					coordA(2*(ind(1)-1)+1:2*ind(1),:)=addA
					coordB(2*(ind(1)-1)+1:2*ind(1),:)=addB					
				end if
			end do
		end do



	end subroutine rectQD_coord	
	 
	function sub2ind(siz,row,col)result(ind)

		implicit none
		integer,intent(in)::siz(2)
		integer,intent(in),dimension(:):: row,col
		
		integer,allocatable,dimension(:):: ind
		
		allocate(ind(size(row)))
		ind=(col-1)*siz(1)+row

	end function sub2ind

	subroutine ind2sub(siz,inde,row,col)

		implicit none
		integer,intent(in)::siz(2)
		integer,intent(in),dimension(:):: inde
		
		integer,allocatable,dimension(:),intent(out):: row,col
		
		allocate(row(size(inde)),col(size(inde)))
		col=(inde-1)/siz(1)+1
		row=inde-(col-1)*siz(1)
		
	end subroutine ind2sub	

	subroutine meshgrid2D(x,y,xx,yy)
		implicit none
		real*8,intent(in):: x(:),y(:)
		real*8,allocatable,intent(out):: xx(:,:),yy(:,:)
		
		integer ilex,iley
		
		ilex=size(x)
		iley=size(y)
		
		allocate(xx(ilex,iley),yy(ilex,iley))
		
		xx=spread(x,2,iley)
		yy=spread(y,1,ilex)
		
	end subroutine meshgrid2D
	
end module tightBinding