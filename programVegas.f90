
program monte
! Computation of all integrals for desired atoms and orbitals in an MoS2 flake

	use monteCarlo
	implicit none
	  
	integer i,nn,kile,oile,j,ilecal,ile,ii,oi,jj,jile,ilex,n,m
	integer Ncalc,nV,x,y,z,w,a,b,c,d,l,mm,ileA,NQ,Nphi
	integer orbvec(3,2),orbs(4),atoms(4),ijkl(4),orbitals(4)
	integer,allocatable:: calc(:),calcAB(:,:)	
	type(parameters), target :: params
	type(dispParameters), target :: dispParams
	integer(fgsl_size_t):: callNum 
	integer(fgsl_size_t),allocatable :: callvec(:) 
	real(fgsl_double) :: xl(6), xu(6)
	real*8 wher,alpha,stru,R(3),limQ0,Rz,Rabs
	real*8,allocatable:: xvec(:),kel(:,:),alphavec(:),kel2(:,:)
	real*8,allocatable:: coordA(:,:),coordB(:,:),coords(:,:,:),pol(:)
	character nam*41,tex*4, pathfile*8,job*14
	logical CorK
	
	real :: start, finish
	
	! SETUP
	pathfile='NN_ints/'
	params%Potential='Co'
	params%alpha=2.2d0	
	xl = -20.0_fgsl_double/bohr ; xu = 20.0_fgsl_double/bohr
	orbvec=reshape([-2,9,2,-1,0,1],[3,2])
	
	! log file
	write(tex,'(A2,A2)')'V_',params%Potential
	write(job,'(A4,A10)') tex, '_vegas.dat'
	open(unit=11,file=pathfile//job,action='write',status='replace')
	write(11,*) 'init:',initNum
	
	
	! build flake of MoS2
	m=3
	n=3
	call rectQD_coord(m,n,coordA,coordB)
	allocate(coords(size(coordA,1),3,2))
	coords(:,1:2,1)=coordA; coords(:,3,1)=0.0_c_double
	coords(:,1:2,2)=coordB; coords(:,3,2)=1.0_c_double*dperp
	
	
	! list all atoms and orbitals to calculate integrals for
	Ncalc=1 ! how many orbitals of interest?
	nV=Ncalc**4
	write(11,*) 'Ncalc:',Ncalc
	write(11,*) 'nV:', nV
	allocate(calcAB(Ncalc,2))
	calcAB(:,1)=[2] ! orbital number list
	calcAB(:,2)=[1]	! atom number list
	write(11,*) calcAB(:,1)
	write(11,*) calcAB(:,2)
	flush(11)
	
	call cpu_time(start)
	
	! how many calculations of increased accuracy?
	ilecal=3
	allocate(callvec(ilecal))
	callvec=[500_fgsl_size_t,2000_fgsl_size_t,5000_fgsl_size_t]
	
	! LOOP to get all elements
	do j=1,nV
		x=(j-1)/Ncalc**3+1
		y=(j-(x-1)*Ncalc**3-1)/Ncalc**2+1
		z=(j-(x-1)*Ncalc**3-(y-1)*Ncalc**2-1)/Ncalc+1
		w=j-(x-1)*Ncalc**3-(y-1)*Ncalc**2-(z-1)*Ncalc
	
		! all possible atom combinations 
		ijkl=[calcAB(x,1),calcAB(y,1),calcAB(z,1),calcAB(w,1)] ! atom number
		atoms=[calcAB(x,2),calcAB(y,2),calcAB(z,2),calcAB(w,2)] ! sublattice number
		
		dispParams%at1=ijkl(1)	  
		dispParams%at4=ijkl(4)	  
		params%Ion1=coords(ijkl(1),:,atoms(1))
		params%Ion4=coords(ijkl(4),:,atoms(4))
		dispParams%at2=ijkl(2)
		dispParams%at3=ijkl(3)
		params%Ion2=coords(ijkl(2),:,atoms(2))
		params%Ion3=coords(ijkl(3),:,atoms(3))
		
		do i=1,3**4
			a=(i-1)/3**3+1
			b=(i-(a-1)*3**3-1)/3**2+1
			c=(i-(a-1)*3**3-(b-1)*3**2-1)/3+1
			d=i-(a-1)*3**3-(b-1)*3**2-(c-1)*3
			
			! all possible orbital combinations
			orbitals=[a,b,c,d]
			orbs=[(orbvec(orbitals(jj),atoms(jj)),jj=1,4)] ! store orbital and sublattice			
			
			dispParams%orbi=orbs(1)
			dispParams%orbl=orbs(4)
			params%orbi=orbs(1)
			params%orbl=orbs(4)
			dispParams%orbj=orbs(2)
			dispParams%orbk=orbs(3)
			params%orbj=orbs(2)
			params%orbk=orbs(3)
				
				do ii=1,ilecal
				
					callNum=fac*callvec(ii)
					print *, "callNum", callNum
					
					call Vegas(xl,xu,callNum,params,dispParams)
					write(11,*) '===============================',char(10)
					flush(11)
					
				end do
		
		end do
	end do
	
	call cpu_time(finish)
	write(11,*) 'TOTAL TIME in secs= ',finish-start
			
	close(11)
	
end program monte
