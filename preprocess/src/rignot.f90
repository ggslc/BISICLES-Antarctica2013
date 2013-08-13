
subroutine injectposdef(b,xb,yb,ewnb,nsnb,a,xa,ya,ewna,nsna)
  !inerpolate b(xb,yb) onto  a(xa,ya) onto 
  !b has finer resolution  a 
  !a,b are positive definite and 0.0 is a missing value

  integer :: ewna , nsna
  real(kind=8), dimension(1:ewna,1:nsna) :: a
  real(kind=8), dimension(1:ewna) :: xa
  real(kind=8), dimension(1:nsna) :: ya
  integer :: ewnb , nsnb
  real(kind=8), dimension(1:ewnb,1:nsnb) :: b
  real(kind=8), dimension(1:ewnb) :: xb
  real(kind=8), dimension(1:nsnb) :: yb
 
  real(kind=8) :: dxa,dxb,ww,tol
  integer ia,ja,ib,jb,ioff,joff,ki,kj
  real(kind=8) , dimension(0:1,0:1) :: w
  
  !assume any data less tha tol is missing
  tol = 1.0d-6

  dxa = xa(2)-xa(1)
  dxb = xb(2)-xb(1)
  
  ib = 2
  do ia = 1,ewna
     !advance ib until xb(ib) > x(ia)
     
     do while ( (ib.le.ewnb).and.(xb(ib)  .le. xa(ia)))
        ib = ib + 1
     end do
     
     jb = 2
     do ja = 1,nsna
        !advance jb
        do while ((jb.le.nsnb).and.(yb(jb) .le. ya(ja)))
           jb = jb + 1
        end do

        if ((ib.le.ewnb).and.(jb.le.nsnb)) then
       
        !now the cell center ib,jb is northeast of ia,ja
        w(0,0) = (xb(ib)-xa(ia))*(yb(jb)-ya(ja)) 
        if ( b(ib-1,jb-1)  .lt. tol)  w(0,0) = 0.0d0

        


        w(1,1) = -(xb(ib-1)-xa(ia))* (-(yb(jb-1)-ya(ja)))
        if ( b(ib,jb)  .lt. tol)  w(1,1) = 0.0d0
        
        w(0,1) = (xb(ib)-xa(ia))* (-(yb(jb-1)-ya(ja)))
        if ( b(ib-1,jb)  .lt. tol)  w(0,1) = 0.0d0

        w(1,0) = -(xb(ib-1)-xa(ia))*(yb(jb)-ya(ja))
        if ( b(ib,jb-1)  .lt. tol)  w(1,0) = 0.0d0

        if (minval(w) .lt. -tol) then
           write(*,*) "negative weight error, w = ",w, & 
                "dx = " ,xb(ib)-xa(ia), "dy = ", yb(jb)-ya(ja)
           stop
        end if

        ww = w(0,0) + w(1,0) + w(0,1) + w(1,1)
        
        if  ( (ww.le.(dxb*dxb + tol)) .and. (ww.gt.tol)) then 
           a(ia,ja) = ( w(0,0) * b(ib-1,jb-1) &
                + w(0,1) * b(ib-1,jb) + w(1,0) * b(ib,jb-1) &
                + w(1,1) * b(ib,jb) ) / ww
        else
           a(ia,ja) = 0.0d0
        end if

        else
           ! we are off the map
           a(ia,ja) = 0.0d0
           
        end if

     end do
  end do


end subroutine injectposdef


program rignot
  use ncio
  implicit none
  
  character(len=512) :: exename, nmlfile
  character(len=512) :: infile,outfile
  namelist /rignotveldata/ infile,outfile

 

  !RIGNOT (900m) data
  integer, parameter :: ewnr = 6223, nsnr = 6223
  real(kind=8), dimension(1:ewnr,1:nsnr) :: umodr,tmpr
  real(kind=8), dimension(1:ewnr) :: xr
  real(kind=8), dimension(1:nsnr) :: yr


  !BEDMAP 2 (1km sized) data
  integer, parameter :: ewn = 6144, nsn = 6144
  real(kind=8), dimension(1:ewn,1:nsn) :: umod,umodc
  real(kind=8), dimension(1:ewn) :: x
  real(kind=8), dimension(1:nsn) :: y
  
  !subset (1km sized) data
  integer, parameter :: ewns = 1024, nsns = 1024
  real(kind=8), dimension(1:ewns,1:nsns) :: umods,umodcs
  real(kind=8), dimension(1:ewns) :: xs
  real(kind=8), dimension(1:nsns) :: ys
  integer, dimension(2) :: lo

  integer i,j,k,a,b,n,m,ilo,ihi,jlo,jhi,argc


  call get_command_argument(0,exename)
  argc = command_argument_count()
  if (argc .ne. 1) then
     write(*,*) "usage ", trim(exename),  &
          " <namelist.nml> "
     stop
  end if

  call get_command_argument(1,nmlfile)
  open(8,file=nmlfile, status='OLD', recl=80, delim='APOSTROPHE')
  !i/o files
  read(8,nml=rignotveldata)  
  
  !x and y-coords the 6144 * 6144 1km BEDMAP2 grid
  x(1) = -3071500.0d0
  do i = 2,ewn
     x(i) = x(i-1) + 1000.0d0
  end do

  y(1) = -3071500.0d0
  do i = 2,nsn
     y(i) = y(i-1) + 1000.0d0
  end do

  umod = 0.0
  umodc = 0.0


  !x and y-coords for  6223*6223 900m rignot data
  xr(1) =  -2800000.0 
  do i = 2,ewnr
     xr(i) = xr(i-1) + 900.0d0
  end do

  yr(1) = -2799800.0d0  
  do i = 2,nsnr
     yr(i) = yr(i-1) + 900.0d0
  end do

  write(*,*) 'reading data'
  !load the x- and y- velocity components, and re-order in y
  call ncloadonenoxy(umodr,infile,"vx",ewnr,nsnr)
  call ncloadonenoxy(tmpr,infile,"vy",ewnr,nsnr)
  tmpr = sqrt(umodr*umodr + tmpr*tmpr)
  do j = 1,nsnr
     umodr(:,j) = tmpr(:,nsnr+1-j)
  end do

  write(*,*) 'interpolating'
  !interpolate
  umod = 0.0d0
  call injectposdef(umodr,xr,yr,ewnr,nsnr,umod,x,y,ewn,nsn)
  

  write(*,*) 'saving'


call ncsaveone(x,y,umod,ewn,nsn,outfile,"umod")


  !load the error estimate, and re-order in y
  !call ncloadonenoxybyte(tmpr,"Antarctica_ice_velocity.nc","err",ewnr,nsnr)
  !do j = 1,nsn
  !   umodr(:,j) = tmpr(:,nsnr+1-j)
  !end do
  ! umodc = 0.0d0
  !call injectposdef(umodr,xr,yr,ewnr,nsnr,umodc,x,y,ewn,nsn)
  !call ncsaveone(x,y,umodc,ewn,nsn,"Antarctica_ice_velocity-err-1km.nc","umodc")
  

end program rignot
