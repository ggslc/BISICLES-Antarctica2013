

subroutine coarsen(nf,xf,yf,nc,xc,yc,af,ac,ncomp)
  implicit none
  integer nf,nc,ncomp,i,j,ii,jj
  real(kind=8), dimension(1:nf) :: xf,yf
  real(kind=8), dimension(1:nc) :: xc,yc
  real(kind=8), dimension(1:nf,1:nf,1:ncomp) :: af
  real(kind=8), dimension(1:nc,1:nc,1:ncomp) :: ac

  if (nc .ne. nf/2) then
     write(*,*) "nc .ne. nf/2 makes no sense"
     stop
  end if

  do j = 1,nc
     jj = 2*j - 1
     xc(j) = 0.5d0*(xf(jj)+xf(jj+1))
     yc(j) = 0.5d0*(yf(jj)+yf(jj+1))
     do i = 1,nc
        ii = 2*i - 1
        ac(i,j,:) = 0.25d0 * & 
             (af(ii,jj,:) + af(ii+1,jj,:) +  af(ii,jj+1,:) + af(ii+1,jj+1,:))
     end do
  end do
 
end subroutine coarsen

subroutine savesubset(x,y,topg,thk,usrf,umod,umodc,mask,nx,ny,ixlo,ixhi,iylo,iyhi,file)
  use ncio
  implicit none
  integer, intent(in) :: ixlo,ixhi,iylo,iyhi,nx,ny
  
  real (kind=8), dimension(1:nx), intent(in) :: x
  real (kind=8), dimension(1:ny), intent(in) :: y
  real (kind=8), dimension(1:nx,1:ny), intent(in) :: topg,thk,usrf,umod,umodc,mask
  character(len=*), intent(in) :: file


  real (kind=8), dimension(ixlo:ixhi,iylo:iyhi) :: a
  real (kind=8), dimension(ixlo:ixhi) :: xs
  real (kind=8), dimension(ixlo:ixhi) :: ys

  xs = x(ixlo:ixhi)
  ys = y(iylo:iyhi) 
  a = topg(ixlo:ixhi,iylo:iyhi)
  call ncsaveone(xs,ys,a,ixhi-ixlo+1,iyhi-iylo+1,file,"topg")

  a = thk(ixlo:ixhi,iylo:iyhi) 
  call ncaddone(xs,ys,a,ixhi-ixlo+1,iyhi-iylo+1,file,"thk")

  a = usrf(ixlo:ixhi,iylo:iyhi) 
  call ncaddone(xs,ys,a,ixhi-ixlo+1,iyhi-iylo+1,file,"usrf")

  a = umod(ixlo:ixhi,iylo:iyhi) 
  call ncaddone(xs,ys,a,ixhi-ixlo+1,iyhi-iylo+1,file,"umod")

  a = umodc(ixlo:ixhi,iylo:iyhi) 
  call ncaddone(xs,ys,a,ixhi-ixlo+1,iyhi-iylo+1,file,"umodc")


  a = mask(ixlo:ixhi,iylo:iyhi) 
  call ncaddone(xs,ys,a,ixhi-ixlo+1,iyhi-iylo+1,file,"mask")

end subroutine savesubset

program subsetbm2
!load the bedmap2 binary topography,thickness,and mask data and produce a netcdf file
 
  use ncio
   implicit none

   character(len=512) :: exename, nmlfile
   character(len=512) :: inthkfile,inusrffile,intopgfile,inmaskfile,invelfile,outfile
   namelist /bedmap2data/ inthkfile,inusrffile,intopgfile,invelfile,inmaskfile,outfile


   real (kind=8), parameter :: rhoi = 918.0d0, rhoo = 1028.0d0, grav = 9.81d0, eps = 1.0e-3
  integer, parameter :: ewn=6667,nsn = 6667,nc=4, unit = 16, & 
       n = 6144, ewmin= 263 ,nsmin = 263 ,ewmax=ewmin+n-1,nsmax=nsmin+n-1
  real(kind=4), dimension(1:ewn,1:nsn) :: sp

  ! velocity data is on a smaller domain
  integer, parameter :: ewnvel = 5602, nsnvel = 5602,  &
       ewvelmin = 534, nsvelmin = 534, ewvelmax=ewvelmin+ewnvel-1, nsvelmax=nsvelmin+nsnvel-1  
  real(kind=4), dimension(1:ewnvel,1:nsnvel) :: spv 

  real(kind=4), dimension(1:1000,1:1000) :: t

  !data on 1 km grid
  real(kind=8), dimension(1:n,1:n) :: usrf,topg,thk,mask,umod,umodc
  real(kind=8), dimension(1:n) :: x
  real(kind=8), dimension(1:n) :: y

  character(len=32) :: file
  real(kind=8) :: dx
  integer i,j,ns,ilo,ihi,jlo,jhi, argc
  

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
  read(8,nml=bedmap2data)

  dx = 1.0e+3
  x(1) = -3333500.0d0 + real(ewmin-1)*dx
  y(1) = -3333500.0d0+  real(nsmin-1)*dx
  
  do i = 2,n
     x(i) = x(i-1) + dx
     y(i) = y(i-1) + dx
  end do
  
  open(unit,file=intopgfile, ACCESS='STREAM', FORM='UNFORMATTED')
  read(unit) sp
  close(unit)
  do i = 0,n-1
     topg(:,n-i) = sp(ewmin:ewmax,i+nsmin)
  end do

  open(unit,file=inthkfile, ACCESS='STREAM', FORM='UNFORMATTED')
  read(unit) sp
  close(unit)
  do i = 0,n-1
     thk(:,n-i) = sp(ewmin:ewmax,i+nsmin)
  end do
  
  where (thk(:,:).lt.eps)
     thk(:,:) = 0.0d0
  end where

  open(unit,file=inmaskfile, ACCESS='STREAM', FORM='UNFORMATTED')
  read(unit) sp
  close(unit)
  do i = 0,n-1
     mask(:,n-i) = sp(ewmin:ewmax,i+nsmin)
  end do
  !set mask to -1 for ice free regions
  where (thk.lt.eps)
      mask = -1
  end where


  open(unit,file=inusrffile, ACCESS='STREAM', FORM='UNFORMATTED')
  read(unit) sp
  close(unit)
  do i = 0,n-1
     usrf(:,n-i) = sp(ewmin:ewmax,i+nsmin)
  end do
  where (usrf.lt.eps)
     usrf = 0.0
  end where

  sp = 0.0e0
  umodc = 0.0
  open(unit,file=invelfile, ACCESS='STREAM', FORM='UNFORMATTED')
  read(unit) spv
  sp(ewvelmin:ewvelmax,nsvelmin:nsvelmax) = spv
  close(unit)
  do i = 0,n-1
     umod(:,n-i) = sp(ewmin:ewmax,i+nsmin)
  end do
  where (umod .gt. 0.0d0)
     umodc = 1.0
  elsewhere
     umodc = 0.0
     umod = 0.0
  end where
  

  call savesubset(x,y,topg,thk,usrf,umod,umodc,mask,n,n,1,n,1,n,outfile)



end program subsetbm2
