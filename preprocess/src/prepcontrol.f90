!prepcontrol.f90
!Read geoemtry data , speed data, sector data 
!Estimate initial friction corigin
!corigin = rho * g * h * |grad(s')|/|u'|
! s' and u' are modified from |s| and |u| to 
!   (i) avoid noisy grad(s)
!   (ii) fill missing data regions 
! crop data to a nx * ny rectanglar region , bottom corner is ixlo,iylo
! relatice to botto left corner of input grids
!Write arrays required by control problem to outctrlfile
!corigin, speed umod , speed error coefficient umodc (1/2*sigma^2), 
!plus optional array qmask >=0 . corigin is set to csticky  and umod set to 0 where qmask = 0


module physparam
  implicit none
  real(kind=8), parameter :: rhoi=918.0d0 , rhoo=1028.0d0, gravity = 9.81d0
  real(kind=8), parameter :: csticky = 1.0d+5 , cslippy = 1.0d+1, cmin = 1.0d+1
  real(kind=8), parameter :: sscale = 8.0d+3 ! scale length used to smooth s before computing gradient
  real(kind=8), parameter :: cm = 1.0d-2, mm = 1.0d-3, tiny = 1.0d-6 ! small numbers
end module physparam

module globaldata

  character(len=512) :: ingeofile, invelfile, insecfile, intempfile, &
       outgeofile, outctrlfile, outsecfile, outtempfile
  
  character(len=16) :: buf
  integer :: argc,nsigma,ixlo,iylo,nx,ny,nxin,nyin,nx4,ny4,nxin4,nyin4,ixlo4,iylo4
  integer, dimension(32) :: qsectors

  namelist /files/  ingeofile, invelfile, insecfile, intempfile, & 
       outgeofile, outctrlfile, outsecfile, outtempfile
  namelist /dims/  ixlo,iylo,nx,ny,nxin,nyin,nsigma
  namelist /quiescent/ qsectors 


end module globaldata

subroutine growmin(a,nx,ny,niter)
  !grow a, copying the minimum value
  integer,  intent(in) :: nx,ny,niter
  real(kind=8), dimension(1:nx,1:ny), intent(inout) :: a
  real(kind=8), dimension(1:nx,1:ny) :: b
  integer i,j,iter
 do iter = 1,niter
    b = a
     do j = 2,ny-1
        do i = 2,nx-1
           b(i,j) = min(a(i-1,j),b(i,j))
           b(i,j) = min(a(i+1,j),b(i,j))
           b(i,j) = min(a(i,j-1),b(i,j))
           b(i,j) = min(a(i,j+1),b(i,j))
        end do
     end do
     a = b
  end do

end subroutine growmin

subroutine growposdef(ug,u,nx,ny,tol,niter)
  !grow a positive definite field u, such that
  !u > tol everywhere 
  integer,  intent(in) :: nx,ny,niter
  real(kind=8), intent(in) :: tol
  real(kind=8), dimension(1:nx,1:ny), intent(inout) :: ug
  real(kind=8), dimension(1:nx,1:ny), intent(in) :: u
  integer i,j,iter,wp,we,ww,ws,wn,wsum,n

  n = 0

  ug = u
  do iter = 1,niter
     do j = 2,ny-1
        do i = 2,nx-1
           if (u(i,j) .lt. tol) then
              if (iter.le.0) then
                 ug(i,j) = max(ug(i,j), ug(i-1,j))
                 ug(i,j) = max(ug(i,j), ug(i+1,j))
                 ug(i,j) = max(ug(i,j), ug(i,j-1))
                 ug(i,j) = max(ug(i,j), ug(i,j+1))
              else
                 if (ug(i-1,j) .gt. tol) then
                    ww = 1.0d0 
                 else
                    ww = 0.0d0
                 end if
                 
                 if (ug(i+1,j) .gt. tol) then
                    we = 1.0d0 
                 else
                    we = 0.0d0
                 end if
                 
                 if (ug(i,j-1) .gt. tol) then
                    ws = 1.0d0 
                 else
                    ws = 0.0d0
                 end if
                 
                 if (ug(i,j+1) .gt. tol) then
                    wn = 1.0d0  
                 else
                    wn = 0.0d0
                 end if
                 
                 wsum = wn + we + ws + ww
                 if (wsum.gt.0.0d0) then
                    ug(i,j) = (ww * ug(i-1,j) & 
                         + we * ug(i+1,j) &
                         + ws * ug(i,j-1) &
                         + wn * ug(i,j+1))/wsum
                    
                 end if
              end if
           end if
        end do
     end do
  end do
  
  where (ug.lt.tol)
     ug = tol
  end where


end subroutine growposdef

subroutine smooth(v,a,dx,ewn,nsn,niter)
  ! smooth a field by solving v(out) - a^2 v(out)" = v(in)
  use mgrelax
  implicit none
  integer :: ewn,nsn
  integer :: iter,niter
  real (kind = 8), dimension(1:ewn,1:nsn) :: v
  real (kind = 8), dimension(0:ewn+1,0:nsn+1) :: umg,vmg,dumg,rmg
  real (kind=8) :: mu, dx, a, resNorm
  
  mu = (a**2) / (dx**2)
  
  rmg = 0.0d0
  vmg = 0.0d0
  vmg(1:ewn,1:nsn) = v
  umg = vmg
  call resid(umg,vmg,rmg,mu,ewn,nsn)
  resNorm = sum(abs(rmg) )
  write(*,*) 'initial res= ', resNorm
 
  do iter = 1, niter
     dumg = 0.0
     call vcycle(dumg,rmg,mu,ewn,nsn,10,4)
     umg = umg + dumg
     call resid(umg,vmg,rmg,mu,ewn,nsn)
     resNorm = sum(abs(rmg) )
     write (*,*) 'MG iter', iter, " norm = ", resNorm
     
  end do
  v = umg(1:ewn,1:nsn) 
  return
end subroutine smooth

subroutine computebtrc(topg, thk, umod, btrc, dx)
  use physparam
  use globaldata
  real(kind=8), intent(in) :: dx
  real(kind=8), dimension(1:nx,1:ny), intent(in) :: thk, topg, umod
  real(kind=8), dimension(1:nx,1:ny), intent(inout) :: btrc
  real(kind=8), dimension(1:nx,1:ny) :: usrf

  integer :: i,j
  real(kind=8) :: dsx,dsy
 
  usrf = max(topg + thk, (1.0d0-rhoi/rhoo)*thk)
  !smooth usrf before numerical differentiation
  call smooth(usrf,sscale,dx,nx,ny,20)
  
  do j = 2,ny -1
     do i = 2,nx-1   
        dsx = 0.5*(usrf(i+1,j) - usrf(i-1,j))/dx
        dsy = 0.5*(usrf(i,j+1) - usrf(i,j-1))/dx
        btrc(i,j) = dsqrt(dsx**2 + dsy**2 + 1.0d-6) * (thk(i,j)+1.0d0) * rhoi * gravity / (umod(i,j) + 1.0d-6)
     end do
  end do

  where (btrc.lt.cslippy)
     btrc = cslippy
  end where

  return
end subroutine computebtrc

module bmmods
  !a bunch of modifications to the bedmap2 geoemtry
   use physparam
  implicit none
 
contains
  subroutine computebox(ncells,ilo,ihi,jlo,jhi,xmin,ymin,xmax,ymax,x,y,nx,ny)
    !find the largest box such that xmin < x < xmax and ymin < y < ymax
    !which fits inside the box x(1) < x < x(nx) and y(1) < y < y(ny)

    integer, intent(out) :: ncells,ilo,ihi,jlo,jhi
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: xmin,xmax,ymin,ymax
    real(kind=8), dimension(1:nx), intent(in) :: x
    real(kind=8), dimension(1:ny), intent(in) :: y

    real(kind=8) :: dx,dy

    dx = x(2) - x(1)
    dy = y(2) - y(1)

    ilo = 0
    jlo = 0
    ihi = 0
    jhi = 0
    ncells = 0
  
    if ( (x(nx).gt.x(1)) .and.(y(ny).gt.y(1)))  then
       if ( (xmax.gt.xmin) .and.(ymax.gt.ymin)) then

          

          if (xmax.ge.x(1))  ilo = max(1, ceiling ( ( xmin - x(1)) / dx))
          if (ymax.ge.y(1))  jlo = max(1, ceiling ( ( ymin - y(1)) / dx))
          if (xmin.le.x(nx))  ihi = min(nx, floor ( ( xmax - x(1)) / dx))
          if (ymin.le.y(ny))  jhi = min(ny, floor ( ( ymax - y(1)) / dx))
          ncells = (ihi-ilo+1)*(jhi-jlo+1)
       end if
    end if

    

    return
  end subroutine computebox


  subroutine removeremoteice(thk,topg,umod,x,y,nx,ny)
    !get rid of the various remote islands of ice around the peninsular
    !and a few other spots.
    
    integer, intent(in) :: nx,ny
    real(kind=8), dimension(1:nx,1:ny), intent(inout) :: thk,topg,umod
    real(kind=8), dimension(1:nx), intent(in) :: x
    real(kind=8), dimension(1:ny), intent(in) :: y
    integer(kind=1), dimension(1:nx,1:ny) :: c

    integer :: i,j,k,ilo,jlo,ihi,jhi,iter,ncells
    real(kind=8) :: dx,r
    dx = x(2) - x(1)
    r = (1-rhoi/rhoo)
    c = 0

    !set c = 1 in some grounded regions.
    !There needs to be at least one point covered, so fail if that 
    !is not the case
    call computebox(ncells,ilo,ihi,jlo,jhi,0.7d+6,-1.5d+6,1.6d+6,1.3d+6,x,y,nx,ny)
    
    if (ncells.gt.0) then
       where  (r * thk(ilo:ihi,jlo:jhi).lt.(thk(ilo:ihi,jlo:jhi) + topg(ilo:ihi,jlo:jhi)))
          c(ilo:ihi,jlo:jhi) = 1
       end where
    end if

    call computebox(ncells,ilo,ihi,jlo,jhi,-1.5d+6,-1.0d+6,2.0d+6,0.5d+6,x,y,nx,ny)
    if (ncells.gt.0) then
       where ((thk(ilo:ihi,jlo:jhi).gt.cm).and. & 
            (r * thk(ilo:ihi,jlo:jhi).lt.(thk(ilo:ihi,jlo:jhi) + topg(ilo:ihi,jlo:jhi))))
          c(ilo:ihi,jlo:jhi) = 1
       end where
    end if

    if (maxval(c).ne.1) then
       write(*,*) "no ice marked as continental, x(1),x(nx),y(1),y(ny) = ", x(1), x(nx), y(1), y(ny)
       stop 1
    end if

    !set c = 1.0 in every cell (i,j) where thk(i,j) > cm
    ! and thk(i+-1,j+-1) = cm

    do iter = 1,10
       k = 0
       do i = 2,nx-1
          do j = 2,ny-1
             if ((c(i,j).ne.1).and.(thk(i,j).gt.cm )) then
                if (c(i,j-1)+c(i,j+1)+c(i+1,j)+c(i-1,j).ge.1) then
                   c(i,j) =   1 
                   k = k + 1
                end if
             end if
             
          end do
          
       end do
       
       do j = 2,ny-1
          do i = 2,nx-1
          if ((c(i,j).ne.1).and.(thk(i,j).gt.cm )) then
             if (c(i,j-1)+c(i,j+1)+c(i+1,j)+c(i-1,j).ge.1) then
                c(i,j) =  1 
                k = k + 1
             end if
          end if
       end do
    end do
    
    do i = 2,nx-1
        do j = ny-1,2,-1
           if ((c(i,j).ne.1).and.(thk(i,j).gt.cm )) then
              if (c(i,j-1)+c(i,j+1)+c(i+1,j)+c(i-1,j).ge.1) then
                 c(i,j) =  1 
                 k = k + 1
              end if
           end if
        end do
     end do  

     
     do j = 2,ny-1
        do i = nx-1,2,-1
          if ((c(i,j).ne.1).and.(thk(i,j).gt.cm )) then
              if (c(i,j-1)+c(i,j+1)+c(i+1,j)+c(i-1,j).ge.1) then
                 c(i,j) =  1 
                 k = k + 1
              end if
           end if
        end do
     end do
     
     write(*,*) "iteration = ",iter," new cells connected = ",k
     if (k .eq. 0) exit
  
  end do


  !set thk = 0 where contact = 0
  where (c.eq.0)
     thk = 0.0d0
  end where

  end subroutine removeremoteice

  subroutine removethintongues(thk,topg,umod,x,y,nx,ny)
    !there are three thin floating tonugues which cause the solvers
    !to converge slowly. Since they don't affect the upstream 
    !dynamics, snip them off

    
    integer, intent(in) :: nx,ny
    real(kind=8), dimension(1:nx,1:ny), intent(inout) :: thk,topg,umod
    real(kind=8), dimension(1:nx), intent(in) :: x
    real(kind=8), dimension(1:ny), intent(in) :: y

    integer :: i,j,ilo,jlo,ihi,jhi,ncells
    real(kind=8) :: dx,r
    dx = x(2) - x(1)
    r = (1-rhoi/rhoo)

    !Tongue 1, East Antarctica near to Lutzow-Holm bay
    call computebox(ncells,ilo,ihi,jlo,jhi,1.35d+6,1.775d+6,1.4d+6,2.075d+6,x,y,nx,ny)
    if (ncells.gt.0) then
       where (r * thk(ilo:ihi,jlo:jhi).gt.(thk(ilo:ihi,jlo:jhi) + topg(ilo:ihi,jlo:jhi))) 
          thk(ilo:ihi,jlo:jhi) = 0.0d0
       end where
    end if

    !Tongue 2, East Antarctica Near to East Adelie land (Mertz tongue)   
    call computebox(ncells,ilo,ihi,jlo,jhi,1.39d+6,-2.14d+6,1.46d+6,-2.05d+6,x,y,nx,ny)
    if (ncells.gt.0) then
       where (r * thk(ilo:ihi,jlo:jhi).gt.(thk(ilo:ihi,jlo:jhi) + topg(ilo:ihi,jlo:jhi))) 
          thk(ilo:ihi,jlo:jhi) = 0.0d0
       end where
    end if

    !Tongue 3, East Antarctic near to McMurdo (Drygalski tongue)
    call computebox(ncells,ilo,ihi,jlo,jhi,0.36d+6,-1.55d+6,0.46d+6,-1.25d+6,x,y,nx,ny)
    if (ncells.gt.0) then
       where (r * thk(ilo:ihi,jlo:jhi).gt.(thk(ilo:ihi,jlo:jhi) + topg(ilo:ihi,jlo:jhi))) 
          thk(ilo:ihi,jlo:jhi) = 0.0d0
       end where
    end if

  end subroutine removethintongues

  subroutine removelakevostok(thk,topg,usrf,x,y,nx,ny)
    integer, intent(in) :: nx,ny
    real(kind=8), dimension(1:nx,1:ny), intent(inout) :: thk,topg,usrf
    real(kind=8), dimension(1:nx), intent(in) :: x
    real(kind=8), dimension(1:ny), intent(in) :: y

    integer :: i,j,ilo,jlo,ihi,jhi,ncells
    real(kind=8) :: dx,r

    dx = x(2) - x(1)
    r = rhoo/rhoi
    call computebox(ncells,ilo,ihi,jlo,jhi, & 
         1.2d+6,-0.4d+6,1.46d+6,-0.3d+6,x,y,nx,ny)
    if (ncells.gt.0) then
       topg(ilo:ihi,jlo:jhi) = usrf(ilo:ihi,jlo:jhi)-thk(ilo:ihi,jlo:jhi)
    end if

  end subroutine removelakevostok

  subroutine fixthwaitesshelf(thk,topg,umod,x,y,nx,ny)
    !bedmap2 has thaites shelf free-floating, but
    !the rignot velocity data has it grounded at the tip
    !Also, remove a hangnail

    integer, intent(in) :: nx,ny
    real(kind=8), dimension(1:nx,1:ny), intent(inout) :: thk,topg,umod
    real(kind=8), dimension(1:nx), intent(in) :: x
    real(kind=8), dimension(1:ny), intent(in) :: y

    integer :: i,j,ilo,jlo,ihi,jhi,ncells
    real(kind=8) :: dx,r

    dx = x(2) - x(1)
    r = rhoo/rhoi

    !ground tip
    !approximate bottom left corner of region to be modified    
    call computebox(ncells,ilo,ihi,jlo,jhi,-1.598d+6,-0.458d+6,-1.586d+6,-0.446d+6,x,y,nx,ny)
    if (ncells.gt.0) then

       do j = jlo,jhi
          do i = ilo,ihi
             if ( (umod(i,j).gt.1.0d0) .and. (umod(i,j).lt.300.0d0)) then
                if (thk(i,j).gt.1.0d0) then
                   if (topg(i,j).gt.-800.0d0) then
                      topg(i,j) = max(topg(i,j),(-300.0d0*0.750d0 + topg(i,j)*0.25))
                      thk(i,j) = max(thk(i,j),-r*topg(i,j) + 20.0d0)
                   end if
                end if
             end if
          end do
       end do
    end if

    !hangnails
    call computebox(ncells,ilo,ihi,jlo,jhi,-1.565d+6,-0.512d+6,-1.545d+6,-0.502d+6,x,y,nx,ny)
    if (ncells.gt.0) then
       write(*,*) ncells,ilo,ihi,jlo,jhi
       thk(ilo:ihi,jlo:jhi) = 0.0d0
    end if

    call computebox(ncells,ilo,ihi,jlo,jhi,-1.602d+6,-0.479d+6,-1.597d+6,-0.473d+6,x,y,nx,ny)
    if (ncells.gt.0) then
       thk(ilo:ihi,jlo:jhi) = 0.0d0
    end if

  end subroutine fixthwaitesshelf

end module bmmods

subroutine prep()
     
  use ncio
  use physparam
  use globaldata
  use bmmods

  implicit none
  real(kind=8), dimension(1:nxin,1:nyin) :: tmp !temporary, used for reading data which might be large
  real(kind=8), dimension(1:nxin) :: xt
  real(kind=8), dimension(1:nyin) :: yt

  real(kind=8), dimension(1:nxin4,1:nyin4) :: tmp4 !temporary, used for reading data which might be large
  real(kind=8), dimension(1:nxin4) :: xt4
  real(kind=8), dimension(1:nyin4) :: yt4


  real(kind=8), dimension(1:nx,1:ny) :: thk, topg, usrf, mask, umod, umodg, umodc, sec, btrc
  real(kind=8), dimension(1:nx) :: x
  real(kind=8), dimension(1:ny) :: y

  real(kind=8), dimension(1:nx4,1:ny4) :: temperature4
  real(kind=8), dimension(1:nx4) :: x4
  real(kind=8), dimension(1:ny4) :: y4

  character(len=12) :: tempname
  integer i,j,ixhi,iyhi,ixhi4,iyhi4
  real(kind=8) :: r
  
  !1km resolution data (geoemtry, speed)
  ixhi = ixlo-1+nx
  iyhi = iylo-1+ny
  
  call ncloadone(xt,yt,tmp,invelfile,"umod",nxin,nyin)
  x = xt(ixlo:ixhi)
  y = yt(iylo:iyhi) 
  umod = tmp(ixlo:ixhi, iylo:iyhi)

  call ncloadone(xt,yt,tmp,ingeofile,"thk",nxin,nyin)
  thk = tmp(ixlo:ixhi, iylo:iyhi)
  
  call ncloadone(xt,yt,tmp,ingeofile,"topg",nxin,nyin)
  topg = tmp(ixlo:ixhi, iylo:iyhi)

  call ncloadone(xt,yt,tmp,ingeofile,"usrf",nxin,nyin)
  usrf = tmp(ixlo:ixhi, iylo:iyhi)

  call ncloadone(xt,yt,tmp,ingeofile,"mask",nxin,nyin)
  mask = tmp(ixlo:ixhi, iylo:iyhi)

  !ensure thk + topg <= usrf... 
  where (thk + topg .gt. usrf + cm)
     thk = usrf - topg
  end where
   
  !try to preserve mask, but limit change in surface
  where ( (mask.gt.0.99) .and. ( abs ( usrf - thk*(1.0d0-rhoi/rhoo)) .le. 30.0))
     topg = min(topg, -20.0 - rhoo/rhoi*thk)
  end where


  !remove lake vostok
  call removelakevostok(thk,topg,usrf,x,y,nx,ny)

  !remove the thin tongues that annoy the solver so much
  call removethintongues(thk,topg,umod,x,y,nx,ny)

  !thaites western (slow) shelf needs to be grounded
  !at its tip.
  call fixthwaitesshelf(thk,topg,umod,x,y,nx,ny)

  

  call ncloadone(xt,yt,tmp,insecfile,"smask",nxin,nyin)
  sec = tmp(ixlo:ixhi, iylo:iyhi)

  !grow the velocity into missing regions befoire estimating corigin
  call growposdef(umodg,umod,nx,ny,1.0d0,100)
  call computebtrc(topg, thk, umodg, btrc, x(2)-x(1))

  do i = 1,32
     where (sec .eq. qsectors(i))
        umod = 0.0d0
        btrc = csticky
     end where

     !remove ice shelves in quiescent regions
     r = (1.0d0 - 918.0d0/1028.d0)
     where ( (sec .eq. qsectors(i)) .and. ( (thk*r).gt.(thk+topg)))
      thk = 0.0d0
     end where
  end do

  !remove any ice that isn't actually part of the AIS
  call removeremoteice(thk,topg,umod,x,y,nx,ny)

  call ncsaveone(x,y,umod,nx,ny,outctrlfile,"umod")
  call ncaddone(x,y,btrc,nx,ny,outctrlfile,"btrc")


  !velocity coef in objective function (1/ (2*sigma^2))
  umodc = 1.0d0 ! assume error is uniform
  where (umod.lt.1.0d0)
     umodc = 0.0d0
  end where
  where (thk.lt.1.0d0)
     umodc = 0.0d0
  end where
  !grow the no-confidence region, avoid fitting edge artifacts in umod
  call growmin(umodc,nx,ny,2)
  call ncaddone(x,y,umodc,nx,ny,outctrlfile,"umodc")

  call ncsaveone(x,y,thk,nx,ny,outgeofile,"thk")
  call ncaddone(x,y,topg,nx,ny,outgeofile,"topg")

  !4km resolution data (temperature)
  ixhi4 = ixlo4-1+nx4
  iyhi4 = iylo4-1+ny4
  
  !not much to do here, just subset the temperature
  do i = 0,nsigma-1
     write(tempname,'("temp",i6.6)') i
     call ncloadone(xt4,yt4,tmp4,intempfile,tempname,nxin4,nyin4)
     temperature4 = tmp4(ixlo4:ixhi4, iylo4:iyhi4)
     if (i.eq.0) then
        x4 = xt4(ixlo4:ixhi4)
        y4 = yt4(iylo4:iyhi4)
        call ncsaveone(x4,y4,temperature4,nx4,ny4,outtempfile,tempname)
     else
        call ncaddone(x4,y4,temperature4,nx4,ny4,outtempfile,tempname)
     end if
  end do


  return
end subroutine prep


program prepcontrol
  use globaldata
  implicit none
 

  character(len=512) :: exename, nmlfile
 

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
  read(8,nml=files)
  !dimensions
  nxin = 6144 !default to bedmap2 sizes
  nyin = 6144
  ixlo = 1
  iylo = 1
  nx = nxin
  ny = nyin
  nsigma = 10
  read(8,nml=dims)


  if (mod(nxin,4).eq.0) then
     nxin4 = nxin/4
  else
     write (*,*) 'nxin must be a factor of 4 but nxin = ',nxin
     stop 1
  end if


  if (mod(nyin,4).eq.0) then
     nyin4 = nyin/4
  else
     write (*,*) 'nyin must be a factor of 4 but nyin = ',nyin
     stop 1
  end if


  if (mod(nx,4).eq.0) then
     nx4 = nx/4
  else
     write (*,*) 'nx must be a factor of 4 but nx = ',nx 
     stop 1
  end if


  if (mod(ny,4).eq.0) then
     ny4 = ny/4
  else
     write (*,*) 'ny must be a factor of 4 but ny = ',ny 
     stop 1
  end if

  if (mod(ixlo-1,4).eq.0) then
     ixlo4 = (ixlo-1)/4 + 1
  else
     write (*,*) 'ixlo - 1 must be a factor of 4 but ixlo = ',ixlo 
     stop 1
  end if

  if (mod(iylo-1,4).eq.0) then
     iylo4 = (iylo-1)/4 + 1
  else
     write (*,*) 'iylo - 1 must be a factor of 4 but iylo = ',iylo 
     stop 1
  end if




  !quescient sectors
  qsectors = -1
  read(8,nml=quiescent)

  close(8)
 
  call prep()

end program prepcontrol
 
