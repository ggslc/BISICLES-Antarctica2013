

subroutine inject(b,xb,yb,ewnb,nsnb,a,xa,ya,ewna,nsna)
 
  integer :: ewna , nsna
  real(kind=8), dimension(1:ewna,1:nsna) :: a
  real(kind=8), dimension(1:ewna) :: xa
  real(kind=8), dimension(1:nsna) :: ya
  integer :: ewnb , nsnb
  real(kind=8), dimension(1:ewnb,1:nsnb) :: b
  real(kind=8), dimension(1:ewnb) :: xb
  real(kind=8), dimension(1:nsnb) :: yb
 
  real(kind=8) :: dxa,dxb
  integer ia,ja,ib,jb,ioff,joff,ioff0,joff0

  dxa = xa(2)-xa(1)
  dxb = xb(2)-xb(1)
  
  ioff0 = 1
  do while (xb(ioff0) .lt. (xa(1) - 0.5*dxa))
     ioff0 = ioff0 + 1
  end do

  joff0 = 1
  do while (yb(joff0) .lt. (ya(1) - 0.5*dxa))
     joff0 = joff0 + 1
  end do


  ioff = ioff0
  

  do ia = 1,ewna
     if (xa(ia).ge.(xb(ioff)-dxa)) then
        do ib = ioff,ioff+5
           if (ib.le.ewnb) then
              joff = joff0
              do ja = 1,nsna
                 if (ya(ja).ge.(yb(joff)-dxa)) then
                    do jb = joff,joff+5
                       if (jb.le.nsnb) then  
                          b(ib, jb) = a(ia,ja)
                       end if
                    end do
                    joff = joff + 5
                 end if
              end do
           end if
        end do
        ioff = ioff + 5
     end if
     
  end do



end subroutine inject

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

end subroutine smooth

module tempmod
   implicit none
   real (kind = 8),parameter ::  maxtemp = 272.0d0, mintemp = 200.0d0, seatemp = 272.0d0
 contains
   
subroutine growtemp(temp,thk,ewn,nsn,niter)
  implicit none
  integer :: ewn,nsn,upn,niter,iter,ew,ns,ww,we,ws,wn,wsum
  real(kind=8), dimension(1:ewn,1:nsn) :: thka,tempa,thk,temp,thkb
 
  
  where (.not.(thk.ge.mintemp))
     temp = seatemp
  end where

  where (.not.(temp.ge.mintemp))
     temp = mintemp
  end where
  thkb = thk

  do iter = 1,niter
     write(*,*) 'growtemp iteration ',iter
     tempa = temp
     thka = thkb
     do ns = 2,nsn-1
        do ew = 2,ewn-1
           if (thk(ew,ns).le.1.0d0) then
            
              
  
              ww = thkb(ew-1,ns)
              we = thkb(ew+1,ns)
              ws = thkb(ew,ns-1)
              wn = thkb(ew,ns+1)

              wsum = ww + we + ws + wn
              if (wsum.gt.1.0d-6) then
                 tempa(ew,ns) = (1.0d0/wsum) * &
                      ( ww*temp(ew-1,ns) + we*temp(ew+1,ns) &
                      +  ws*temp(ew,ns-1) + wn*temp(ew,ns+1))
              end if

              thka(ew,ns) = 0.25 * &
                   ( thkb(ew-1,ns) + thkb(ew+1,ns) &
                  +  thkb(ew,ns-1) + thkb(ew,ns+1))
              
           end if

        end do
     end do
    
     temp = tempa
     thkb=thka
  end do
  return
end subroutine growtemp


subroutine sanetemp(temp,ewn,nsn,upn,tlambda,dx)
  implicit none
  integer :: ewn, nsn, upn, i
  real (kind = 8) :: tlambda,dx
  real (kind = 8), dimension(1:ewn,1:nsn,1:upn) :: temp
  

  where (.not.(temp.gt.mintemp))
     temp = seatemp
  end where


  do i = 2,upn
     where (temp(:,:,i).lt.temp(:,:,i-1))
        temp(:,:,i) = temp(:,:,i-1)
     end where
     where (temp(:,:,i).gt.(5.0 + temp(:,:,i-1)))
        temp(:,:,i) = temp(:,:,i-1) + 5.0
     end where
  end do

 
  where (temp .gt. maxtemp)
     temp = maxtemp
  end where

  do i = 1,upn
     write(*,*) "smoothing temperature in layer ", i
     call smooth(temp(:,:,i),tlambda,dx,ewn,nsn,5)
  end do
end subroutine sanetemp
end module tempmod

program albtobed
  use ncio
  use tempmod
  implicit none
  
  character(len=64)  albmap
  
  integer, parameter :: upn = 10
  character(len=10), dimension(upn) :: tempnames
  character(len=10) :: tempname

  !ALBMAP (5km sized) data
  integer, parameter :: ewna = 1280, nsna = 1280
  real(kind=8), dimension(1:ewna,1:nsna) :: tmpa, thka
  real(kind=8), dimension(1:ewna) :: xa,xat
  real(kind=8), dimension(1:nsna) :: ya,yat


  real(kind=8), dimension(1:1160,1:1120) :: smasks !sectors ahve some missing data relative to the 1280*1280 grid

  !BEDMAP2 (1km sized) data
  integer, parameter :: ewn = 6144, nsn = 6144
  real(kind=8), dimension(1:ewn,1:nsn) :: tmp
  real(kind=8), dimension(1:ewn) :: x
  real(kind=8), dimension(1:nsn) :: y

  !BEDMAP2 coarsened (2km sized) data
  integer, parameter :: ewnc = 3072, nsnc = 3072
  real(kind=8), dimension(1:ewnc,1:nsnc) :: tmpc
  real(kind=8), dimension(1:ewnc) :: xc
  real(kind=8), dimension(1:nsnc) :: yc

  !BEDMAP2 twice coarsened (4km sized) data
  integer, parameter :: ewncc = 1536, nsncc = 1536
  real(kind=8), dimension(1:ewncc,1:nsncc) :: tmpcc
  real(kind=8), dimension(1:ewncc) :: xcc
  real(kind=8), dimension(1:nsncc) :: ycc
  real(kind=8), dimension(1:ewncc,1:nsncc,1:upn) :: tempcc

  integer, parameter :: ewmin= 263 ,nsmin = 263 !start of 6144*6144 BEDMAP2
  !  data relative to base BEDMAP2 6777*6777 data

  integer i,j,k
  real(kind=8) dxa,dx,dxc,dxcc
  dxa = 5.0d+3 ! ALBMAP resolution
  dx = 1.0d+3  ! BEDMAP2 resolution
  dxc = 2.0d0 * dx
  dxcc = 2.0d0 * dxc

  !cell centers on 1km grid
  x(1) = -3333500.0d0 + real(ewmin-1)*dx
  y(1) = -3333500.0d0 + real(nsmin-1)*dx

  do i = 2,ewn
     x(i) = x(i-1) + dx
  end do
  
  do i = 2,nsn
     y(i) = y(i-1) + dx
  end do
  
  

  !albmap input file, copied from Antarctica2012 work. HAS 
  albmap = "ALBMAP_i2s_4BISICLES.nc"
  

  !Frank Pattyn temperature data on ALBMAP  
  do i = 1,upn
     if (i .lt.10) then
        write(tempname,'("temp",i1)') i
     else
        write(tempname,'("temp",i2)') i
     end if
     call ncloadone(xa,ya,tmpa,albmap,tempname,ewna,nsna)
     if (i .eq. 1) then
        call ncloadone(xa,ya,thka,albmap,"thick",ewna,nsna)
     end if

     call growtemp(tmpa,thka,ewna,nsna,20)


     call inject(tmp,x,y,ewn,nsn,tmpa,xa,ya,ewna,nsna)

     if (i.eq.1) then
        call ncsaveone(xa,ya,tmpa,ewna,nsna,"test.nc","stemp")
     end if

     !coarsen(nf,xf,yf,nc,xc,yc,af,ac,ncomp)
     call coarsen( ewn, x, y, ewnc, xc, yc, tmp, tmpc, 1)
     call coarsen( ewnc, xc, yc, ewncc, xcc, ycc, tmpc, tmpcc, 1)
     tempcc(:,:,i) = tmpcc

     write(tempnames(i),'("temp",i6.6)') i-1
  end do

  call sanetemp(tempcc,ewncc,nsncc,upn,dxcc*4.0,dxcc)
  
  call ncsaven(xcc,ycc,tempcc,ewncc,nsncc,upn,"Antarctica-temperature-4km.nc",tempnames)
  
  !sector data on ALBMAP 
  call ncloadonenoxy(smasks,'sectormask_highres_v3.nc','sectors',1160,1120)
  tmpa = 0.0
  tmpa(1:1160,1:1120) = smasks
  call inject(tmp,x,y,ewn,nsn,tmpa,xa,ya,ewna,nsna)
  call ncsaveone(x,y,tmp,ewn,nsn,"Antarctica_sectors-1km.nc","smask")


end program albtobed
