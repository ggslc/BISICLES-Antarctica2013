module ncio
  use netcdf
  implicit none
contains
  subroutine nccheck(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop 2
    end if
  end subroutine nccheck

  subroutine ncaddone(x,y,z,n,m,file,field)
    
    character(len=*),intent(in) :: file,field
    integer,intent(in) :: n,m
    real(kind=8), dimension(1:n) ,intent(in):: x
    real(kind=8), dimension(1:m) ,intent(in):: y
    real(kind=8), dimension(1:n,1:m) ,intent(in):: z
    integer i,j,nc_id, var_id, cell_dim_id(2)

    call nccheck( nf90_open(file,  NF90_WRITE, nc_id) )
    call nccheck( nf90_redef(nc_id) )
    call nccheck(nf90_inq_dimid(nc_id,"x",cell_dim_id(1)))
    call nccheck(nf90_inq_dimid(nc_id,"y",cell_dim_id(2)))

    !cell centered array defintions
    call nccheck( nf90_def_var(nc_id, field, nf90_real8, cell_dim_id, var_id) )

    !end of definition section
    call nccheck( nf90_enddef(nc_id) )

    !write cell centred arrays
    call nccheck( nf90_inq_varid(nc_id, field, var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , z) )
    call nccheck( nf90_close(nc_id) )

    return
  end subroutine ncaddone

  subroutine ncsaveone(x,y,z,n,m,file,field)

    character(len=*),intent(in) :: file,field
    integer,intent(in) :: n,m
    real(kind=8), dimension(1:n) ,intent(in):: x
    real(kind=8), dimension(1:m) ,intent(in):: y
    real(kind=8), dimension(1:n,1:m) ,intent(in):: z
    integer i,j,nc_id, var_id, cell_dim_id(2)

    call nccheck( nf90_create(file,  NF90_CLOBBER, nc_id) )
    
    !cell centred dimensions
    call nccheck( nf90_def_dim(nc_id, "x", n , cell_dim_id(1)) )
    call nccheck( nf90_def_dim(nc_id, "y", m , cell_dim_id(2)) )
    

    !dimension definitions
    call nccheck( nf90_def_var(nc_id, "x", nf90_real8, cell_dim_id(1), var_id))
    call nccheck( nf90_put_att(nc_id, var_id, "units", "m") )
    call nccheck( nf90_def_var(nc_id, "y", nf90_real8, cell_dim_id(2), var_id))
    call nccheck( nf90_put_att(nc_id, var_id, "units", "m") )
    
    !cell centered array defintions
    call nccheck( nf90_def_var(nc_id, field, nf90_real8, cell_dim_id, var_id) )

    !end of definition section
    call nccheck( nf90_enddef(nc_id) )

    !cell centred x and y
    call nccheck( nf90_inq_varid(nc_id, "x", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , x) )
    call nccheck( nf90_inq_varid(nc_id, "y", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , y) )

    !cell centred arrays
    call nccheck( nf90_inq_varid(nc_id, field, var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , z) )

    call nccheck( nf90_close(nc_id) )

  end subroutine ncsaveone

  subroutine ncloadone(x,y,a,file,field,ewn,nsn)
    integer ewn,nsn
    real(kind=8), dimension(1:ewn) :: x
    real(kind=8), dimension(1:nsn) :: y
    real(kind=8), dimension(1:ewn,1:nsn) :: a
   
    character(len=*),intent(in) :: file, field
    integer  var_id, nc_id
    
    call nccheck( nf90_open(file, NF90_NOWRITE, nc_id) )
    
    call nccheck( nf90_inq_varid(nc_id, "x", var_id) )
    call nccheck( nf90_get_var(nc_id, var_id , x) )
    call nccheck( nf90_inq_varid(nc_id, "y", var_id) )
    call nccheck( nf90_get_var(nc_id, var_id , y) )
    call nccheck( nf90_inq_varid(nc_id, field, var_id) )
    call nccheck( nf90_get_var(nc_id, var_id , a) )
    call nccheck( nf90_close(nc_id) )
 
    return
  end subroutine ncloadone


  subroutine ncsaven(x,y,c,ewn,nsn,nc,file,field)
    !write a 2D netcdf file 'file' with 1 <= ic <= nc fields named field(ic)
    !containing data f(x,y) stored in c(:,:,ic)

    character(len=*), intent(in) :: file
    character(len=*), dimension(nc), intent(in) :: field
    integer,intent(in) :: ewn,nsn,nc
    real(kind=8), dimension(1:ewn) ,intent(in):: x
    real(kind=8), dimension(1:nsn) ,intent(in):: y
    real(kind=8), dimension(1:ewn,1:nsn, 1:nc) ,intent(in):: c
    integer i,j,nc_id, var_id, cell_dim_id(2), ic

    call nccheck( nf90_create(file,  NF90_CLOBBER, nc_id) )
    
    !cell centred dimensions
    call nccheck( nf90_def_dim(nc_id, "x", ewn , cell_dim_id(1)) )
    call nccheck( nf90_def_dim(nc_id, "y", nsn , cell_dim_id(2)) )
    
    !dimension definitions
    call nccheck( nf90_def_var(nc_id, "x", nf90_real8, cell_dim_id(1), var_id))
    call nccheck( nf90_put_att(nc_id, var_id, "units", "m") )
    call nccheck( nf90_def_var(nc_id, "y", nf90_real8, cell_dim_id(2), var_id))
    call nccheck( nf90_put_att(nc_id, var_id, "units", "m") )
    
    !cell centered array defintions
    do ic = 1,nc
       call nccheck( nf90_def_var(nc_id, field(ic), nf90_real8, cell_dim_id, var_id) )
    end do
    !end of definition section
    call nccheck( nf90_enddef(nc_id) )

    !cell centred x and y
    call nccheck( nf90_inq_varid(nc_id, "x", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , x) )
    call nccheck( nf90_inq_varid(nc_id, "y", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , y) )

    !cell centred arrays
    do ic = 1,nc
       call nccheck( nf90_inq_varid(nc_id, field(ic), var_id) )
       call nccheck( nf90_put_var(nc_id, var_id , c(:,:,ic)) )
    end do
    call nccheck( nf90_close(nc_id) )

  end subroutine ncsaven


  subroutine ncloadonenoxy(a,file,field,ewn,nsn)
    integer ewn,nsn
    
    real(kind=8), dimension(1:ewn,1:nsn) :: a
   
    character(len=*),intent(in) :: file, field
    integer  var_id, nc_id
    
    call nccheck( nf90_open(file, NF90_NOWRITE, nc_id) )
    

    call nccheck( nf90_inq_varid(nc_id, field, var_id) )
    call nccheck( nf90_get_var(nc_id, var_id , a) )
    call nccheck( nf90_close(nc_id) )
 
    return
  end subroutine ncloadonenoxy

end module ncio
