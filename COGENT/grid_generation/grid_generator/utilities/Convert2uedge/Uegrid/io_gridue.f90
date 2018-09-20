subroutine read_gridue(&
     nxm,nym,ixpt1,ixpt2,iysptrx1,&
     rm,zm,psi,br,bz,bpol,bphi,b)

  implicit none

  character*80 fname
  character*80 runidg
  integer*4 nxm, nym, ixpt1,ixpt2,iysptrx1
  double precision, dimension(nxm+2,nym+2,5), intent(IN OUT) :: rm,zm,psi,br,bz,bpol,bphi,b

  INTEGER*4 nunit, ix, iy, n
  INTEGER*4 iunit, ios

  fname="gridue"
  !!runidg="iogridue"
  !!print *, "In Fortran read_gridue..."


  iunit=1
  !!open (iunit, file='gridue', form='formatted', iostat=ios, status='read')
  open (iunit, file=fname, form='formatted', iostat=ios)

  if (ios .ne. 0) then
     print *, "**** Cannot open ", fname
     STOP
  else

     read(iunit,1999) nxm,nym,ixpt1,ixpt2,iysptrx1
     read(iunit,2000)


     !!allocate(rm(1:nxm+2,1:nym+2,5))
     read(iunit,2001) (((rm(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     print *, "done with rm..."
     read(iunit,2000)


     !!allocate(zm(1:nxm+2,1:nym+2,5))
     read(iunit,2001) (((zm(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     print *, "done with zm..."
     read(iunit,2000)


     !!allocate(psi(1:nxm+2,1:nym+2,5))
     read(iunit,2001) (((psi(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     print *, "done with psi..."
     read(iunit,2000)


     !!allocate(br(1:nxm+2,1:nym+2,5))
     read(iunit,2001) (((br(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     print *, "done with br..."
     read(iunit,2000)


     !!allocate(bz(1:nxm+2,1:nym+2,5))
     read(iunit,2001) (((bz(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     print *, "done with bz..."
     read(iunit,2000)


     !!allocate(bpol(1:nxm+2,1:nym+2,5))
     read(iunit,2001) (((bpol(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     print *, "done with bpol..."
     read(iunit,2000)


     !!allocate(bphi(1:nxm+2,1:nym+2,5))
     read(iunit,2001) (((bphi(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     print *, "done with bphi..."
     read(iunit,2000)


     !!allocate(b(1:nxm+2,1:nym+2,5))
     read(iunit,2001) (((b(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     print *, "done with b..."
     read(iunit,2002) runidg

     close (iunit)

     !!print *, 'Reading file ', fname, ' with runidg:  ', runidg
     !!print *, "file ", fname
     write(*,*) 'Reading file: ', fname, ' with runidg: ', runidg
     !!write(*,*)
          
  endif
     


1999 format(5i4)  
2000 format()

     !!  ifelse([WORDSIZE],64,\
     !!2001 format(1p3e23.15)
     !!  ,\
     !!2001 format(1p3d23.15)
     !!  )\
     
2001 format(1p3d23.15)
2002 format(a60)
     

!
!
!
!  print *, "Exiting read_gridue..."
end subroutine read_gridue




subroutine read_gridue_dims(nxm,nym)

  implicit none

  character*80 fname
  character*80 runidg
  logical res

  INTEGER*4, intent(OUT) :: nxm, nym
  INTEGER*4 ixpt1,ixpt2,iysptrx1
  INTEGER*4 nunit, ix, iy, n
  INTEGER*4 iunit, ios

  fname="gridue"
  !!runidg="iogridue"


  print *, "In Fortran read_gridue_dims ..."


  inquire (file=fname, exist=res)
  if (res .eqv. .false.) then
     print *, 'File', fname, ' does not exist!'
     stop
  else
     print *, "" !!'File', fname, ' does exist!'
  endif


  !!open (iunit, file='gridue', form='formatted', iostat=ios)
  !!print *, "#1: In Fortran opening gridue ..., ios=", ios
  !!close (iunit)



  iunit=1
  !!open (iunit, file='gridue', form='formatted', iostat=ios, status="new")
  open (iunit, file=fname, form='formatted', iostat=ios)

  !!print *, "In Fortran opening ", fname, "..., ios=", ios

  if (ios .ne. 0) then
     print *, "**** Cannot open ", fname
     STOP
  else

     read(iunit,1999) nxm,nym,ixpt1,ixpt2,iysptrx1
     close (iunit)

     write(*,*) 'Reading file :', fname !!, ' with runidg: ', runidg
     write(*,*)
          
  endif
     


1999 format(5i4)  
2000 format()

     !!  ifelse([WORDSIZE],64,\
     !!2001 format(1p3e23.15)
     !!  ,\
     !!2001 format(1p3d23.15)
     !!  )\
     
2001 format(1p3d23.15)
2002 format(a60)
     
end subroutine read_gridue_dims



subroutine write_gridue (&
     nxm,nym,ixpt1,ixpt2,iysptrx1,&
     rm,zm,psi,br,bz,bpol,bphi,b)


  implicit none

  character*80 fname
  character*80 runidg

  INTEGER*4 nxm, nym, ixpt1,ixpt2,iysptrx1
  double precision, dimension(nxm+2,nym+2,5), intent (IN) :: rm,zm,psi,br,bz,bpol,bphi,b

  INTEGER*4 nunit, ix, iy, n
  INTEGER*4 iunit, ios

  fname="gridue"
  runidg="iogridue"


  iunit=1
  open (iunit, file=fname, form='formatted', iostat=ios, status='replace')

  if (ios .ne. 0) then
     print *, "**** Cannot open ", fname
     STOP
  else

     !print *, "In Fortran: rm(1,1,1)=", rm(1,1,1)
     !print *, "In Fortran: zm(1,1,1)=", zm(1,1,1)

     !do ix=1,nxm+2
     !   do iy=1,nym+2
     !      print *, "ix,iy, rm(ix,iy,0)", ix, iy, rm(ix,iy,1)
     !   enddo
     !enddo
     !!STOP

     write(iunit,1999) nxm,nym,ixpt1,ixpt2,iysptrx1
     write(iunit,2000)
     write(iunit,2001) (((rm(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((zm(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((psi(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((br(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((bz(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((bpol(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((bphi(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2000)
     write(iunit,2001) (((b(ix,iy,n),ix=1,nxm+2),iy=1,nym+2),n=1,5)
     write(iunit,2002) runidg

     close (iunit)

     write(*,*) 'Wrote file "', fname, '" with runidg:  ', runidg
     write(*,*)
          
  endif
     


1999 format(5i4)  
2000 format()

     !!  ifelse([WORDSIZE],64,\
     !!2001 format(1p3e23.15)
     !!  ,\
     !!2001 format(1p3d23.15)
     !!  )\
     
2001 format(1p3d23.15)
2002 format(a60)
     
end subroutine write_gridue
