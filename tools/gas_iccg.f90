! SPDX-License-Identifier: Apache-2.0
! Single-rank surroundings for the complete unchanged cgsol.F90 (REAL4).
module global
 use mpi
 implicit none
 integer,parameter :: pc=1,sc=2,nc=3,wc=4,ec=5,lc=6,rc=7,bpc=8,nv=1
 integer :: lower=1,upper,l,m,n,nijk,info,mout=0,scrn=0
 real :: dvol
 real,parameter :: big=1e30,small=1e-30
 logical,parameter :: master=.true.
 logical :: periodic(3)
end module
module arrays
 use global,only: nijk,lower,upper
 implicit none
 real,allocatable,target :: work(:),f(:)
 real,allocatable :: ajc(:),field(:),coef(:,:)
 real :: resid(1),rnorm(1),tol(1),solver_error(1)
 integer,allocatable :: io(:),jo(:),ko(:)
 integer :: icycl(1),nfo(4),ncycl(1)
end module
module exchange
 use global
 implicit none
 contains
 function assign_pointer(value,lo,hi) result(ptr)
   integer,intent(in) :: lo,hi
   real,target,intent(inout) :: value(lo:hi)
   real,pointer :: ptr(:)
   ptr=>value
 end function
 subroutine pbsrhl(values,mode)
   real,intent(inout) :: values(lower:upper)
   integer,intent(in) :: mode
   integer :: x,y,z,a,at(3),src(3),dims(3),i,j
   dims=[l+1,m+1,n+1]
   do z=1,dims(3)
    do y=1,dims(2)
     do x=1,dims(1)
      at=[x,y,z];src=at
      if(all(at>=2).and.all(at<=dims-1)) cycle
      i=x+dims(1)*((y-1)+dims(2)*(z-1))
      values(i)=0
      do a=1,3
       if(at(a)==1) then
        if(.not.periodic(a)) exit
        src(a)=dims(a)-1
       endif
       if(at(a)==dims(a)) then
        if(.not.periodic(a)) exit
        src(a)=2
       endif
      enddo
      if(a<=3) cycle
      j=src(1)+dims(1)*((src(2)-1)+dims(2)*(src(3)-1))
      values(i)=values(j)
     enddo
    enddo
   enddo
 end subroutine
end module
module extras
end module
subroutine boffin_stop(filename,line)
 character(*),intent(in) :: filename
 integer,intent(in) :: line
 write(0,*) 'reference stop ',filename,line
 stop 3
end subroutine
program probe
 use global
 use arrays
 implicit none
 integer :: nx,ny,nz,x,y,z,i,ios,mask(3)
 call MPI_Init(info)
 do
  read(*,*,iostat=ios) nx,ny,nz,mask
  if(ios/=0) exit
  l=nx+1;m=ny+1;n=nz+1;nijk=(nx+2)*(ny+2)*(nz+2);upper=nijk
  periodic=mask/=0
  allocate(work(4*nijk),f(nijk),ajc(nijk),field(nijk),coef(8,nijk))
  allocate(io(nx+2),jo(ny+2),ko(nz+2))
  work=0;f=0;field=0;coef=0;ajc=1
  io=[(x-1,x=1,nx+2)]
  jo=[((y-1)*(nx+2),y=1,ny+2)]
  ko=[((z-1)*(nx+2)*(ny+2),z=1,nz+2)]
  nfo=[0,nijk,2*nijk,3*nijk];rnorm=1;tol=2e-6;icycl=2000
  do z=2,n
   do y=2,m
    do x=2,l
     i=x+jo(y)+ko(z)
     read(*,*) ajc(i),coef(pc,i),coef(sc,i),coef(nc,i),coef(wc,i), &
       coef(ec,i),coef(lc,i),coef(rc,i),coef(bpc,i)
    enddo
   enddo
  enddo
  call cgsol
  write(*,'(i8,1x,es25.17)') ncycl(1),solver_error(1)
  do z=2,n
   do y=2,m
    do x=2,l
     i=x+jo(y)+ko(z)
     write(*,'(10(es25.17,1x))') f(i),work(nijk+i),coef(:,i)
    enddo
   enddo
  enddo
  deallocate(work,f,ajc,field,coef,io,jo,ko)
 enddo
 call MPI_Finalize(info)
end program
