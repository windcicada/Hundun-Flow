! SPDX-License-Identifier: Apache-2.0
! Serial driver for the complete, unchanged reference courant.F90 routine.
module global
 integer,parameter :: lower=1,upper=32,l=2,m=2,n=2
 integer :: nv
 real :: dtim
end module
module arrays
 real :: ajc(32),b11(32),b22(32),b33(32),rho(32),cou(3),diffnum(3),pec(3)
 integer :: jo(3)=[0,3,6],ko(3)=[0,9,18]
end module
module exchange
 integer,parameter :: mpi_in_place=0,mpi_real=1,mpi_max=2,mpi_comm_world=3
 integer :: info
contains
 subroutine mpi_allreduce(in_place,value,count,datatype,operation,comm,error)
  integer,intent(in)::in_place,count,datatype,operation,comm
  real,intent(inout)::value
  integer,intent(out)::error
  error=0 ! MPI_COMM_SELF identity; this driver compares the local operator.
 end subroutine
end module
program probe
 use global
 use arrays
 implicit none
 real::r,v,flux(6),gi(32),gj(32),gk(32),gam(32)
 integer::ios
 do
  read(*,*,iostat=ios)r,v,dtim,flux
  if(ios/=0)exit
  rho=r;ajc=v;b11=0;b22=0;b33=0;gam=0;gi=0;gj=0;gk=0
  gi(14)=flux(1);gi(15)=flux(2)
  gj(14)=flux(3);gj(17)=flux(4)
  gk(14)=flux(5);gk(23)=flux(6)
  call courant(gi,gj,gk,gam)
  write(*,'(es25.17)')maxval(cou)
 enddo
end program
