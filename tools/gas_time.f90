! SPDX-License-Identifier: Apache-2.0
! Serial driver for complete, unchanged cmod.F90 and step.F90.
module global
 integer,parameter :: lower=1,upper=32,l=2,m=2,n=2
 integer,parameter :: pc=1,sc=2,nc=3,wc=4,ec=5,lc=6,rc=7,bpc=8,nvdp=2
 integer :: nv
 real :: dtim
end module
module arrays
 real :: coef(8,32),fold(64),drhodt(32)
 integer :: jo(3)=[0,3,6],ko(3)=[0,9,18],nfo(2)=[0,32]
end module
module exchange
end module
module sgs_pdf
end module
program probe
 use global
 use arrays
 implicit none
 real :: density(32),storage,old(7),spatial(8),rate
 integer :: ios,mode
 do
  read(*,*,iostat=ios) mode,storage,dtim,old,spatial,rate
  if(ios/=0) exit
  coef=0;fold=0;density=storage;drhodt=rate;nv=1
  coef(:,14)=spatial
  fold(14)=old(1);fold(13)=old(2);fold(15)=old(3)
  fold(11)=old(4);fold(17)=old(5);fold(5)=old(6);fold(23)=old(7)
  if(mode==1) call cmod
  if(mode==3) nv=nvdp
  call step(density)
  write(*,'(8(es25.17,1x))')coef(:,14)
 enddo
end program
