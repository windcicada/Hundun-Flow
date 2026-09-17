! SPDX-License-Identifier: Apache-2.0
! Serial surrounding driver for complete, unchanged statistics.F90.
! The serial reduction and Cartesian gradient are fixture infrastructure;
! statistics, Dynamic_Cphi and test_avg execute their original full bodies.
module global
 implicit none
 integer,parameter :: lower=1,upper=64,l=3,m=3,n=3,lp1=4,mp1=4,np1=4
 integer,parameter :: nvu=1,nvv=2,nvw=3,nvf=4,nf=4,nfield=2
 integer :: nv,dyn_C_count=0,ph=1,phmax=1,mout=6,scrn=99,info
 real :: dtim,atime=0,tim=0,small=1e-30,dyn_C_time=0,frequency=1
 logical :: master=.false.,pdf_kappa=.true.,phase_averaging=.false.
 character(16) :: species_output='mass_fraction'
end module
module arrays
 use global,only:lower,upper
 implicit none
 integer :: io(4)=[0,1,2,3],jo(4)=[0,4,8,12],ko(4)=[0,16,32,48],nfo(64)
 real,target :: f(4096)
 real :: fstat(4096),ftau(4096),fschem(4096),p_mean(64),p(64),x(64),y(64),z(64),ajc(64)
 real :: gam_sgs(64),rho(64),visc(64),eps(64),q(64),dfdx(64),dfdy(64),dfdz(64)
 real :: dyn_LM(3,64),dyn_M2(3,64)
 contains
 function assign_pointer(a,lo,hi) result(ptr)
 integer,intent(in)::lo,hi
 real,target,intent(inout)::a(hi-lo+1)
 real,pointer::ptr(:)
 ptr(lo:hi)=>a
 end function
end module
module chemistry
 implicit none
 integer,parameter :: nsp=6,nsc=7
 integer :: isp,jfuel=3,dstep_ww0=0,init_ww0=0
 character(8) :: names(7)=[character(8)::'O2','N2','FUEL','H2O','CO2','OH','H']
 real :: tim_sp(6,64),tim_flow(64),kappa(7,64),temp_i(64),temp_k(64),prev_rdot(6,64),arr_eta(64)
 real :: sum_w(6,64),sum_w0(6,64),fsc(7,64),sumn(64),temp(64),qdot_rad(64)
 real(kind=8) :: wm(7)=[32d0,28d0,16d0,18d0,44d0,17d0,1d0]
end module
module exchange
 implicit none
 integer,parameter :: MPI_REAL=0,mpi_max=1,mpi_comm_world=0
 contains
 subroutine mpi_allreduce(a,b,c,t,o,comm,ierr)
 real,intent(in)::a
 real,intent(out)::b
 integer,intent(in)::c,t,o,comm
 integer,intent(out)::ierr
 b=a;ierr=0
 end subroutine
end module
module extras
 real :: deta1dx,deta2dx,deta3dx,deta1dy,deta2dy,deta3dy,deta1dz,deta2dz,deta3dz
 real :: phase_average(1,4096),phase_time(1)
 integer :: phn(1)
end module
module sgs_pdf
 real(kind=8) :: rdot(6,64),rdot_mean(6,64)
 real :: field_hdot(0:2,64)
end module
program probe
 use global
 use arrays
 use chemistry
 use extras
 use sgs_pdf
 implicit none
 integer :: reset,ios,i,j,k,c,s,steps
 real :: eta,pdf(6),psr(6),composition(6),dissipation,velocity
 character(8) :: mode
 do i=1,64
   nfo(i)=(i-1)*64
 enddo
 call get_command_argument(1,mode)
 if(trim(mode)=='mix')then
   call mix_probe
   stop
 endif
 do k=1,4
 do j=1,4
 do i=1,4
   c=i+jo(j)+ko(k);x(c)=i;y(c)=j;z(c)=k
 enddo
 enddo
 enddo
 do
  read(*,*,iostat=ios) reset,dtim,eta,dissipation,velocity,pdf,psr,composition
  if(ios/=0)exit
  if(reset/=0)then
   atime=0;tim=0;dyn_C_time=0;dyn_C_count=0;dstep_ww0=0;init_ww0=0
   sum_w=0;sum_w0=0;kappa=1;fstat=0;ftau=0;fschem=0;p_mean=0;rdot_mean=0
   dyn_LM=2;dyn_M2=0;tim_sp=0;tim_flow=0;steps=0
  endif
  ajc=1;rho=1;visc=1e-4;eps=dissipation;q=1;gam_sgs=0;p=101325
  temp=600;qdot_rad=0;field_hdot=0;f=0;f(1:64)=velocity
  arr_eta=eta;sumn=sum(composition);fsc(7,:)=0
  do s=1,6
   rdot(s,:)=pdf(s);prev_rdot(s,:)=psr(s);fsc(s,:)=composition(s)
   f(nfo(nf+s)+1:nfo(nf+s)+64)=composition(s)
  enddo
  call statistics
  steps=steps+1;tim=tim+dtim;c=2+jo(2)+ko(2)
  write(*,'(4(i0,1x),32(es25.17,1x))')steps,dstep_ww0,init_ww0,dyn_C_count, &
   kappa(1:6,c),sum_w(:,c),sum_w0(:,c),tim_sp(:,c),tim_flow(c),dyn_LM(1,c)
 enddo
 contains
 subroutine mix_probe
 integer :: status,cell,ii,jj,kk
 real :: scalar(64)
 ! Isolate the complete source's mixture-fraction channel. Other channels
 ! start with zero products, so their invalid-ratio markers remain -1.
 names(6)='TRACE';jfuel=0
 do
   read(*,*,iostat=status) rho,ajc,scalar
   if(status/=0)exit
   f=0;f(nfo(nvf)+1:nfo(nvf)+64)=scalar
   fsc=0;dyn_LM=0;dyn_M2=0;temp_i=0;temp_k=0;dtim=1;dyn_C_time=0
   call Dynamic_Cphi
   do kk=2,3
   do jj=2,3
   do ii=2,3
     cell=ii+jo(jj)+ko(kk)
     write(*,'(2(es25.17,1x))')dyn_M2(1,cell),dyn_LM(1,cell)
   enddo
   enddo
   enddo
 enddo
 end subroutine
end program
subroutine pbsrhl(a,components)
 use global,only:lower,upper
 implicit none
 integer,intent(in)::components
 real,intent(inout)::a(lower:upper,components)
 ! Fixture inputs include complete constant face/corner states. Statistics
 ! does not require remote donors in this single-block algebra/cadence probe.
end subroutine
subroutine gradient(a,gx,gy,gz)
 use arrays,only:jo,ko
 implicit none
 real,intent(in)::a(64)
 real,intent(out)::gx(64),gy(64),gz(64)
 integer::i,j,k,c,ia,ib,ja,jb,ka,kb
 do k=1,4
 do j=1,4
 do i=1,4
  ia=max(1,i-1);ib=min(4,i+1);ja=max(1,j-1);jb=min(4,j+1);ka=max(1,k-1);kb=min(4,k+1)
  c=i+jo(j)+ko(k)
  gx(c)=(a(ib+jo(j)+ko(k))-a(ia+jo(j)+ko(k)))/(ib-ia)
  gy(c)=(a(i+jo(jb)+ko(k))-a(i+jo(ja)+ko(k)))/(jb-ja)
  gz(c)=(a(i+jo(j)+ko(kb))-a(i+jo(j)+ko(ka)))/(kb-ka)
 enddo
 enddo
 enddo
end subroutine
