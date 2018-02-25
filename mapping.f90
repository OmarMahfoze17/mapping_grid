!############################################################################
! This Code is made to map the field of INCOMPACT3D of a specific grid size 
! to a differnet grid.
! The code uses 2nd order polynomial interpolation
! 
! The code is made by Omar Mahfoze, E-mail: omar.mahfoze15@imperial.ac.uk
program mapping
!
!############################################################################

  implicit none
   
!!! ============== File 1 ===================================
  integer :: i,j,k,i1,j1,k1,j1_start,i_var
  integer :: i1_start,k1_start,dny1,dnx1
  integer(8) ::count
  integer :: istret1,ncly1
  integer, parameter :: nx1=128 , ny1=129, nz1=256,N_VAR=10 !13
  real(8),dimension(nx1, ny1, nz1) :: ux1,uy1,uz1
  real(8),dimension(:),allocatable :: y1
  real(8),dimension(nx1) :: x1
  real(8),dimension(nz1) :: z1
  real(8) ::xlx1,yly1,zlz1,dx1,dz1,t1,t2,beta1
  real(8) ::ay0,ay1,ay2,ax0,ax1,ax2,az0,az1,az2,za0
!================== File 2  ================================
  integer :: i2,j2,k2,dnx2,dny2,istret2,ncly2
  integer,parameter:: nx2=128, ny2=129, nz2=256
  real(8),dimension(nx2, ny2, nz2) :: ux2,uy2,uz2
  real(8),dimension(:),allocatable :: y2
  real(8),dimension(nx2) :: x2
  real(8),dimension(nz2) :: z2
  real(8) ::xlx2,yly2,zlz2,dx2,dz2,dy2,beta2  
  integer :: np=1 ! polynomial order  NOT USED
  logical :: exist
  !character(len=100) :: fileName1='/home/om1014/PhD/INCOMPACT3D/TBL/TBL_PRO_xlx_500/sauve.dat'
  character(len=100) :: fileName1='/home/om1014/PhD/INCOMPACT3D/Channel/channel/sauve_orig.dat'
  !character(len=100) :: fileName2='sauve_mapped_2.dat'
  character(len=100) :: fileName2='/home/om1014/PhD/INCOMPACT3D/Channel/channel/sauve_in.dat'
!!! ============== File 1 ===================================
  allocate(y1(ny1))
  xlx1=12.5663706144 !500.
  yly1=2. !40.
  zlz1=2. !40.
  istret1=2 !3
  beta1=0.259065151 !1.3
  ncly1=2
  dx1=xlx1/(nx1-1)
  dz1=zlz1/(nz1-1)
  DO I=1,nx1
     x1(i)=dx1*(i-1)
  enddo
   !open (15,file='yp1.dat',form='formatted',status='unknown')
   !do j=1,ny1
      !read(15,*) y1(j)
   !enddo
  ! close(15)
  DO K=1,nz1
     z1(k)=dz1*(k-1)
  enddo

call stretching(y1,yly1,ny1,beta1,istret1,ncly1,'yp_1.dat')  ! get Y1 values and save them in a file
!!!==========================================================
!!! ============== File 2 ===================================   
  allocate(y2(ny2))
  xlx2=xlx1
  yly2=yly1
  zlz2=zlz1
  istret2=istret1
  beta2=beta1
  ncly2=ncly1
  dx2=xlx2/(nx2-1)
  dy2=yly2/(ny2-1)
  dz2=zlz2/(nz2-1)
  DO I=1,nx2
     x2(i)=dx2*(i-1)
  enddo
  !DO j=1,ny2
     !y2(j)=dy2*(j-1)
  !enddo
  DO K=1,nz2
     z2(k)=dz2*(k-1)
  enddo

call stretching(y2,yly2,ny2,beta2,istret2,ncly2,'yp_2.dat') ! get Y2 values and save them in a file

!!!==========================================================

!! =========== TEST THE DIMENSIONS OF THE BOUNDARY ================
! The BC determines the grid numbers (odd/even), so the mappped and the original BC sould similar.
if (mod(ny2,2)/=mod(ny1,2)) THEN 
   print *, 'Dimension Mismatch NY1= ',NY1,',  NY2= ',NY2 
   stop
ENDIF
if (mod(nx2,2)/=mod(nX1,2)) THEN
   print *, 'Dimension Mismatch NX1= ',NX1,',  NX2= ',NX2 
   stop
ENDIF
!!==================================================================

 OPEN(10,FILE=fileName1,FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
DO i_var=1,N_VAR !! THIS THE OUTER LOOP TO SWETCH BETWEEN VARIABLS

if (i_var==N_VAR .and. N_VAR>1) then
   if (mod(ny1,2)/=0) then
       dny1=1; dny2=1
   endif
   if (mod(nx1,2)/=0) then
       dnx1=1; dnx2=1
    endif
else
   dny1=0; dny2=0; dnx1=0; dnx2=0
endif

print *, 'Reading the variable number : ',i_var
ux1=0.
ux2=1.
  DO K=1,nz1
     DO J=1,ny1-dny1
        DO I=1,nx1-dnx1
           READ(10,REC=COUNT) ux1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO

print *, 'COUNT = ',COUNT


call cpu_time(t1)
k1_start=2
do k2=1,nz2

j1_start=2
do j2=1,ny2-dny2
i1_start=2
do i2=1,nx2-dnx2
   do k1=k1_start,nz1
   do j1=j1_start,ny1-dny1
   do i1=i1_start,nx1-dnx1
         if (k1/=1 .and. k1/=nz1 .and. j1/=1 .and. j1/=ny1 .and. i1/=1 .and. i1/=nx1) then
         if (y1(j1-1)<=y2(j2) .and. y1(j1)>=y2(j2) .and.&
             x1(i1-1)<=x2(i2) .and. x1(i1)>=x2(i2) .and.&
             z1(k1-1)<=z2(k2) .and. z1(k1)>=z2(k2)) then

             ay0=(y2(j2)-y1(j1))*(y2(j2)-y1(j1+1))/((y1(j1-1)-y1(j1))*(y1(j1-1)-y1(j1+1))) 
             ay1=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1+1))/((y1(j1)-y1(j1-1))*(y1(j1)-y1(j1+1))) 
             ay2=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1))/((y1(j1+1)-y1(j1-1))*(y1(j1+1)-y1(j1)))

             ax0=(x2(i2)-x1(i1))*(x2(i2)-x1(i1+1))/((x1(i1-1)-x1(i1))*(x1(i1-1)-x1(i1+1))) 
             ax1=(x2(i2)-x1(i1-1))*(x2(i2)-x1(i1+1))/((x1(i1)-x1(i1-1))*(x1(i1)-x1(i1+1))) 
             ax2=(x2(i2)-x1(i1-1))*(x2(i2)-x1(i1))/((x1(i1+1)-x1(i1-1))*(x1(i1+1)-x1(i1)))

             az0=(z2(k2)-z1(k1))*(z2(k2)-z1(k1+1))/((z1(k1-1)-z1(k1))*(z1(k1-1)-z1(k1+1))) 
             az1=(z2(k2)-z1(k1-1))*(z2(k2)-z1(k1+1))/((z1(k1)-z1(k1-1))*(z1(k1)-z1(k1+1))) 
             az2=(z2(k2)-z1(k1-1))*(z2(k2)-z1(k1))/((z1(k1+1)-z1(k1-1))*(z1(k1+1)-z1(k1)))

             ux2(i2,j2,k2)=((ay0*ux1(i1-1,j1-1,k1-1)+ay1*ux1(i1-1,j1,k1-1)+ay2*ux1(i1-1,j1+1,k1-1))*ax0+&
			  (ay0*ux1(i1,j1-1,k1-1)+ay1*ux1(i1,j1,k1-1)+ay2*ux1(i1,j1+1,k1-1))*ax1+&
                          (ay0*ux1(i1-1,j1-1,k1-1)+ay1*ux1(i1-1,j1,k1-1)+ay2*ux1(i1-1,j1+1,k1-1))*ax2)*az0+&
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ((ay0*ux1(i1-1,j1-1,k1)+ay1*ux1(i1-1,j1,k1)+ay2*ux1(i1-1,j1+1,k1))*ax0+&
			  (ay0*ux1(i1,j1-1,k1)+ay1*ux1(i1,j1,k1)+ay2*ux1(i1,j1+1,k1))*ax1+&
                          (ay0*ux1(i1-1,j1-1,k1)+ay1*ux1(i1-1,j1,k1)+ay2*ux1(i1-1,j1+1,k1))*ax2)*az1+&
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ((ay0*ux1(i1-1,j1-1,k1+1)+ay1*ux1(i1-1,j1,k1+1)+ay2*ux1(i1-1,j1+1,k1+1))*ax0+&
			  (ay0*ux1(i1,j1-1,k1+1)+ay1*ux1(i1,j1,k1+1)+ay2*ux1(i1,j1+1,k1+1))*ax1+&
                          (ay0*ux1(i1-1,j1-1,k1+1)+ay1*ux1(i1-1,j1,k1+1)+ay2*ux1(i1-1,j1+1,k1+1))*ax2)*az2
             !ux2(i2,j2,k2)=ux1(i1-1,j1-1,k1-1)
             k1_start=k1
             j1_start=j1
             i1_start=i1
             goto 10
         endif

         elseif (i1==nx1 .and. k1/=nz1) then !#######################################################################
         if (y1(j1-1)<=y2(j2) .and. y1(j1)>=y2(j2) .and.&
             x1(i1-1)<=x2(i2) .and. x1(i1)>=x2(i2) .and.&
             z1(k1-1)<=z2(k2) .and. z1(k1)>=z2(k2)) then

             ay0=(y2(j2)-y1(j1))*(y2(j2)-y1(j1+1))/((y1(j1-1)-y1(j1))*(y1(j1-1)-y1(j1+1))) 
             ay1=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1+1))/((y1(j1)-y1(j1-1))*(y1(j1)-y1(j1+1))) 
             ay2=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1))/((y1(j1+1)-y1(j1-1))*(y1(j1+1)-y1(j1)))

             ax0=(x2(i2)-x1(i1-1))*(x2(i2)-x1(i1 ))/((x1(i1-2)-x1(i1-1))*(x1(i1-2)-x1(i1 ))) 
             ax1=(x2(i2)-x1(i1-2))*(x2(i2)-x1(i1 ))/((x1(i1-1)-x1(i1-2))*(x1(i1-1)-x1(i1 ))) 
             ax2=(x2(i2)-x1(i1-2))*(x2(i2)-x1(i1-1))/((x1(i1 )-x1(i1-2))*(x1(i1 )-x1(i1-1)))

             az0=(z2(k2)-z1(k1))*(z2(k2)-z1(k1+1))/((z1(k1-1)-z1(k1))*(z1(k1-1)-z1(k1+1))) 
             az1=(z2(k2)-z1(k1-1))*(z2(k2)-z1(k1+1))/((z1(k1)-z1(k1-1))*(z1(k1)-z1(k1+1))) 
             az2=(z2(k2)-z1(k1-1))*(z2(k2)-z1(k1))/((z1(k1+1)-z1(k1-1))*(z1(k1+1)-z1(k1)))

             ux2(i2,j2,k2)=((ay0*ux1(i1-2,j1-1,k1-1)+ay1*ux1(i1-2,j1,k1-1)+ay2*ux1(i1-2,j1+1,k1-1))*ax0+&
			  (ay0*ux1(i1-1,j1-1,k1-1)+ay1*ux1(i1-1,j1,k1-1)+ay2*ux1(i1-1,j1+1,k1-1))*ax1+&
                          (ay0*ux1(i1-2,j1-1,k1-1)+ay1*ux1(i1-2,j1,k1-1)+ay2*ux1(i1-2,j1+1,k1-1))*ax2)*az0+&
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ((ay0*ux1(i1-2,j1-1,k1)+ay1*ux1(i1-2,j1,k1)+ay2*ux1(i1-2,j1+1,k1))*ax0+&
			  (ay0*ux1(i1-1,j1-1,k1)+ay1*ux1(i1-1,j1,k1)+ay2*ux1(i1-1,j1+1,k1))*ax1+&
                          (ay0*ux1(i1-2,j1-1,k1)+ay1*ux1(i1-2,j1,k1)+ay2*ux1(i1-2,j1+1,k1))*ax2)*az1+&
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ((ay0*ux1(i1-2,j1-1,k1+1)+ay1*ux1(i1-2,j1,k1+1)+ay2*ux1(i1-2,j1+1,k1+1))*ax0+&
			  (ay0*ux1(i1-1,j1-1,k1+1)+ay1*ux1(i1-1,j1,k1+1)+ay2*ux1(i1-1,j1+1,k1+1))*ax1+&
                          (ay0*ux1(i1-2,j1-1,k1+1)+ay1*ux1(i1-2,j1,k1+1)+ay2*ux1(i1-2,j1+1,k1+1))*ax2)*az2
             !ux2(i2,j2,k2)=ux1(i1,j1,k1)
	     if (isnan(ax0) .or. isnan(ax1) .or. isnan(ax2) .or.isnan(ay0) .or. isnan(ay1) &
                  .or. isnan(ay2) .or. isnan(az0) .or. isnan(az1) .or. isnan(az2)) then
	         print *, '-----------------------------------------'
	         print *, i1,j1,k1
                 print *, ax0 ,ax1 , ax2
	         print *, ay0 ,ay1 , ay2
	         print *, az0 ,az1 , az2
	         print *, 'XXXXXXXX',(z1(k1+1)-z1(k1-1)),(z1(k1+1)-z1(k1))
             endif

             k1_start=k1
             j1_start=j1
             i1_start=i1
             goto 10
         endif
         elseif (j1==ny1 .and. k1/=nz1) then !#######################################################################
         if (y1(j1-1)<=y2(j2) .and. y1(j1)>=y2(j2) .and.&
             x1(i1-1)<=x2(i2) .and. x1(i1)>=x2(i2) .and.&
             z1(k1-1)<=z2(k2) .and. z1(k1)>=z2(k2)) then

             ay0=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1-1+1))/((y1(j1-1-1)-y1(j1-1))*(y1(j1-1-1)-y1(j1-1+1))) 
             ay1=(y2(j2)-y1(j1-1-1))*(y2(j2)-y1(j1-1+1))/((y1(j1-1)-y1(j1-1-1))*(y1(j1-1)-y1(j1-1+1))) 
             ay2=(y2(j2)-y1(j1-1-1))*(y2(j2)-y1(j1-1))/((y1(j1-1+1)-y1(j1-1-1))*(y1(j1-1+1)-y1(j1-1)))

             ax0=(x2(i2)-x1(i1))*(x2(i2)-x1(i1+1))/((x1(i1-1)-x1(i1))*(x1(i1-1)-x1(i1+1))) 
             ax1=(x2(i2)-x1(i1-1))*(x2(i2)-x1(i1+1))/((x1(i1)-x1(i1-1))*(x1(i1)-x1(i1+1))) 
             ax2=(x2(i2)-x1(i1-1))*(x2(i2)-x1(i1))/((x1(i1+1)-x1(i1-1))*(x1(i1+1)-x1(i1)))

             az0=(z2(k2)-z1(k1))*(z2(k2)-z1(k1+1))/((z1(k1-1)-z1(k1))*(z1(k1-1)-z1(k1+1))) 
             az1=(z2(k2)-z1(k1-1))*(z2(k2)-z1(k1+1))/((z1(k1)-z1(k1-1))*(z1(k1)-z1(k1+1))) 
             az2=(z2(k2)-z1(k1-1))*(z2(k2)-z1(k1))/((z1(k1+1)-z1(k1-1))*(z1(k1+1)-z1(k1)))

             ux2(i2,j2,k2)=((ay0*ux1(i1-1,j1-1-1,k1-1)+ay1*ux1(i1-1,j1-1,k1-1)+ay2*ux1(i1-1,j1-1+1,k1-1))*ax0+&
			  (ay0*ux1(i1,j1-1-1,k1-1)+ay1*ux1(i1,j1-1,k1-1)+ay2*ux1(i1,j1-1+1,k1-1))*ax1+&
                          (ay0*ux1(i1-1,j1-1-1,k1-1)+ay1*ux1(i1-1,j1-1,k1-1)+ay2*ux1(i1-1,j1-1+1,k1-1))*ax2)*az0+&
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ((ay0*ux1(i1-1,j1-1-1,k1)+ay1*ux1(i1-1,j1-1,k1)+ay2*ux1(i1-1,j1-1+1,k1))*ax0+&
			  (ay0*ux1(i1,j1-1-1,k1)+ay1*ux1(i1,j1-1,k1)+ay2*ux1(i1,j1-1+1,k1))*ax1+&
                          (ay0*ux1(i1-1,j1-1-1,k1)+ay1*ux1(i1-1,j1-1,k1)+ay2*ux1(i1-1,j1-1+1,k1))*ax2)*az1+&
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ((ay0*ux1(i1-1,j1-1-1,k1+1)+ay1*ux1(i1-1,j1-1,k1+1)+ay2*ux1(i1-1,j1-1+1,k1+1))*ax0+&
			  (ay0*ux1(i1,j1-1-1,k1+1)+ay1*ux1(i1,j1-1,k1+1)+ay2*ux1(i1,j1-1+1,k1+1))*ax1+&
                          (ay0*ux1(i1-1,j1-1-1,k1+1)+ay1*ux1(i1-1,j1-1,k1+1)+ay2*ux1(i1-1,j1-1+1,k1+1))*ax2)*az2
             !ux2(i2,j2,k2)=ux1(i1,j1,k1)
	     if (isnan(ax0) .or. isnan(ax1) .or. isnan(ax2) .or.isnan(ay0) .or. isnan(ay1) &
                  .or. isnan(ay2) .or. isnan(az0) .or. isnan(az1) .or. isnan(az2)) then
	     print *, '-----------------------------------------'
	     print *, i1,j1,k1
             print *, ax0 ,ax1 , ax2
	     print *, ay0 ,ay1 , ay2
	     print *, az0 ,az1 , az2
	     print *, '>>>>>YYYYYY>>>>',(z1(k1+1)-z1(k1-1)),(z1(k1+1)-z1(k1))

             endif
             !if (ux1(I1,J1,K1)/=ux2(I2,J2,K2)) then 
                 !print *, ux1(I1,J1,K1),ux2(I2,J2,K2),I1,J1,K1
             !endif
             k1_start=k1
             j1_start=j1
             i1_start=i1
             goto 10
         endif
         elseif (k1==nz1) then !#######################################################################
         if (y1(j1-1)<=y2(j2) .and. y1(j1)>=y2(j2) .and.&
             x1(i1-1)<=x2(i2) .and. x1(i1)>=x2(i2) .and.&
             z1(k1-1)<=z2(k2) .and. z1(k1)>=z2(k2)) then

             ay0=(y2(j2)-y1(j1))*(y2(j2)-y1(j1+1))/((y1(j1-1)-y1(j1))*(y1(j1-1)-y1(j1+1))) 
             ay1=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1+1))/((y1(j1)-y1(j1-1))*(y1(j1)-y1(j1+1))) 
             ay2=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1))/((y1(j1+1)-y1(j1-1))*(y1(j1+1)-y1(j1)))

             ax0=(x2(i2)-x1(i1))*(x2(i2)-x1(i1+1))/((x1(i1-1)-x1(i1))*(x1(i1-1)-x1(i1+1))) 
             ax1=(x2(i2)-x1(i1-1))*(x2(i2)-x1(i1+1))/((x1(i1)-x1(i1-1))*(x1(i1)-x1(i1+1))) 
             ax2=(x2(i2)-x1(i1-1))*(x2(i2)-x1(i1))/((x1(i1+1)-x1(i1-1))*(x1(i1+1)-x1(i1)))

             az0=(z2(k2)-z1(k1-1))*(z2(k2)-z1(k1-1+1))/((z1(k1-1-1)-z1(k1-1))*(z1(k1-1-1)-z1(k1-1+1))) 
             az1=(z2(k2)-z1(k1-1-1))*(z2(k2)-z1(k1-1+1))/((z1(k1-1)-z1(k1-1-1))*(z1(k1-1)-z1(k1-1+1))) 
             az2=(z2(k2)-z1(k1-1-1))*(z2(k2)-z1(k1-1))/((z1(k1-1+1)-z1(k1-1-1))*(z1(k1-1+1)-z1(k1-1)))

             ux2(i2,j2,k2)=((ay0*ux1(i1-1,j1-1,k1-1-1)+ay1*ux1(i1-1,j1,k1-1-1)+ay2*ux1(i1-1,j1+1,k1-1-1))*ax0+&
			  (ay0*ux1(i1,j1-1,k1-1-1)+ay1*ux1(i1,j1,k1-1-1)+ay2*ux1(i1,j1+1,k1-1-1))*ax1+&
                          (ay0*ux1(i1-1,j1-1,k1-1-1)+ay1*ux1(i1-1,j1,k1-1-1)+ay2*ux1(i1-1,j1+1,k1-1-1))*ax2)*az0+&
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ((ay0*ux1(i1-1,j1-1,k1-1)+ay1*ux1(i1-1,j1,k1-1)+ay2*ux1(i1-1,j1+1,k1-1))*ax0+&
			  (ay0*ux1(i1,j1-1,k1-1)+ay1*ux1(i1,j1,k1-1)+ay2*ux1(i1,j1+1,k1-1))*ax1+&
                          (ay0*ux1(i1-1,j1-1,k1-1)+ay1*ux1(i1-1,j1,k1-1)+ay2*ux1(i1-1,j1+1,k1-1))*ax2)*az1+&
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ((ay0*ux1(i1-1,j1-1,k1-1+1)+ay1*ux1(i1-1,j1,k1-1+1)+ay2*ux1(i1-1,j1+1,k1-1+1))*ax0+&
			  (ay0*ux1(i1,j1-1,k1-1+1)+ay1*ux1(i1,j1,k1-1+1)+ay2*ux1(i1,j1+1,k1-1+1))*ax1+&
                          (ay0*ux1(i1-1,j1-1,k1-1+1)+ay1*ux1(i1-1,j1,k1-1+1)+ay2*ux1(i1-1,j1+1,k1-1+1))*ax2)*az2
             !ux2(i2,j2,k2)=ux1(i1,j1,k1)
             k1_start=k1
             j1_start=j1
             i1_start=i1
             goto 10
         endif
         endif
    enddo
    enddo
    enddo
10 enddo 
enddo
enddo

call cpu_time(t2)
print *, 'time = ',t2-t1
!ux2(nx2,:,:)=ux1(nx1,:,:)
!ux2(:,ny2,:)=ux1(:,ny1,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (i_var ==1) then 
    call system('rm '//fileName2)
    open(20, file=fileName2, & 
                 FORM='UNFORMATTED', action="write",access='stream')
else
    open(20,file=fileName2,status='old',action='write',&
         access='stream',form='unformatted',position='append')
endif

print *, 'Writting the variable number : ',i_var
  DO K=1,nz2
     DO J=1,ny2-dny2
        DO I=1,nx2-dnx2
	    if (isnan(ux2(I,J,K))) then 
		print*, 'NAN NAN NAN NAN NAN NAN NAN NAN NAN NAN NAN NAN '
		print*, i,k,k
		stop
	    endif
	    write (20) ux2(I,J,K)
            !if (abs((ux1(I,J,K)-ux2(I,J,K))/ux1(i,j,k))>100) then 
            !if (isnan(ux2(I,J,K))) then 
                !print *, ux1(I,J,K),ux2(I,J,K),I,J,K
            !endif
        ENDDO
     ENDDO
  ENDDO
 close(20)
ENDDO !! THIS THE OUTER LOOP TO SWETCH BETWEEN VARIABLS
 
 close(10)
end program mapping
!

!*******************************************************************
subroutine stretching(yp,yly,ny,beta,istret,ncly,fileName)
! This subroutine is taken form INCOMAPACT3D code.
! https://www.incompact3d.com/
! It determine the grid points locations when stretching is used.
 
!*******************************************************************
! 

implicit none

real(8) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
real(8) :: alpha,beta,pi,yly
integer :: j,istret,ncly,ny
real(8),dimension(:), allocatable :: yeta,yetai,ypi
real(8),intent(out) :: yp(*)
character(len=*) :: fileName
allocate(yeta(ny))
allocate(yetai(ny))
allocate(ypi(ny))
pi=acos(-1.)
yinf=-yly/2.
den=2.*beta*yinf
xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
alpha=abs(xnum/den)
xcx=1./beta/alpha

if (alpha.ne.0.) then
   if (istret.eq.1) yp(1)=0.
   if (istret.eq.2) yp(1)=0.
   if (istret.eq.1) yeta(1)=0.
   if (istret.eq.2) yeta(1)=-1./2.
   if (istret.eq.3) yp(1)=0.
   if (istret.eq.3) yeta(1)=-1./2.
   do j=2,ny!-1
      if (istret==1) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yeta(j)=((j-1.)*(1./2./ny)-0.5)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=((j-1.)*(1./2./(ny-1.))-0.5)
      endif
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yeta(j))+1.
      if ((ncly.ne.0).and.(j==ny).and.(istret==1)) then
         xnum1=0.
      else
         xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
      endif
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf
         if (yeta(j).eq.0.5) yp(j)=0.-yinf
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst+yly
         if (yeta(j).eq.0.5) yp(j)=0.+yly
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yeta(j).lt.0.5) yp(j)=(xnum1-cst+yly)*2.
         if (yeta(j).eq.0.5) yp(j)=(0.+yly)*2.
         if (yeta(j).gt.0.5) yp(j)=(xnum1+cst+yly)*2.
      endif
   enddo

endif
!if (nrank==0) then
!do j=1,ny
!print *,j,yp(j),yeta(j)
!enddo
!endif
!stop
if (alpha.eq.0.) then
   yp(1)=-1.e10
   do j=2,ny
      yeta(j)=(j-1.)*(1./ny)
      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
   enddo
endif
if (alpha.ne.0.) then
   do j=1,ny
      if (istret==1) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./2./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./2./(ny-1.))-0.5
      endif
      
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yetai(j))+1.
      xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst-yinf
         if (yetai(j).eq.0.5) ypi(j)=0.-yinf
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst+yly
         if (yetai(j).eq.0.5) ypi(j)=0.+yly
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yetai(j).lt.0.5) ypi(j)=(xnum1-cst+yly)*2.
         if (yetai(j).eq.0.5) ypi(j)=(0.+yly)*2.
         if (yetai(j).gt.0.5) ypi(j)=(xnum1+cst+yly)*2.
      endif
   enddo
endif
if (alpha.eq.0.) then
   ypi(1)=-1.e10
   do j=2,ny
      yetai(j)=(j-1.)*(1./ny)
      ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
   enddo
endif


open(10,file=fileName, form='formatted')
do j=1,ny
write(10,*)yp(j)
enddo
 close(10)



end subroutine stretching
