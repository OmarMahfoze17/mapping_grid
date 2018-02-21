!############################################################################
!
program mapping
!
!############################################################################
  implicit none
  integer :: i,j,k,i1,j1,k1,i2,j2,k2,j1_start,i_var
  integer :: i1_start,k1_start,count,dny1,dny2,dnx1,dnx2
  integer, parameter :: nx1=128 , ny1=129, nz1=256, N_VAR=10
  real(8),dimension(nx1, ny1, nz1) :: ux1,uy1,uz1
  real(8),dimension(ny1) :: y1
  real(8),dimension(nx1) :: x1
  real(8),dimension(nz1) :: z1
  real(8) ::xlx1,yly1,zlz1,dx1,dz1,t1,t2 
  real(8) ::ay0,ay1,ay2,ax0,ax1,ax2,az0,az1,az2,za0
!=========================================================
  integer,parameter:: nx2=128, ny2=129, nz2=128
  real(8),dimension(nx2, ny2, nz2) :: ux2,uy2,uz2
  real(8),dimension(ny2) :: y2
  real(8),dimension(nx2) :: x2
  real(8),dimension(nz2) :: z2
  real(8) ::xlx2,yly2,zlz2,dx2,dz2,dy2 
  integer :: np=1 ! polynomial order
  logical :: exist
  character(len=100) :: fileName1='sauve1.dat'
  character(len=100) :: fileName2='/home/om1014/PhD/INCOMPACT3D/Channel/channel/sauve.dat'

  xlx1=12.56637061
  yly1=2.
  zlz1=4.188790205

  xlx2=12.56637061
  yly2=2.
  zlz2=4.188790205
  dx1=xlx1/(nx1-1)
  dz1=zlz1/(nz1-1)
  DO I=1,nx1
     x1(i)=dx1*(i-1)
  enddo
   open (15,file='yp1.dat',form='formatted',status='unknown')
   do j=1,ny1
      read(15,*) y1(j)
   enddo
   close(15)
  DO K=1,nz1
     z1(k)=dz1*(k-1)
  enddo


  dx2=xlx2/(nx2-1)
  dy2=yly2/(ny2-1)
  dz2=zlz2/(nz2-1)
  DO I=1,nx2
     x2(i)=dx2*(i-1)
  enddo
  DO j=1,ny2
     y2(j)=dy2*(j-1)
  enddo
  DO K=1,nz2
     z2(k)=dz2*(k-1)
  enddo


!! =========== TEST THE DIMENSIONS OF THE BOUNDARY ================
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

print *, 'i_VAR = ',i_var
  DO K=1,nz1
     DO J=1,ny1-dny1
        DO I=1,nx1-dnx1
           READ(10,REC=COUNT) ux1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO




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
             
             k1_start=k1
             j1_start=j1
             i1_start=i1
             goto 10
         endif

         elseif (i1==nx1) then !#######################################################################
         if (y1(j1-1)<=y2(j2) .and. y1(j1)>=y2(j2) .and.&
             x1(i1-1)<=x2(i2) .and. x1(i1)>=x2(i2) .and.&
             z1(k1-1)<=z2(k2) .and. z1(k1)>=z2(k2)) then

             ay0=(y2(j2)-y1(j1))*(y2(j2)-y1(j1+1))/((y1(j1-1)-y1(j1))*(y1(j1-1)-y1(j1+1))) 
             ay1=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1+1))/((y1(j1)-y1(j1-1))*(y1(j1)-y1(j1+1))) 
             ay2=(y2(j2)-y1(j1-1))*(y2(j2)-y1(j1))/((y1(j1+1)-y1(j1-1))*(y1(j1+1)-y1(j1)))

             ax0=(x2(i2)-x1(i1-1))*(x2(i2)-x1(i1))/((x1(i1-2)-x1(i1-1))*(x1(i1-2)-x1(i1))) 
             ax1=(x2(i2)-x1(i1-2))*(x2(i2)-x1(i1))/((x1(i1-1)-x1(i1-2))*(x1(i1-1)-x1(i1))) 
             ax2=(x2(i2)-x1(i1-2))*(x2(i2)-x1(i1-1))/((x1(i1)-x1(i1-2))*(x1(i1)-x1(i1-1)))

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
             
             k1_start=k1
             j1_start=j1
             i1_start=i1
             goto 10
         endif
         elseif (j1==ny1) then !#######################################################################
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

  DO K=1,nz2
     DO J=1,ny2-dny2
        DO I=1,nx2-dnx2
	    write (20) ux2(I,J,K)
        ENDDO
     ENDDO
  ENDDO
 close(20)
ENDDO !! THIS THE OUTER LOOP TO SWETCH BETWEEN VARIABLS
 
 close(10)
end program mapping
!
