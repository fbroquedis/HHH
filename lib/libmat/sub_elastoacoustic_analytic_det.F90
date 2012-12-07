!!http://people.sc.fsu.edu/~jburkardt/f_src/specfun/specfun.f90

subroutine elastoacoustic_analytic_det(k,kp,ks,R0,R1,Nsomme,PhiEdge)
  Use m_condbord    
  Use m_mat
  Use m_mesh
  Use m_gen
  implicit none
  !Pour Mumps
!  INCLUDE 'mpif.h'
!  INCLUDE 'zmumps_struc.h'
  !Fin pour mumps
  real*8 :: somme,k,kp,ks,R0,R1,coor_x,coor_y,r,theta,V12(2),V23(2),V31(2),gamma13,gamma23,theta0,condi,condi1,res
  integer :: i,l, Nsomme, NCALC,int_j,KK,type_cla,j,N, Iloc,lp,ls,nsomme2
  real*8,dimension(:),allocatable :: J0,Y0,J1,Y1,J0p,Y0p,J1p,Y1p,J0s,Y0s,J1s,Y1s,Jr,Yr,Jrp,Jrs
  complex*16,dimension(:),allocatable :: Hankel0,Hankel1,coeffA,coeffB,coeffC,coeffD, rhssystcoeff,coeff
  complex*16,dimension(:,:),allocatable :: matsystcoeff,matsystcoeffbis

  complex*16 :: tmp,c,d,dd,alphab,betab,ur,utheta,pr,d1,d2,c1,c2,a1,a2,b1,b2
  INTEGER :: PhiEdge(3,Order+1),Node(3)

  INTEGER ::Lwork,INFO,IPIV(4)
  real *8 :: WORK
  complex*16, dimension(:),allocatable ::det1,det2,det3
  !Pour Mumps	
!  TYPE (ZMUMPS_STRUC) id
  INTEGER IERR,nproc
!  integer, dimension(MPI_STATUS_SIZE) :: status
	nsomme2=nsomme
!write(6,*) 'heho',nsomme,nsomme2
  N=4
!!$  id%SYM = 0
!!$  id%PAR = 1
!!$  id%JOB = -1
!!$  id%icntl(4)=4
!!$  id%comm=MPI_COMM_WORLD
!!$  CALL ZMUMPS(id)
!!$  !fin pour mumps
!!$  IF ( id%MYID .eq. 0 ) THEN
!!$     id%N = N
!!$     id%NZ = N**2
!!$     ALLOCATE( id%IRN ( id%NZ ) )
!!$     ALLOCATE( id%JCN ( id%NZ ) )
!!$     ALLOCATE( id%A( id%NZ ) )
!!$     ALLOCATE( id%RHS ( id%N ) )
!!$  END IF

  allocate(J0(Nsomme+2))
  allocate(Y0(Nsomme+2))
  allocate(J0p(Nsomme+2))
  allocate(J0s(Nsomme+2))
  allocate(Hankel0(Nsomme+2))
  allocate(J1(Nsomme+2))
  allocate(Y1(Nsomme+2))
  allocate(Hankel1(Nsomme+2))

  allocate(Jr(Nsomme+2))
  allocate(Yr(Nsomme+2))
  allocate(Jrp(Nsomme+2))
  allocate(Jrs(Nsomme+2))

  allocate(coeffA(Nsomme+1))
  allocate(coeffB(Nsomme+1))
  allocate(coeffC(Nsomme+1))
  allocate(coeffD(Nsomme+1))
  allocate(det1(Nsomme+1))
  allocate(det2(Nsomme+1))
  allocate(det3(Nsomme+1))
  allocate(matsystcoeff(N,N))
  allocate(matsystcoeffbis(N,N))
  allocate(rhssystcoeff(N))

  NCALC=0

  call RJBESL(k*R0,0.D0,Nsomme+2,J0,NCALC)	
  call RYBESL(k*R0,0.D0,Nsomme+2,Y0,NCALC)

!!$  write(6,*) 'kR0', k*R0
!!$  write(6,*) 'Bessel', J0(1), Y0(1), J0(2), Y0(2), J0(5), Y0(5)

  Hankel0 = cmplx(J0,Y0)

  call RJBESL(kp*R0,0.D0,Nsomme+2,J0p,NCALC)	
  call RJBESL(ks*R0,0.D0,Nsomme+2,J0s,NCALC)	


  call RJBESL(k*R1,0.D0,Nsomme+2,J1,NCALC)	
  call RYBESL(k*R1,0.D0,Nsomme+2,Y1,NCALC)

  Hankel1 = cmplx(J1,Y1)

  coeffA=0D0
  coeffB=0D0
  coeffC=0D0
  coeffD=0D0
  condi1=.5D0/(omega**2*rho(1))
	condi=1D0!1/1e-13
!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Calcul des coefficients 
!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i=1,Nsomme+1

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Matrice du système
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
     if (i==1) then
        rhssystcoeff(1) = (-k*(cmplx(0,1))**(i-1)*((i-1)/(k*R0)*J0(i)-J0(i+1)))
        rhssystcoeff(3) = (-(cmplx(0,1))**(i-1)*J0(i))
     else
        rhssystcoeff(1) = (-2*k*(cmplx(0,1))**(i-1)*((i-1)/(k*R0)*J0(i)-J0(i+1)))
        rhssystcoeff(3) = (-2*(cmplx(0,1))**(i-1)*J0(i))
     endif

     rhssystcoeff(2) = 0.d0
     rhssystcoeff(4) = 0.d0

matsystcoeff=0D0
 matsystcoeff(1,1) =0D0
     matsystcoeff(1,1) =(k*((i-1)/(k*R0)*Hankel0(i)-Hankel0(i+1)))
     matsystcoeff(1,2) = (k*((i-1)/(k*R0)*conjg(Hankel0(i))-conjg(Hankel0(i+1))))
     matsystcoeff(1,3) = -omega**2*rho(1)/R0*((i-1)*J0p(i)-kp*R0*J0p(i+1))
     matsystcoeff(1,4) = -omega**2*rho(1)/R0*(i-1)*J0s(i)
     matsystcoeff(2,1) = (k*((i-1)/(k*R1)*Hankel1(i)-Hankel1(i+1)) - cmplx(0,1)*k*Hankel1(i))
     matsystcoeff(2,2) =(k*((i-1)/(k*R1)&
          &*conjg(Hankel1(i))-conjg(Hankel1(i+1))) - cmplx(0,1)*k*conjg(Hankel1(i)))
     matsystcoeff(2,3) = 0.d0
     matsystcoeff(2,4) = 0.d0
     matsystcoeff(3,1) = (Hankel0(i))
     matsystcoeff(3,2) = (conjg(Hankel0(i)))
     matsystcoeff(3,3) = 2.d0*Cij(1,3,3)/(R0**2)*(((i-1)**2+(i-1)-1/2.d0*ks**2*R0**2)*J0p(i)-& 
          & kp*R0*(2.d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1)))
!!     write(6,*) 'heho',k,ks,kp,Cij(1,3,3),J0p(i)
     matsystcoeff(3,4) = 2.d0*Cij(1,3,3)/(R0**2)*((i-1)*(-i*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1))))
!    matsystcoeff(3,3) = (((i-1)**2+(i-1)-1/2.d0*ks**2*R0**2)*J0p(i)-& 
!          & kp*R0*(2.d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1)))
!!     write(6,*) 'heho',k,ks,kp,Cij(1,3,3),J0p(i)
!     matsystcoeff(3,4) = ((i-1)*(-i*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1))))
     matsystcoeff(4,1) = 0.d0
     matsystcoeff(4,2) = 0.d0
     matsystcoeff(4,3) = 2.d0*Cij(1,3,3)/(R0**2)*(-(i-1)*(-i*J0p(i)+&
&kp*R0*(2d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1))))
     matsystcoeff(4,4) = 2.d0*Cij(1,3,3)/(R0**2)*(-((i-1)**2+(i-1)-&
&1/2.d0*ks**2*R0**2)*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1)))
!!$     matsystcoeff(4,3) = (-(i-1)*(-i*J0p(i)+&
!!$&kp*R0*(2d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1))))
!!$     matsystcoeff(4,4) = (-((i-1)**2+(i-1)-&
!!$&1/2.d0*ks**2*R0**2)*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1)))

matsystcoeff(:,3:4)=condi1*matsystcoeff(:,3:4)
!!matsystcoeff(3:4,:)=matsystcoeff(3:4,:)*100D0
!!$if(i.eq.1) then
!!$write(6,*) '1',abs(matsystcoeff(1,:))
!!$write(6,*) '2',abs(matsystcoeff(2,:))
!!$write(6,*) '3',abs(matsystcoeff(3,:))
!!$write(6,*) '4',abs(matsystcoeff(4,:))
!!$end if
!!$ matsystcoeff(3:4,3:4)=condi1*matsystcoeff(3:4,3:4)
!!$     matsystcoeff(1:2,1:2)=matsystcoeff(1:2,1:2)/condi1/condi
!!$  
!!$ matsystcoeff(1:2,3:4)=matsystcoeff(1:2,3:4)/condi
!!$     matsystcoeff(:,4)=condi*matsystcoeff(:,4)
!!$
!!$rhssystcoeff(1:2)=rhssystcoeff(1:2)/condi/condi1
!!$
matsystcoeffbis=matsystcoeff
matsystcoeff(1,:)=matsystcoeff(1,:)/dsqrt(sum(abs(matsystcoeff(1,:))**2))
matsystcoeff(2,:)=matsystcoeff(2,:)/dsqrt(sum(abs(matsystcoeff(2,:))**2))
matsystcoeff(3,:)=matsystcoeff(3,:)/dsqrt(sum(abs(matsystcoeff(3,:))**2))
matsystcoeff(4,:)=matsystcoeff(4,:)/dsqrt(sum(abs(matsystcoeff(4,:))**2))

!matsystcoeffbis=matsystcoeff
    call ZGESV( 4, 1, matsystcoeff, 4, IPIV, rhssystcoeff, 4, INFO )
det1(i)=matsystcoeff(1,1)*matsystcoeff(2,2)*matsystcoeff(3,3)*matsystcoeff(4,4)
matsystcoeff=0.D0
matsystcoeff(1:2,1:2)=matsystcoeffbis(3:4,3:4)
matsystcoeff(1,:)=matsystcoeff(1,:)/dsqrt(sum(abs(matsystcoeff(1,:))**2))
matsystcoeff(2,:)=matsystcoeff(2,:)/dsqrt(sum(abs(matsystcoeff(2,:))**2))
matsystcoeff(3,:)=matsystcoeff(3,:)/dsqrt(sum(abs(matsystcoeff(3,:))**2))
matsystcoeff(4,:)=matsystcoeff(4,:)/dsqrt(sum(abs(matsystcoeff(4,:))**2))
!!$if(I.eq.1) then
!!$write(6,*)'3',matsystcoeff(1,1:2)
!!$write(6,*)'4',matsystcoeff(2,1:2)
!!$end if
    call ZGESV( 2, 1, matsystcoeff(1:2,1:2), 2, IPIV(1:2), rhssystcoeff(1:2), 2, INFO )
det3(i)=matsystcoeff(1,1)*matsystcoeff(2,2)
matsystcoeff=0D0
matsystcoeff(2:3,2:3)=matsystcoeffbis(3:4,3:4)
matsystcoeff(1,1)=matsystcoeffbis(1,1)
matsystcoeff(1,2:3)=matsystcoeffbis(1,3:4)
matsystcoeff(2:3,1)=matsystcoeffbis(3:4,1)
matsystcoeff(1,:)=matsystcoeff(1,:)/dsqrt(sum(abs(matsystcoeff(1,:))**2))
matsystcoeff(2,:)=matsystcoeff(2,:)/dsqrt(sum(abs(matsystcoeff(2,:))**2))
matsystcoeff(3,:)=matsystcoeff(3,:)/dsqrt(sum(abs(matsystcoeff(3,:))**2))
matsystcoeff(4,:)=matsystcoeff(4,:)/dsqrt(sum(abs(matsystcoeff(4,:))**2))
    call ZGESV( 3, 1, matsystcoeff(1:3,1:3), 3, IPIV(1:3), rhssystcoeff(1:3), 3, INFO )
det2(i)=matsystcoeff(1,1)*matsystcoeff(2,2)*matsystcoeff(3,3)


!!$     matsystcoeff(4,3) = (-(i-1)*(-i*J0p(i)+&
!!$&kp*R0*(2d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1))))
!!$     matsystcoeff(4,4) = (-((i-1)**2+(i-1)-&
!!$&1/2.d0*ks**2*R0**2)*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1)))
!!$!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!$!!! Second membre du système
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
     if (i==1) then
        rhssystcoeff(1) = (-k*(cmplx(0,1))**(i-1)*((i-1)/(k*R0)*J0(i)-J0(i+1)))
        rhssystcoeff(3) = (-(cmplx(0,1))**(i-1)*J0(i))
     else
        rhssystcoeff(1) = (-2*k*(cmplx(0,1))**(i-1)*((i-1)/(k*R0)*J0(i)-J0(i+1)))
        rhssystcoeff(3) = (-2*(cmplx(0,1))**(i-1)*J0(i))
     endif
!!$
!!$     rhssystcoeff(2) = 0.d0
!!$     rhssystcoeff(4) = 0.d0
!!$
 
!!$     coeffA(i)=rhssystcoeff(1)
!!$     coeffB(i)=rhssystcoeff(2)
!!$     coeffC(i)=condi1*rhssystcoeff(3)
!!$     coeffD(i)=condi1*rhssystcoeff(4)

!!$     d1 = 2.d0*Cij(1,3,3)/(R0**2)*(-(i-1)*(-i*J0p(i)+&
!!$&kp*R0*(2d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1))))
!!$     d2 = 2.d0*Cij(1,3,3)/(R0**2)*(-((i-1)**2+(i-1)-&
!!$&1/2.d0*ks**2*R0**2)*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1)))
!!$     c1 = (k*((i-1)/(k*R1)*Hankel1(i)-Hankel1(i+1)) - cmplx(0,1)*k*Hankel1(i))
!!$     c2 =(k*((i-1)/(k*R1)&
!!$          &*conjg(Hankel1(i))-conjg(Hankel1(i+1))) - cmplx(0,1)*k*conjg(Hankel1(i)))
!!$
!!$     matsystcoeff(1,1) =(k*((i-1)/(k*R0)*Hankel0(i)-Hankel0(i+1)))
!!$     matsystcoeff(1,1) = matsystcoeff(1,1)-(k*((i-1)/(k*R0)*conjg(Hankel0(i))-conjg(Hankel0(i+1))))*c1/c2
!!$     matsystcoeff(1,2) = -omega**2*rho(1)/R0*((i-1)*J0p(i)-kp*R0*J0p(i+1))
!!$     matsystcoeff(1,2) =  matsystcoeff(1,2)+omega**2*rho(1)/R0*(i-1)*J0s(i)*d1/d2
!!$!     matsystcoeff(2,3) = 0.d0
!!$     matsystcoeff(2,1) = (Hankel0(i))
!!$     matsystcoeff(2,1) = matsystcoeff(2,1)-(conjg(Hankel0(i)))*c1/c2
!!$     matsystcoeff(2,2) = 2.d0*Cij(1,3,3)/(R0**2)*(((i-1)**2+(i-1)-1/2.d0*ks**2*R0**2)*J0p(i)-& 
!!$          & kp*R0*(2.d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1)))
!!$     matsystcoeff(2,2) = matsystcoeff(2,2)-2.d0*Cij(1,3,3)/(R0**2)&
!!$          &*((i-1)*(-i*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1))))*d1/d2
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!$!!! Second membre du système
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$     if (i==1) then
!!$        rhssystcoeff(1) = (-k*(cmplx(0,1))**(i-1)*((i-1)/(k*R0)*J0(i)-J0(i+1)))
!!$        rhssystcoeff(2) = (-(cmplx(0,1))**(i-1)*J0(i))
!!$     else
!!$        rhssystcoeff(1) = (-2*k*(cmplx(0,1))**(i-1)*((i-1)/(k*R0)*J0(i)-J0(i+1)))
!!$        rhssystcoeff(2) = (-2*(cmplx(0,1))**(i-1)*J0(i))
!!$     endif
!!$
!!$
!!$     matsystcoeff(2,2)=condi1*matsystcoeff(2,2)
!!$     matsystcoeff(1,1)=matsystcoeff(1,1)/condi1/condi
!!$     matsystcoeff(1,2)=matsystcoeff(1,2)/condi
!!$     matsystcoeff(:,4)=condi*matsystcoeff(:,4)
!!$rhssystcoeff(1)=rhssystcoeff(1)/condi/condi1
!!$     a1= matsystcoeff(1,1)
!!$     a2=matsystcoeff(1,2)
!!$     b1=matsystcoeff(2,1)
!!$     b2=matsystcoeff(2,2)
!!$     coeffA(i)=(b2*rhssystcoeff(1)-a2*rhssystcoeff(2))/(a1*b2-a2*b1)
!!$     coeffB(i)=-coeffA(i)*c1/c2
!!$     coeffC(i)=condi1*(-b1*rhssystcoeff(1)+a1*rhssystcoeff(2))/(a1*b2-a2*b1)
!!$!     coeffC(i)=condi1*rhssystcoeff(2)
!!$     coeffD(i)=-coeffC(i)*d1/d2
!!$     call ZGESV( 2, 1, matsystcoeff(1:2,1:2), 2, IPIV(1:2), rhssystcoeff(1:2), 2, INFO )
!!$     coeffA(i)=rhssystcoeff(1)
!!$     coeffB(i)=rhssystcoeff(1)*c1/c2
!!$     coeffC(i)=condi1*rhssystcoeff(2)
!!$     coeffD(i)=-condi1*rhssystcoeff(2)*d1/d2

!!$!!! Résolution du système avec Mumps

!!$     Iloc=0
!!$     id%RHS(1:N)=rhssystcoeff(1:N)
!!$     DO j=1,N
!!$        DO l=1,N
!!$           Iloc=Iloc+1
!!$           id%IRN(Iloc)=j
!!$           id%JCN(Iloc)=l
!!$           id%A(Iloc)=matsystcoeff(j,l)
!!$        ENDDO
!!$     ENDDO
!!$     id%JOB = 6
!!$     CALL ZMUMPS(id)

!!$    write(6,*)'rhs', rhssystcoeff

!!$     coeffA(i)=id%rhs(1)
!!$     coeffB(i)=id%rhs(2)
!!$     coeffC(i)=id%rhs(3)
!!$     coeffD(i)=id%rhs(4)

!!$     write(6,*) 'mat', matmul(matsystcoeff,rhssystcoeff)
!!$write(6,*) 'i', i
!!$write(6,*) 'coeffA', coeffA(i)
!!$write(6,*) 'coeffB', coeffB(i)
!!$write(6,*) 'coeffC', coeffC(i)
!!$write(6,*) 'coeffD', coeffD(i)
!matsystcoeff(3:4,3)=matsystcoeff(3:4,3)*maxval(abs(matsystcoeff(3:4,4)))/maxval(abs(matsystcoeff(3:4,3)))
!!$write(6,*) '3',matsystcoeff(3,3:4)
!!$write(6,*) '4',matsystcoeff(4,3:4)

!!$det3(i)=(matsystcoeff(3,3)*matsystcoeff(4,4)-matsystcoeff(3,4)*matsystcoeff(4,3))
!!$det2(i)=(matsystcoeff(1,1)-matsystcoeff(1,3))*det3(i)
!!$det1(i)=(matsystcoeff(1,2)-matsystcoeff(2,2))*det2(i)
!!$det3(i)=det3(i)/dsqrt(sum(abs(matsystcoeff(3,3:4))**2))/dsqrt(sum(abs(matsystcoeff(4,3:4))**2))
!!$det2(i)=det2(i)/dsqrt(sum(abs(matsystcoeff(1,1:4))**2))/&
!!$    &dsqrt(sum(abs(matsystcoeff(3,1:4))**2))/dsqrt(sum(abs(matsystcoeff(4,1:4))**2))
!!$det1(i)=det1(i)/dsqrt(sum(abs(matsystcoeff(1,1:4))**2))/dsqrt(sum(abs(matsystcoeff(2,1:4))**2))/&
!!$    &dsqrt(sum(abs(matsystcoeff(3,1:4))**2))/dsqrt(sum(abs(matsystcoeff(4,1:4))**2))
!!$det3(i)=(matsystcoeff(3,3)*matsystcoeff(4,4)-matsystcoeff(3,4)*matsystcoeff(4,3))/&
!!$     &dsqrt(sum(abs(matsystcoeff(3,3:4))**2))/dsqrt(sum(abs(matsystcoeff(4,3:4))**2))
!!$det2(i)=(matsystcoeff(1,1)-matsystcoeff(1,3))*(matsystcoeff(3,3)*matsystcoeff(4,4)-matsystcoeff(3,4)*matsystcoeff(4,3))/&
!!$     &dsqrt(sum(abs(matsystcoeff(1,1:4))**2))/&
!!$     &dsqrt(sum(abs(matsystcoeff(3,1:4))**2))/dsqrt(sum(abs(matsystcoeff(4,1:4))**2))
!!$det1(i)=(matsystcoeff(1,2)-matsystcoeff(2,2))*det2(i)/&
!!$     &dsqrt(sum(abs(matsystcoeff(2,1:4))**2))
!!$det3(i)=(matsystcoeff(3,3)*matsystcoeff(4,4)-matsystcoeff(3,4)*matsystcoeff(4,3))
!!$det2(i)=matsystcoeff(1,1)*det3(i)-matsystcoeff(3,1)*matsystcoeff(1,3)*matsystcoeff(4,4)&
!!$     &+matsystcoeff(3,1)*matsystcoeff(1,4)*matsystcoeff(4,3)
!!$det1(i)=matsystcoeff(1,2)*(matsystcoeff(1,1)*matsystcoeff(3,3)*matsystcoeff(4,4)+&
!!$     &matsystcoeff(3,1)*matsystcoeff(4,3)*matsystcoeff(1,4)-&
!!$     &matsystcoeff(3,1)*matsystcoeff(1,3)*matsystcoeff(4,4)-&
!!$     &matsystcoeff(1,1)*matsystcoeff(3,4)*matsystcoeff(4,3))-matsystcoeff(2,2)*det2(i)
!det3(i)=matsystcoeff(3,3)*matsystcoeff(4,4)
!det2(i)=det3(i)*matsystcoeff(1,1)

  ENDDO
write(32,*) k*r0,abs(det3(1)),abs(det3(2)),abs(det3(3)),abs(det3(4)),abs(det3(5)),abs(det3(6))
write(31,*) k*r0,abs(det2(1)),abs(det2(2)),abs(det2(3)),abs(det2(4)),abs(det2(5)),abs(det2(6))
write(30,*) k*r0,abs(det1(1)),abs(det1(2)),abs(det1(3)),abs(det1(4)),abs(det1(5)),abs(det1(6))

end subroutine elastoacoustic_analytic_det

