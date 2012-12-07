!!http://people.sc.fsu.edu/~jburkardt/f_src/specfun/specfun.f90
subroutine elastoacoustic_exact_asymp_Vero(k,kp,ks,R0,R1,Nsomme,PhiEdge,epsilon)
  Use m_condbord    
  Use m_mat
  Use m_mesh
  Use m_gen
  implicit none
  !Pour Mumps
  !Fin pour mumps
  real*8 :: somme,k,kp,ks,R0,R1,coor_x,coor_y,r,theta,V12(2),V23(2),V31(2),gamma13,gamma23,theta0,condi,condi1,res
  integer :: i,l, Nsomme, NCALC,int_j,KK,type_cla,j,N, Iloc,lp,ls,nsomme2,nsomme3
  real*8,dimension(:),allocatable :: J0,Y0,J1,Y1,J0p,Y0p,J1p,Y1p,J0s,Y0s,J1s,Y1s,Jr,Yr,Jrp,Jrs
  complex*16,dimension(:),allocatable :: Hankel0,Hankel1,coeffA,coeffB,coeffC,coeffD, rhssystcoeff,coeff
  complex*16,dimension(:,:),allocatable :: matsystcoeff,matsystcoeffinv

  complex*16 :: tmp,c,d,dd,alphab,betab,ur,utheta,pr,alphac,betac,deltac
  INTEGER :: PhiEdge(3,Order+1),Node(3)

  INTEGER ::Lwork,INFO,IPIV(4)
  real *8 :: WORK,epsilon,heps
  heps=1+epsilon/R0
	nsomme2=nsomme 
nsomme2=nsomme
nsomme3=nsomme
write(6,*) 'heho',nsomme,nsomme2,k
  N=2
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
  allocate(matsystcoeff(N,N))
  allocate(rhssystcoeff(N))

  NCALC=0

  call RJBESL(k*R0,0.D0,Nsomme+2,J0,NCALC)	
  call RYBESL(k*R0,0.D0,Nsomme+2,Y0,NCALC)


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

  condi=1!omega**2!1D0/(omega**2*rho(1)/k)
	condi1=1!1e-13
!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Calcul des coefficients 
!!!!!!!!!!!!!!!!!!!!!!!!!!
 gamma13=2.6789385347d0
     gamma23=1.3541179394d0
alphac=-(0D0,1D0)*k-1/(2D0*R0)+(6/R0)**(1/3D0)*gamma23/gamma13*(1/2D0-(0,1)*sqrt(3D0)/2D0)*k**(2/3D0) 
betac=-(6/R0)**(1/3D0)*gamma23/gamma13*(1/2D0-(0,1)*sqrt(3D0)/2D0)*k**(2/3D0)*(-(0D0,1D0)*k+1/(2D0*R0))+k**2
betac=betac/alphac 
  DO i=1,Nsomme+1
deltac=(betac-(i-1)**2/R0**2/alphac)
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Matrice du système
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$     matsystcoeff(1,1) =(k*((i-1)/(k*R0)*Hankel0(i)-Hankel0(i+1)))&
!!$          &- (k*((i-1)/(k*R0)*conjg(Hankel0(i))-conjg(Hankel0(i+1))))
!!$     matsystcoeff(1,2) = -omega**2*rho(1)/R0*((i-1)*J0p(i)-kp*R0*J0p(i+1))
!!$     matsystcoeff(1,3) = -omega**2*rho(1)/R0*(i-1)*J0s(i)
!!$!     matsystcoeff(2,1) = 1D0
!!$!     matsystcoeff(2,2) = 1D0
!!$!     matsystcoeff(2,3) = 0.d0
!!$!     matsystcoeff(2,4) = 0.d0
!!$     matsystcoeff(2,1) = (Hankel0(i))-(conjg(Hankel0(i)))
!!$     matsystcoeff(2,2) = 2.d0*Cij(1,3,3)/(R0**2)*(((i-1)**2+(i-1)-1/2.d0*ks**2*R0**2)*J0p(i)-& 
!!$          & kp*R0*(2.d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1)))
!!$     matsystcoeff(2,3) = 2.d0*Cij(1,3,3)/(R0**2)*((i-1)*(-i*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1))))
!!$     matsystcoeff(3,1) = 0.d0
!!$!     matsystcoeff(4,2) = 0.d0
!!$     matsystcoeff(3,2) = 2.d0*Cij(1,3,3)/(R0**2)*(-(i-1)*(-i*J0p(i)+&
!!$&kp*R0*(2d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1))))
!!$     matsystcoeff(3,3) = 2.d0*Cij(1,3,3)/(R0**2)*(-((i-1)**2+(i-1)-&
!!$&1/2.d0*ks**2*R0**2)*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1)))
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
!!$!     rhssystcoeff(2) = 0.d0
!!$     rhssystcoeff(3) = 0.d0
!!$
!!$     matsystcoeff(2:3,2:3)=condi1*matsystcoeff(2:3,2:3)
!!$     matsystcoeff(1,1)=matsystcoeff(1,1)/condi1/condi
!!$     matsystcoeff(1,2:3)=matsystcoeff(1,2:3)/condi
!!$!     matsystcoeff(:,4)=condi*matsystcoeff(:,4)
!!$rhssystcoeff(1:2)=rhssystcoeff(1:2)/condi/condi1
!!$
!!$     call ZGESV( 3, 1, matsystcoeff, 3, IPIV, rhssystcoeff, 3, INFO )
!!$     coeffA(i)=rhssystcoeff(1)
!!$     coeffB(i)=-rhssystcoeff(1)
!!$     coeffC(i)=condi1*rhssystcoeff(2)
!!$     coeffD(i)=condi1*rhssystcoeff(3)

!     matsystcoeff(1,1) =(k*((i-1)/(k*R0)*Hankel0(i)-Hankel0(i+1)))
!     matsystcoeff(1,2) = (k*((i-1)/(k*R0)*conjg(Hankel0(i))-conjg(Hankel0(i+1))))
!     matsystcoeff(1,3) = -omega**2*rho(1)/R0*((i-1)*J0p(i)-kp*R0*J0p(i+1))
!     matsystcoeff(1,4) = -omega**2*rho(1)/R0*(i-1)*J0s(i)
!     matsystcoeff(2,1) = Hankel1(i)
!     matsystcoeff(2,2) = conjg(Hankel1(i))
!     matsystcoeff(2,3) = 0.d0
!     matsystcoeff(2,4) = 0.d0
!     matsystcoeff(3,1) = (Hankel0(i))
!     matsystcoeff(3,2) = (conjg(Hankel0(i)))
     matsystcoeff(1,1) = 2.d0*Cij(1,3,3)/(R0**2)*(((i-1)**2+(i-1)-1/2.d0*ks**2*R0**2)*J0p(i)-& 
          & kp*R0*(2.d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1)))
     matsystcoeff(1,1) = matsystcoeff(1,1)*deltac 
     matsystcoeff(1,1) = matsystcoeff(1,1)+omega**2*rho(1)/R0&
          &*((i-1)*J0p(i)-kp*R0*J0p(i+1))
     matsystcoeff(1,2) = 2.d0*Cij(1,3,3)/(R0**2)*((i-1)*(-i*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1))))
 matsystcoeff(1,2) =  matsystcoeff(1,2)*deltac
     matsystcoeff(1,2) =  matsystcoeff(1,2)+omega**2*rho(1)/R0*(i-1)*J0s(i)
!     matsystcoeff(4,1) = 0.d0
!     matsystcoeff(4,2) = 0.d0
     matsystcoeff(2,1) = 2.d0*Cij(1,3,3)/(R0**2)*(-(i-1)*(-i*J0p(i)+&
&kp*R0*(2d0*(i-1)/(kp*R0)*J0p(i)-J0p(i+1))))
     matsystcoeff(2,2) = 2.d0*Cij(1,3,3)/(R0**2)*(-((i-1)**2+(i-1)-&
&1/2.d0*ks**2*R0**2)*J0s(i)+ks*R0*(2.d0*(i-1)/(ks*R0)*J0s(i)-J0s(i+1)))

!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!! Second membre du système
!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (i==1) then
        rhssystcoeff(1) =-cmplx(0,1)**(i-1)*J0(i)
        rhssystcoeff(1) =rhssystcoeff(1)*deltac
        
        rhssystcoeff(1) =rhssystcoeff(1)+k*(cmplx(0,1))**(i-1)*((i-1)/(k*R0)*J0(i)-J0(i+1))
     else
        rhssystcoeff(1) =-2D0*cmplx(0,1)**(i-1)*J0(i)
        rhssystcoeff(1) =rhssystcoeff(1)*deltac
        
        rhssystcoeff(1) =rhssystcoeff(1)+2D0*k*(cmplx(0,1))**(i-1)*((i-1)/(k*R0)*J0(i)-J0(i+1))
     endif

     rhssystcoeff(2) = 0.d0
!!     rhssystcoeff(4) = 0.d0

!!$     matsystcoeff(3:4,3:4)=condi1*matsystcoeff(3:4,3:4)
!!$     matsystcoeff(1:2,1:2)=matsystcoeff(1:2,1:2)/condi1/condi
!!$     matsystcoeff(1:2,3:4)=matsystcoeff(1:2,3:4)/condi
!     matsystcoeff(:,4)=condi*matsystcoeff(:,4)
!!rhssystcoeff(1:2)=rhssystcoeff(1:2)/condi/condi1

     call ZGESV( N, 1, matsystcoeff, N, IPIV, rhssystcoeff, N, INFO )
!!$     coeffA(i)=rhssystcoeff(1)
!!$     coeffB(i)=rhssystcoeff(2)
!!$     coeffC(i)=condi1*rhssystcoeff(3)
!!$     coeffD(i)=condi1*rhssystcoeff(4)
 coeffC(i)=rhssystcoeff(1)
     coeffD(i)=rhssystcoeff(2)

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
!!$if(i.eq.1) then
!!$write(6,*) 'i', i
!!$write(6,*) 'coeffA', coeffA(i)
!!$write(6,*) 'coeffB', coeffB(i)
!!$write(6,*) 'coeffC', coeffC(i)
!!$write(6,*) 'coeffD', coeffD(i)
!!$else
!!$write(6,*) 'i', i
!!$write(6,*) 'coeffA', coeffA(i),abs((coeffA(i)-coeffA(i-1))/coeffA(i-1))
!!$write(6,*) 'coeffB', coeffB(i),abs((coeffB(i)-coeffB(i-1))/coeffB(i-1))
!!$write(6,*) 'coeffC', coeffC(i),abs((coeffC(i)-coeffC(i-1))/coeffC(i-1))
!!$write(6,*) 'coeffD', coeffD(i),abs((coeffD(i)-coeffD(i-1))/coeffD(i-1))
!!$endif
  ENDDO
write(6,*) 'toto'
  P_analytic=0.d0
  Ux_analytic=0.d0
  Uy_analytic=0.d0
!!$write(6,*) CoeffC
!!$write(6,*) CoeffD
!!$stop


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calcul de la solution exacte p 
!!!!!!!!!!!!!!!!!!!!!!!!!!
        res=10.D0
        l=1

!!$  Do J=1,Nflu+Nflusol
!!$ ! Do J=Nflu,Nflu+Nflusol
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calcul de la solution exacte p aux ddl d'ordre 1
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$     DO int_j=1,3
!!$        coor_x =dsqrt(2D0)/2D0*0.01D0! Coor(Tri(J,int_j),1)
!!$        coor_y =dsqrt(2D0)/2D0*0.01D0! Coor(Tri(J,int_j),2)
!!$        coor_x =Coor(Tri(J,int_j),1)
!!$        coor_y =Coor(Tri(J,int_j),2)
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$           P_analytic((J-1)*Nphi+int_j)=P_analytic((J-1)*Nphi+int_j)+pr
!!$
!!$     enddo
!!$
!!$     SELECT CASE(Order)
!!$
!!$     CASE(2)
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calcul de l'onde p en ajoutant les ddl de l'ordre 2
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$        Node = Tri(J,:)
!!$
!!$        !!Computation of vectors V12,V23 and V31
!!$        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
!!$        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
!!$        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
!!$        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
!!$        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
!!$        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)
!!$
!!$        coor_x=Coor(Tri(J,1),1)+V12(1)/Order
!!$        coor_y=Coor(Tri(J,1),2)+V12(2)/Order
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+PhiEdge(3,2))=P_analytic((J-1)*Nphi+PhiEdge(3,2))+pr
!!$
!!$        coor_x=Coor(Tri(J,2),1)+V23(1)/Order
!!$        coor_y=Coor(Tri(J,2),2)+V23(2)/Order
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+PhiEdge(1,2))=P_analytic((J-1)*Nphi+PhiEdge(1,2))+pr
!!$
!!$        coor_x=Coor(Tri(J,3),1)+V31(1)/Order
!!$        coor_y=Coor(Tri(J,3),2)+V31(2)/Order
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+PhiEdge(2,2))=P_analytic((J-1)*Nphi+PhiEdge(2,2))+pr
!!$
!!$
!!$     CASE(3)
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calcul de l'onde p en ajoutant les ddl de l'ordre 3 supplémentaires à l'ordre 1
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$        Node = Tri(J,:)
!!$        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
!!$        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
!!$        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
!!$        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
!!$        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
!!$        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)
!!$
!!$        DO KK=1,Order-1
!!$           coor_x=Coor(Tri(J,1),1)+KK*V12(1)/Order
!!$           coor_y=Coor(Tri(J,1),2)+KK*V12(2)/Order
!!$           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$           call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$           call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$           call calc_theta(coor_x,coor_y,theta)
!!$
!!$           pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$                &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$              P_analytic((J-1)*Nphi+PhiEdge(3,KK+1))=P_analytic((J-1)*Nphi+PhiEdge(3,KK+1))+pr
!!$
!!$        ENDDO
!!$
!!$
!!$        DO KK=1,Order-1
!!$           coor_x=Coor(Tri(J,2),1)+KK*V23(1)/Order
!!$           coor_y=Coor(Tri(J,2),2)+KK*V23(2)/Order
!!$           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$           call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$           call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$           call calc_theta(coor_x,coor_y,theta)
!!$
!!$           pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$                &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$              P_analytic((J-1)*Nphi+PhiEdge(1,KK+1))=P_analytic((J-1)*Nphi+PhiEdge(1,KK+1))+pr
!!$
!!$        ENDDO
!!$
!!$        DO KK=1,Order-1
!!$           coor_x=Coor(Tri(J,3),1)+KK*V31(1)/Order
!!$           coor_y=Coor(Tri(J,3),2)+KK*V31(2)/Order
!!$           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$           call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$           call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$           call calc_theta(coor_x,coor_y,theta)
!!$
!!$           pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$                &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$              P_analytic((J-1)*Nphi+PhiEdge(2,KK+1))=P_analytic((J-1)*Nphi+PhiEdge(2,KK+1))+pr
!!$
!!$
!!$        ENDDO
!!$
!!$        coor_x=Coor(Tri(J,1),1)+V12(1)/3.-V31(1)/3.
!!$        coor_y=Coor(Tri(J,1),2)+V12(2)/3.-V31(2)/3.
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+Nphi)=P_analytic((J-1)*Nphi+Nphi)+pr
!!$
!!$
!!$        ENDSELECT
!!$
!!$     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calcul de u exact
!!!!!!!!!!!!!!!!!!!!!!!!!!

  Do J=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu+Nsol
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calcul de u exact aux ddl d'ordre 1
!!!!!!!!!!!!!!!!!!!!!!!!!!

     I=((J-1)-Nflu-Nflusol)*Nphi
     DO int_j=1,3

        coor_x = Coor(Tri(J,int_j),1)
        coor_y = Coor(Tri(J,int_j),2)
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I+int_j) =Ux_analytic(I+int_j)+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I+int_j) =Uy_analytic(I+int_j)+ur*dsin(theta)+utheta*dcos(theta)
        enddo
     SELECT CASE(Order)
     CASE(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calcul de u aux ddl d'ordre 2
!!!!!!!!!!!!!!!!!!!!!!!!!!


        I=((J-1)-Nflu-Nflusol)*Nphi  

        Node = Tri(J,:)

        !!Computation of vectors V12,V23 and V31
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)

        coor_x=Coor(Tri(J,1),1)+V12(1)/Order
        coor_y=Coor(Tri(J,1),2)+V12(2)/Order
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I+PhiEdge(3,2)) =Ux_analytic(I+PhiEdge(3,2))+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I+PhiEdge(3,2)) =Uy_analytic(I+PhiEdge(3,2))+ur*dsin(theta)+utheta*dcos(theta)

        coor_x=Coor(Tri(J,2),1)+V23(1)/Order
        coor_y=Coor(Tri(J,2),2)+V23(2)/Order
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I+PhiEdge(1,2)) =Ux_analytic(I+PhiEdge(1,2))+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I+PhiEdge(1,2)) =Uy_analytic(I+PhiEdge(1,2))+ur*dsin(theta)+utheta*dcos(theta)

        coor_x=Coor(Tri(J,3),1)+V31(1)/Order
        coor_y=Coor(Tri(J,3),2)+V31(2)/Order
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I+PhiEdge(2,2)) =Ux_analytic(I+PhiEdge(2,2))+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I+PhiEdge(2,2)) =Uy_analytic(I+PhiEdge(2,2))+ur*dsin(theta)+utheta*dcos(theta)


     CASE(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calcul de u aux ddl d'ordre 3
!!!!!!!!!!!!!!!!!!!!!!!!!!

        Node = Tri(J,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)

        DO KK=1,Order-1
           I=((J-1)-Nflu-Nflusol)*Nphi+PhiEdge(3,KK+1)
           coor_x=Coor(Tri(J,1),1)+KK*V12(1)/Order
           coor_y=Coor(Tri(J,1),2)+KK*V12(2)/Order
           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
           call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
           call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
           call calc_theta(coor_x,coor_y,theta)

           ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
           utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

              Ux_analytic(I) =Ux_analytic(I)+ur*dcos(theta)-utheta*dsin(theta)

              Uy_analytic(I) =Uy_analytic(I)+ur*dsin(theta)+utheta*dcos(theta)


        ENDDO

        DO KK=1,Order-1
           I=((J-1)-Nflu-Nflusol)*Nphi+PhiEdge(1,KK+1)
           coor_x=Coor(Tri(J,2),1)+KK*V23(1)/Order
           coor_y=Coor(Tri(J,2),2)+KK*V23(2)/Order
           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
           call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
           call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
           call calc_theta(coor_x,coor_y,theta)

           ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
           utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

              Ux_analytic(I) =Ux_analytic(I)+ur*dcos(theta)-utheta*dsin(theta)

              Uy_analytic(I) =Uy_analytic(I)+ur*dsin(theta)+utheta*dcos(theta)

        ENDDO

        DO KK=1,Order-1
           I=((J-1)-Nflu-Nflusol)*Nphi+PhiEdge(2,KK+1)
           coor_x=Coor(Tri(J,3),1)+KK*V31(1)/Order
           coor_y=Coor(Tri(J,3),2)+KK*V31(2)/Order
           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
           call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
           call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
           call calc_theta(coor_x,coor_y,theta)

           ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
           utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

              Ux_analytic(I) =Ux_analytic(I)+ur*dcos(theta)-utheta*dsin(theta)

              Uy_analytic(I) =Uy_analytic(I)+ur*dsin(theta)+utheta*dcos(theta)


        ENDDO

        I=((J-1)-Nflu-Nflusol)*Nphi+Nphi
        coor_x=Coor(Tri(J,1),1)+V12(1)/3.-V31(1)/3.
        coor_y=Coor(Tri(J,1),2)+V12(2)/3.-V31(2)/3.
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I) =Ux_analytic(I)+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I) =Uy_analytic(I)+ur*dsin(theta)+utheta*dcos(theta)

     ENDSELECT

  enddo


res=10D0
!toto
  do while (res.gt.1D-15.and.(l.lt.Nsomme3+1)) 
res=0D0
           l=l+1
!!$  Do J=1,Nflu+Nflusol
!!$ ! Do J=Nflu,Nflu+Nflusol
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calcul de la solution exacte p aux ddl d'ordre 1
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$     DO int_j=1,3
!!$        coor_x =dsqrt(2D0)/2D0*0.01D0! Coor(Tri(J,int_j),1)
!!$        coor_y =dsqrt(2D0)/2D0*0.01D0! Coor(Tri(J,int_j),2)
!!$        coor_x =Coor(Tri(J,int_j),1)
!!$        coor_y =Coor(Tri(J,int_j),2)
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+int_j)=P_analytic((J-1)*Nphi+int_j)+pr
!!$           res=max(res,abs(pr/P_analytic((J-1)*Nphi+int_j)))
!!$!           write(6,*)  res,l,coor_x,coor_y
!!$        enddo
!!$
!!$
!!$     SELECT CASE(Order)
!!$
!!$     CASE(2)
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calcul de l'onde p en ajoutant les ddl de l'ordre 2
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$        Node = Tri(J,:)
!!$
!!$        !!Computation of vectors V12,V23 and V31
!!$        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
!!$        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
!!$        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
!!$        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
!!$        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
!!$        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)
!!$
!!$        coor_x=Coor(Tri(J,1),1)+V12(1)/Order
!!$        coor_y=Coor(Tri(J,1),2)+V12(2)/Order
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+PhiEdge(3,2))=P_analytic((J-1)*Nphi+PhiEdge(3,2))+pr
!!$           res=max(res,abs(pr/P_analytic((J-1)*Nphi+PhiEdge(3,2))))
!!$     
!!$
!!$        coor_x=Coor(Tri(J,2),1)+V23(1)/Order
!!$        coor_y=Coor(Tri(J,2),2)+V23(2)/Order
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+PhiEdge(1,2))=P_analytic((J-1)*Nphi+PhiEdge(1,2))+pr
!!$           res=max(res,abs(pr/P_analytic((J-1)*Nphi+PhiEdge(1,2))))
!!$
!!$        coor_x=Coor(Tri(J,3),1)+V31(1)/Order
!!$        coor_y=Coor(Tri(J,3),2)+V31(2)/Order
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+PhiEdge(2,2))=P_analytic((J-1)*Nphi+PhiEdge(2,2))+pr
!!$           res=max(res,abs(pr/P_analytic((J-1)*Nphi+PhiEdge(2,2))))
!!$
!!$
!!$     CASE(3)
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calcul de l'onde p en ajoutant les ddl de l'ordre 3 supplémentaires à l'ordre 1
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$        Node = Tri(J,:)
!!$        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
!!$        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
!!$        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
!!$        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
!!$        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
!!$        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)
!!$
!!$        DO KK=1,Order-1
!!$           coor_x=Coor(Tri(J,1),1)+KK*V12(1)/Order
!!$           coor_y=Coor(Tri(J,1),2)+KK*V12(2)/Order
!!$           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$           call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$           call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$           call calc_theta(coor_x,coor_y,theta)
!!$
!!$           pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$                &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$              P_analytic((J-1)*Nphi+PhiEdge(3,KK+1))=P_analytic((J-1)*Nphi+PhiEdge(3,KK+1))+pr
!!$              res=max(res,abs(pr/P_analytic((J-1)*Nphi+PhiEdge(3,KK+1))))
!!$
!!$        ENDDO
!!$
!!$
!!$        DO KK=1,Order-1
!!$           coor_x=Coor(Tri(J,2),1)+KK*V23(1)/Order
!!$           coor_y=Coor(Tri(J,2),2)+KK*V23(2)/Order
!!$           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$           call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$           call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$           call calc_theta(coor_x,coor_y,theta)
!!$
!!$           pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$                &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$              P_analytic((J-1)*Nphi+PhiEdge(1,KK+1))=P_analytic((J-1)*Nphi+PhiEdge(1,KK+1))+pr
!!$              res=max(res,abs(pr/P_analytic((J-1)*Nphi+PhiEdge(1,KK+1))))
!!$
!!$        ENDDO
!!$
!!$        DO KK=1,Order-1
!!$           coor_x=Coor(Tri(J,3),1)+KK*V31(1)/Order
!!$           coor_y=Coor(Tri(J,3),2)+KK*V31(2)/Order
!!$           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$           call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$           call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$           call calc_theta(coor_x,coor_y,theta)
!!$
!!$           pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$                &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$              P_analytic((J-1)*Nphi+PhiEdge(2,KK+1))=P_analytic((J-1)*Nphi+PhiEdge(2,KK+1))+pr
!!$              res=max(res,abs(pr/P_analytic((J-1)*Nphi+PhiEdge(2,KK+1))))
!!$
!!$
!!$        ENDDO
!!$
!!$        coor_x=Coor(Tri(J,1),1)+V12(1)/3.-V31(1)/3.
!!$        coor_y=Coor(Tri(J,1),2)+V12(2)/3.-V31(2)/3.
!!$        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
!!$        call RJBESL(k*r,0.D0,Nsomme3+1,Jr,NCALC)	
!!$        call RYBESL(k*r,0.D0,Nsomme3+1,Yr,NCALC)
!!$        call calc_theta(coor_x,coor_y,theta)
!!$
!!$        pr=(coeffA(l)*cmplx(Jr(l),Yr(l))+&
!!$             &coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
!!$
!!$           P_analytic((J-1)*Nphi+Nphi)=P_analytic((J-1)*Nphi+Nphi)+pr
!!$           res=max(res,abs(pr/P_analytic((J-1)*Nphi+Nphi)))
!!write(6,*) l,res
!!$
!!$
!!$     ENDSELECT
!!$
!!$  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calcul de u exact
!!!!!!!!!!!!!!!!!!!!!!!!!!

  Do J=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu+Nsol
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calcul de u exact aux ddl d'ordre 1
!!!!!!!!!!!!!!!!!!!!!!!!!!

     I=((J-1)-Nflu-Nflusol)*Nphi
     DO int_j=1,3

        coor_x = Coor(Tri(J,int_j),1)
        coor_y = Coor(Tri(J,int_j),2)
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I+int_j) =Ux_analytic(I+int_j)+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I+int_j) =Uy_analytic(I+int_j)+ur*dsin(theta)+utheta*dcos(theta)
           res=max(res,2*dsqrt(abs(ur)**2+abs(utheta)**2)/&
                &dsqrt(abs(Ux_analytic(I+int_j))**2+&
                &abs(Uy_analytic(I+int_j))**2))
     enddo
     SELECT CASE(Order)
     CASE(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calcul de u aux ddl d'ordre 2
!!!!!!!!!!!!!!!!!!!!!!!!!!


        I=((J-1)-Nflu-Nflusol)*Nphi  

        Node = Tri(J,:)

        !!Computation of vectors V12,V23 and V31
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)

        coor_x=Coor(Tri(J,1),1)+V12(1)/Order
        coor_y=Coor(Tri(J,1),2)+V12(2)/Order
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I+PhiEdge(3,2)) =Ux_analytic(I+PhiEdge(3,2))+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I+PhiEdge(3,2)) =Uy_analytic(I+PhiEdge(3,2))+ur*dsin(theta)+utheta*dcos(theta)
           res=max(res,2*dsqrt(abs(ur)**2+abs(utheta)**2)/&
                &dsqrt(abs(Ux_analytic(I+PhiEdge(3,2)))**2+&
                &abs(Uy_analytic(I+PhiEdge(3,2)))**2))

        coor_x=Coor(Tri(J,2),1)+V23(1)/Order
        coor_y=Coor(Tri(J,2),2)+V23(2)/Order
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I+PhiEdge(1,2)) =Ux_analytic(I+PhiEdge(1,2))+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I+PhiEdge(1,2)) =Uy_analytic(I+PhiEdge(1,2))+ur*dsin(theta)+utheta*dcos(theta)
           res=max(res,2*dsqrt(abs(ur)**2+abs(utheta)**2)/&
                &dsqrt(abs(Ux_analytic(I+PhiEdge(1,2)))**2+&
                &abs(Uy_analytic(I+PhiEdge(1,2)))**2))

        coor_x=Coor(Tri(J,3),1)+V31(1)/Order
        coor_y=Coor(Tri(J,3),2)+V31(2)/Order
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I+PhiEdge(2,2)) =Ux_analytic(I+PhiEdge(2,2))+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I+PhiEdge(2,2)) =Uy_analytic(I+PhiEdge(2,2))+ur*dsin(theta)+utheta*dcos(theta)
           res=max(res,2*dsqrt(abs(ur)**2+abs(utheta)**2)/&
                &dsqrt(abs(Ux_analytic(I+PhiEdge(2,2)))**2+&
                &abs(Uy_analytic(I+PhiEdge(2,2)))**2))


     CASE(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calcul de u aux ddl d'ordre 3
!!!!!!!!!!!!!!!!!!!!!!!!!!

        Node = Tri(J,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)

        DO KK=1,Order-1
           I=((J-1)-Nflu-Nflusol)*Nphi+PhiEdge(3,KK+1)
           coor_x=Coor(Tri(J,1),1)+KK*V12(1)/Order
           coor_y=Coor(Tri(J,1),2)+KK*V12(2)/Order
           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
           call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
           call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
           call calc_theta(coor_x,coor_y,theta)

           ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
           utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

              Ux_analytic(I) =Ux_analytic(I)+ur*dcos(theta)-utheta*dsin(theta)

              Uy_analytic(I) =Uy_analytic(I)+ur*dsin(theta)+utheta*dcos(theta)
              res=max(res,2*dsqrt(abs(ur)**2+abs(utheta)**2)/&
                   &dsqrt(abs(Ux_analytic(I))**2+&
                   &abs(Uy_analytic(I))**2))


        ENDDO

        DO KK=1,Order-1
           I=((J-1)-Nflu-Nflusol)*Nphi+PhiEdge(1,KK+1)
           coor_x=Coor(Tri(J,2),1)+KK*V23(1)/Order
           coor_y=Coor(Tri(J,2),2)+KK*V23(2)/Order
           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
           call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
           call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
           call calc_theta(coor_x,coor_y,theta)

           ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
           utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

              Ux_analytic(I) =Ux_analytic(I)+ur*dcos(theta)-utheta*dsin(theta)

              Uy_analytic(I) =Uy_analytic(I)+ur*dsin(theta)+utheta*dcos(theta)
              res=max(res,2*dsqrt(abs(ur)**2+abs(utheta)**2)/&
                   &dsqrt(abs(Ux_analytic(I))**2+&
                   &abs(Uy_analytic(I))**2))

        ENDDO

        DO KK=1,Order-1
           I=((J-1)-Nflu-Nflusol)*Nphi+PhiEdge(2,KK+1)
           coor_x=Coor(Tri(J,3),1)+KK*V31(1)/Order
           coor_y=Coor(Tri(J,3),2)+KK*V31(2)/Order
           r=dsqrt(coor_x*coor_x+coor_y*coor_y)
           call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
           call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
           call calc_theta(coor_x,coor_y,theta)

           ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
           utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

              Ux_analytic(I) =Ux_analytic(I)+ur*dcos(theta)-utheta*dsin(theta)

              Uy_analytic(I) =Uy_analytic(I)+ur*dsin(theta)+utheta*dcos(theta)
              res=max(res,2*dsqrt(abs(ur)**2+abs(utheta)**2)/&
                   &dsqrt(abs(Ux_analytic(I))**2+&
                   &abs(Uy_analytic(I))**2))


        ENDDO

        I=((J-1)-Nflu-Nflusol)*Nphi+Nphi
        coor_x=Coor(Tri(J,1),1)+V12(1)/3.-V31(1)/3.
        coor_y=Coor(Tri(J,1),2)+V12(2)/3.-V31(2)/3.
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(kp*r,0.D0,Nsomme2+2,Jrp,NCALC)	
        call RJBESL(ks*r,0.D0,Nsomme2+2,Jrs,NCALC)	
        call calc_theta(coor_x,coor_y,theta)

        ur=(coeffC(l)*((l-1)/(kp*r)*Jrp(l)-Jrp(l+1))*kp+coeffD(l)*(l-1)/r*Jrs(l))*dcos((l-1)*theta)
        utheta=(-coeffC(l)*(l-1)/r*Jrp(l)-coeffD(l)*((l-1)/(ks*r)*Jrs(l)-Jrs(l+1))*ks)*dsin((l-1)*theta)

           Ux_analytic(I) =Ux_analytic(I)+ur*dcos(theta)-utheta*dsin(theta)

           Uy_analytic(I) =Uy_analytic(I)+ur*dsin(theta)+utheta*dcos(theta)
           res=max(res,2*dsqrt(abs(ur)**2+abs(utheta)**2)/&
                &dsqrt(abs(Ux_analytic(I))**2+&
                &abs(Uy_analytic(I))**2))

     ENDSELECT
enddo
write(6,*) 'Modes', res,l,Nsomme
end do

  deallocate(J0)
  deallocate(Y0)
  deallocate(J0p)
  deallocate(J0s)
  deallocate(Jr)
  deallocate(Yr)
  deallocate(Jrp)
  deallocate(Jrs)
  deallocate(matsystcoeff)
  deallocate(rhssystcoeff)
  deallocate(coeffA)
  deallocate(coeffB)
  deallocate(coeffC)
  deallocate(coeffD)
end subroutine elastoacoustic_exact_asymp_Vero
