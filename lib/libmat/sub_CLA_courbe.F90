SUBROUTINE sub_CLA_courbe(R1)
  Use m_mat
  Use m_gen
!  Use m_cla
  Use m_condbord
  Use m_mesh

  implicit none

  INTEGER :: I,INeigh,Node(3),J,K,L,L1,L2,Q,JJ,NbnoCla,J1,J2
  INTEGER,allocatable :: PhiEdge(:,:)
  REAL*8,allocatable::PhiIPhiJ(:,:)
  REAL*8 :: V(3,2),TestCla
  REAL*8,allocatable :: Phi(:,:),Ainter(:,:,:)
  REAL*8,allocatable :: GradPhi1D(:,:),GradPhi2D(:,:,:)
  REAL*8 :: res,R1,theta1,theta2,d,pi,pt_courbe_x(Order+1),pt_courbe_y(Order+1),dtheta
  pi=2.d0*dasin(1.d0)
  allocate(Phi(1+Order,NGL1D))
  allocate(GradPhi1D(Order+1,NGL1D))
  allocate(PhiEdge(3,1+Order))

  SELECT CASE(Order)
  CASE(1)

     CALL Phi1DOrder1(PtGL1D,NGL1D,Phi)
     CALL Gradphi1Dorder1(PtGL1D,NGL1D,GradPhi1D)
     !Functions on Edge 23
     PhiEdge(1,1)=2
     PhiEdge(1,2)=3
     !Functions on Edge 31
     PhiEdge(2,1)=3
     PhiEdge(2,2)=1
     !Functions on Edge 12
     PhiEdge(3,1)=1
     PhiEdge(3,2)=2

  CASE(2)
     CALL Phi1DOrder2(PtGL1D,NGL1D,Phi)
     CALL Gradphi1Dorder2(PtGL1D,NGL1D,GradPhi1D)
!Functions on Edge 23
PhiEdge(1,1)=2
PhiEdge(1,2)=5
PhiEdge(1,3)=3
!Functions on Edge 31
PhiEdge(2,1)=3
PhiEdge(2,2)=6
PhiEdge(2,3)=1
!Functions on Edge 12
PhiEdge(3,1)=1
PhiEdge(3,2)=4
PhiEdge(3,3)=2

  CASE(3)
     CALL Phi1DOrder3(PtGL1D,NGL1D,Phi)
     CALL Gradphi1Dorder3(PtGL1D,NGL1D,GradPhi1D)
!Functions on Edge 23
PhiEdge(1,1)=2
PhiEdge(1,2)=6
PhiEdge(1,3)=7
PhiEdge(1,4)=3
!Functions on Edge 31
PhiEdge(2,1)=3
PhiEdge(2,2)=8
PhiEdge(2,3)=9
PhiEdge(2,4)=1
!Functions on Edge 12
PhiEdge(3,1)=1
PhiEdge(3,2)=4
PhiEdge(3,3)=5
PhiEdge(3,4)=2
  END SELECT

!!Int(PhiIPhiJ) (1D)
allocate(PhiIPhiJ(1+Order,1+Order))
DO I=1,1+Order
   DO J=1,1+Order
      PhiIPhiJ(I,J)=sum(Phi(I,:)*Phi(J,:)*wGL1D)
   END DO
END DO
NBcla=0
DO I=1,Ntri
   DO J=1,3
      IF (Neigh(I,J)==-3) THEN
         NbCLA=NbCLA+1
exit
      END IF
   ENDDO
ENDDO
allocate(TriCla(NbCLA))
!allocate(Clatri(Ntri))
write(6,*) Ntri
write(6,*) NbCLA
write(6,*) Ntri-NbCla
allocate(BCla(NbCla,Nphi,Nphi))
allocate(BCla2(NbCla,Nphi,Nphi))
BCla=0.D0
BCla2=0.D0
NbCla=0
NbnoCla=0
DO I=1,NTri
   Node = Tri(I,:)
   !!Computation of vectors V12,V23 and V31
   V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
   V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
   V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
   V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
   V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
   V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)
   TestCla=0
   DO JJ=1,3
!!! Searching Neighbor
      INeigh=Neigh(I,JJ)
!!! If there is a Neighbor
      IF(Ineigh.eq.-3) then
         SELECT CASE(JJ)
         CASE(1)
            J1=2
            J2=3
         CASE(2)
            J1=3
            J2=1
         CASE(3)
            J1=1
            J2=2
         END SELECT
pt_courbe_x(1)=Coor(Node(J1),1)
pt_courbe_x(order+1)=Coor(Node(J2),1)
pt_courbe_y(1)=Coor(Node(J1),2)
pt_courbe_y(order+1)=Coor(Node(J2),2)
call calc_theta(coor(Node(J1),1), coor(Node(J1),2),theta1)
call calc_theta(coor(Node(J2),1), coor(Node(J2),2),theta2)
dtheta=theta2-theta1
if (dtheta.gt.pi) then
dtheta=dtheta-2*pi
elseif (dtheta.lt.-pi) then
dtheta=dtheta+2*pi
end if
   DO K=2,order
      pt_courbe_x(k)= R1*dcos(theta1+(k-1)*dtheta/real(order,8))
      pt_courbe_y(k)= R1*dsin(theta1+(k-1)*dtheta/real(order,8))
   ENDDO

If (TestCla==0) then
            TestCla=1
            NbCla=NbCla+1
            TriCla(NbCla)=I
           ! Clatri(I)=Nbcla
         END If
         !! [u][v]
         DO J=1,Order+1
            DO K=1,Order+1
                res=0D0
                DO L=1,NGL1D
                   d=sum(pt_courbe_x*gradphi1D(:,L))**2
                   d=d+sum(pt_courbe_y*gradphi1D(:,L))**2
                   res=res+Phi(J,L)*Phi(K,L)*wGL1D(L)*sqrt(d)
                ENDDO
!write(6,*) res,sqrt(sum(V(JJ,:)**2))
               L=PhiEdge(JJ,K)
               Q=PhiEdge(JJ,J)
               Bcla(NbCla,Q,L)=Bcla(NbCla,Q,L)+res
            END DO
         END DO
      END IF
   END DO
   IF (TestCla==0) then
      NbnoCla=NbnoCla+1
   END IF

END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!    
!!!!!!!!! Computation of Minv*C
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
If (helmholtz.eq.0) then
DO I=1,Nbcla
   Bcla(I,1:Nphi,1:Nphi)=Matmul(Minv,Bcla(I,1:Nphi,1:Nphi))/DFVec(TriCla(I))*mu(TriCla(I))
ENDDO
end If
!write(6,*) NbnoCla+Nbcla,Ntri

!	open(21,File="FILM/matrice_CLA.dat")
!		DO I=1,NbCla
!		write(21,*) I,Bcla(I,:,:)
!		ENDDO
!	close(21)

 

END SUBROUTINE sub_CLA_courbe
