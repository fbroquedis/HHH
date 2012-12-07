SUBROUTINE sub_defmatG
  Use m_mat
  Use m_gen
  Use m_condbord
  Use m_mesh

  implicit none

  INTEGER :: I,INeigh,Node(3),J,K,L,L1,L2,Q,JJ,NbnoCla
  INTEGER,allocatable :: PhiEdge(:,:)
  REAL*8,allocatable::PhiIPhiJ(:,:)
  REAL*8 :: V(3,2),TestCla
  REAL*8,allocatable :: Phi(:,:),Ainter(:,:,:)
  REAL*8,allocatable :: GradPhi1D(:,:,:,:),GradPhi2D(:,:,:)
  REAL*8 :: res
  allocate(Phi(1+Order,NGL1D))
  allocate(PhiEdge(3,1+Order))

  SELECT CASE(Order)
  CASE(1)

     CALL Phi1DOrder1(PtGL1D,NGL1D,Phi)
 
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
allocate(G(NbCla,Nphi,Nphi))
G=0.D0
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
         If (TestCla==0) then
            TestCla=1
            NbCla=NbCla+1
         END If
         !! [u][v]
         DO J=1,Order+1
            DO K=1,Order+1
               res=PhiIPhiJ(J,K)
               L=PhiEdge(JJ,K)
               Q=PhiEdge(JJ,J)
               G(NbCla,Q,L)=G(NbCla,Q,L)+res*sqrt(sum(V(JJ,:)**2))
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
if (helmholtz.eq.0)then
DO I=1,Nbcla
   G(I,1:Nphi,1:Nphi)=Matmul(Minv,G(I,1:Nphi,1:Nphi))/DFVec(TriCla(I))
ENDDO
end if
END SUBROUTINE sub_defmatG
