SUBROUTINE sub_defmatD
  Use m_mat
  Use m_gen
  Use m_mesh
  Use m_condbord
  implicit none
  INTEGER :: I,INeigh,Node(3),J,K,L,L1,noeud,JJ,testcla,Q
  REAL*8 :: Jfinv,N(2,2),V(3,2)
  REAL*8,allocatable::DerivPhiIDerivPhiJ(:,:)
  REAL*8,allocatable::DerivPhiI(:,:)
  REAL*8,allocatable :: DerivPhi1D(:,:),DerivPhi0D(:,:,:)
  REAL*8 :: norm_point(2),h1,res
  INTEGER,allocatable :: PhiEdge(:,:)
  
  allocate(PhiEdge(3,1+Order))
  allocate(DerivPhi1D(1+Order,NGL1D))
  allocate(DerivPhi0D(2,1+Order,1))

  SELECT CASE(Order)
  CASE(1)

     CALL DerivPhi1DOrder1(PtGL1D,NGL1D,DerivPhi1D)
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
     CALL DerivPhi1DOrder2(PtGL1D,NGL1D,DerivPhi1D)     
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
     CALL DerivPhi1DOrder3(PtGL1D,NGL1D,DerivPhi1D)     
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


allocate(DerivPhiIDerivPhiJ(1+Order,1+Order))
DO I=1,1+Order
   DO J=1,1+Order
      DerivPhiIDerivPhiJ(I,J)=sum(DerivPhi1D(I,:)*DerivPhi1D(J,:)*wGL1D)
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
allocate(Dcla(Nbcla,Nphi,Nphi))
Dcla=0.D0
Nbcla=0
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
           ! TriCla(NbCla)=I
         END If
           DO J=1,1+Order
               DO K=1,1+Order
                  res=DerivPhiIDerivPhiJ(J,K)
                  L=PhiEdge(JJ,K)
                  Q=PhiEdge(JJ,J)
                  Dcla(Nbcla,Q,L)=Dcla(Nbcla,Q,L)+res*sqrt(sum(V(JJ,:)**2))
		end do
	  end do
	end if
   END DO
END DO

allocate(DerivPhiI(2,1+Order))
DO L=1,2
   DO I=1,1+Order
            DerivPhiI(L,I)= sum(DerivPhi0D(L,I,:)*1)  
   END DO
END DO

DO I=1,Nbcla
   Dcla(I,1:Nphi,1:Nphi)=Matmul(Minv,Dcla(I,1:Nphi,1:Nphi))/DFVec(TriCla(I))
ENDDO
END SUBROUTINE sub_defmatD
