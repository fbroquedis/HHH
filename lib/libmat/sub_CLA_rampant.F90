SUBROUTINE sub_CLA_rampant
  Use m_mat
  Use m_gen
!  Use m_cla
  Use m_condbord
  Use m_mesh

  implicit none

  INTEGER :: I,INeigh,Node(3),J,K,L,L1,L2,Q,JJ,NbnoCla,tri_glob
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


allocate(BCla(Narete_du_bord,Nphi,Nphi))

BCla=0.D0


DO I=1,Narete_du_bord
   tri_glob=tri_bord(I)
   Node = Tri(tri_glob,:)
   !!Computation of vectors V12,V23 and V31
   V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
   V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
   V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
   V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
   V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
   V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)

   DO JJ=1,3
!!! Searching Neighbor
      INeigh=Neigh(tri_glob,JJ)
!!! If there is a Neighbor
      IF(Ineigh.eq.-3) then
         DO J=1,Order+1
            DO K=1,Order+1
               res=PhiIPhiJ(J,K)
               L=PhiEdge(JJ,K)
               Q=PhiEdge(JJ,J)
               Bcla(I,Q,L)=Bcla(I,Q,L)+res*sqrt(sum(V(JJ,:)**2))
            END DO
         END DO
      END IF
   END DO
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!    
!!!!!!!!! Computation of Minv*C
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(6,*) Bcla(1,:,:)
If (helmholtz.eq.0) then
DO I=1,Nbcla
   Bcla(I,1:Nphi,1:Nphi)=Matmul(Minv,Bcla(I,1:Nphi,1:Nphi))/DFVec(TriCla(I))
ENDDO
end If
!write(6,*) NbnoCla+Nbcla,Ntri

 

	open(21,File="FILM/matrice_CLA_ramp.dat")
		DO I=1,Narete_du_bord
		write(21,*) I,Bcla(I,:,:)
		ENDDO
	close(21)

write(6,*) 'BCLA', Bcla(31,1,:) 
write(6,*) 'BCLA', Bcla(31,2,:) 
write(6,*) 'BCLA', Bcla(31,3,:) 
write(6,*) 'BCLA', Bcla(31,4,:) 
 



END SUBROUTINE sub_CLA_rampant
