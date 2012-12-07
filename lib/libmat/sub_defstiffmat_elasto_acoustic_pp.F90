SUBROUTINE sub_defstiffmat_elasto_acoustic_pp
  Use m_mat
  Use m_gen
  Use m_mesh
  implicit none
  INTEGER :: I,INeigh,Node(3),J,K,L,L1,Q,Arete,Arete_neigh
  INTEGER,allocatable :: PhiEdge(:,:)
  REAL*8 :: Jfinv(2,2),V(3,2),Test(3),cmax,V_voisin(3,2)
  real*8 :: Jfinv_voisin(2,2)
  REAL*8,allocatable::GradPhiIGradPhiJ(:,:,:,:)
  REAL*8,allocatable::GradPhiIPhiJ(:,:,:,:)
  REAL*8,allocatable::PhiIPhiJ(:,:)
  REAL*8,allocatable :: Phi(:,:),Ainter(:,:,:)
  REAL*8,allocatable :: GradPhi1D(:,:,:,:),GradPhi2D(:,:,:)
  REAL*8 :: norm(3,2),h1,h2,hmin
  complex*16::res
  allocate(Phi(1+Order,NGL1D))
  allocate(PhiEdge(3,1+Order))
  allocate(GradPhi1D(3,Nphi,2,NGL1D))
  allocate(GradPhi2D(Nphi,2,NGL2D))
  SELECT CASE(Order)
  CASE(1)

     CALL GradPhiOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     !On edge 1 (23)
     CALL GradPhiOrder1(1-PtGL1D,PtGL1D,NGL1D,GradPhi1D(1,:,:,:))
     !On edge 2 (31)
     CALL GradPhiOrder1(0*PtGL1D,1-PtGL1D,NGL1D,GradPhi1D(2,:,:,:))
     !On edge 3 (12)
     CALL GradPhiOrder1(PtGL1D,0*PtGL1D,NGL1D,GradPhi1D(3,:,:,:))
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
     CALL GradPhiOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
!On edge 1 (23)
     CALL GradPhiOrder2(1.-PtGL1D,PtGL1D,NGL1D,GradPhi1D(1,:,:,:))
!On edge 2 (31)
     CALL GradPhiOrder2(0.*PtGL1D,1.-PtGL1D,NGL1D,GradPhi1D(2,:,:,:))
!On edge 3 (12)
     CALL GradPhiOrder2(PtGL1D,0.*PtGL1D,NGL1D,GradPhi1D(3,:,:,:))
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
     CALL GradPhiOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
!On edge 1 (23)
     CALL GradPhiOrder3(1.-PtGL1D,PtGL1D,NGL1D,GradPhi1D(1,:,:,:))
!On edge 2 (31)
     CALL GradPhiOrder3(0.*PtGL1D,1.-PtGL1D,NGL1D,GradPhi1D(2,:,:,:))
!On edge 1 (12)
     CALL GradPhiOrder3(PtGL1D,0.*PtGL1D,NGL1D,GradPhi1D(3,:,:,:))
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

!!Int(GradphiI.GradphiJ) (2D)
allocate(GradPhiIGradPhiJ(Nphi,Nphi,2,2))
!! GradPhiIGradPhiJ(I,J,K,L) =dPhi_I/dx_K*dPhi_J/dx_L
DO I=1,Nphi
   DO J=1,Nphi
      DO K=1,2
         DO L=1,2
            GradPhiIGradPhiJ(I,J,K,L)=sum(GradPhi2D(I,K,:)*GradPhi2D(J,L,:)*wGL2D)
         END DO
      END DO
   END DO
END DO
allocate(GradPhiIPhiJ(3,Nphi,1+Order,2))
! GradPhiIPhiJ(L,I,J,K)=dphi_I/dx_K*phi_J on face L
DO L=1,3
   DO I=1,Nphi
      DO J=1,1+Order
         DO K=1,2
            GradPhiIPhiJ(L,I,J,K)= sum(GradPhi1D(L,I,K,:)*Phi(J,:)*wGL1D)  
         END DO
      END DO
   END DO
END DO
!!Int(PhiIPhiJ) (1D)
allocate(PhiIPhiJ(1+Order,1+Order))
! GradPhiIPhiJ(I,J)=phi_I*phi_J
DO I=1,1+Order
   DO J=1,1+Order
      PhiIPhiJ(I,J)=sum(Phi(I,:)*Phi(J,:)*wGL1D)
   END DO
END DO
!initialize h to avoid some problems :-)
h=sqrt((Coor(Tri(1,2),1)-Coor(Tri(1,1),1))**2+(Coor(Tri(1,2),2)-Coor(Tri(1,1),2))**2)

DO I=1,Nflu+Nflusol
   Node = Tri(I,:)
   !!Computation of vectors V12,V23 and V31
   V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
   V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
   V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
   V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
   V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
   V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)
 CALL sub_ray_circle_ins(Coor(Node(1),:),Coor(Node(2),:),Coor(Node(3),:),h1) 
!!$   h1=dsqrt(sum(V(1,:)**2))
!!$   h1=min(h1,dsqrt(sum(V(2,:)**2)))
!!$   h1=min(h1,dsqrt(sum(V(3,:)**2)))
   h=min(h,h1)
   !write(6,*) h,Veloc(I)
   CFL=min(CFL,h/dsqrt(mu(I)/rho(I)))
   !write(6,*) CFL
   !!Computation of normal vectors n1,n2 and n3
   norm(1,1)=V(1,2)
   norm(1,2)=-V(1,1)
   norm(1,:)=norm(1,:)/dsqrt(norm(1,1)**2+norm(1,2)**2)
   norm(2,1)=V(2,2)
   norm(2,2)=-V(2,1)
   norm(2,:)=norm(2,:)/dsqrt(norm(2,1)**2+norm(2,2)**2)
   norm(3,1)=V(3,2)
   norm(3,2)=-V(3,1)
   norm(3,:)=norm(3,:)/dsqrt(norm(3,1)**2+norm(3,2)**2)
   
   !! Computation of JFinv and DF
   DFVEC(I)=-V(3,1)*V(2,2)+V(3,2)*V(2,1)
   
   if(Dfvec(I).lt.0.D0) then
      write(6,*) 'Error, Dfvec(I)<0'
      stop
   end if
!!! Attention : Jfinv est la transposée de JF^-1
     Jfinv(1,1)=-V(2,2)
     Jfinv(2,2)=V(3,1)
     Jfinv(1,2)=-V(3,2)
     Jfinv(2,1)=V(2,1)
     Jfinv=(Jfinv/DFVEC(I))

!! Computation of int(grad u .grad v)
     DO J=1,Nphi
        DO K=1,Nphi
           res=(Jfinv(1,1)**2+Jfinv(2,1)**2)*GradPhiIGradPhiJ(J,K,1,1)
           res=res+(Jfinv(1,1)*Jfinv(1,2)+Jfinv(2,1)*Jfinv(2,2))*(GradPhiIGradPhiJ(J,K,1,2)+GradPhiIGradPhiJ(J,K,2,1))
           res=res+(Jfinv(1,2)**2+Jfinv(2,2)**2)*GradPhiIGradPhiJ(J,K,2,2)
           A_pp(I,J,K)=res*DFVEC(I)/rho(I)
        END DO
     END DO
     L1=Nphi
     DO Arete=1,3
!!! Searching Neighbor
     INeigh=Neigh(I,Arete)
!!! If there is a Neighbor
     IF((Ineigh.gt.0).and.(Ineigh.le.Nflu+Nflusol)) then
     V_voisin(1,1)=Coor(Tri(Ineigh,3),1)-Coor(Tri(Ineigh,2),1)
     V_voisin(1,2)=Coor(Tri(Ineigh,3),2)-Coor(Tri(Ineigh,2),2)
     V_voisin(2,1)=Coor(Tri(Ineigh,1),1)-Coor(Tri(Ineigh,3),1)
     V_voisin(2,2)=Coor(Tri(Ineigh,1),2)-Coor(Tri(Ineigh,3),2)
     V_voisin(3,1)=Coor(Tri(Ineigh,2),1)-Coor(Tri(Ineigh,1),1)
     V_voisin(3,2)=Coor(Tri(Ineigh,2),2)-Coor(Tri(Ineigh,1),2)


!!! Compute hmin=min(h)
 CALL sub_ray_circle_ins(Coor(Tri(INeigh,1),:),Coor(Tri(INeigh,2),:),Coor(Tri(INeigh,3),:),h2)
!!$     h2=dsqrt(sum(V_voisin(1,:)**2))
!!$     h2=min(h2,dsqrt(sum(V_voisin(2,:)**2)))
!!$     h2=min(h2,dsqrt(sum(V_voisin(3,:)**2)))
     hmin=min(h1,h2)
     cmax=max(1D0/rho(I),1D0/rho(Ineigh))




!!! Attention : Jfinv est la transposée de JF^-1
     Jfinv_voisin(1,1)=-V_voisin(2,2)
     Jfinv_voisin(2,2)=V_voisin(3,1)
     Jfinv_voisin(1,2)=-V_voisin(3,2)
     Jfinv_voisin(2,1)=V_voisin(2,1)
     Jfinv_voisin=Jfinv_voisin/abs(-V_voisin(3,1)*V_voisin(2,2)+V_voisin(3,2)*V_voisin(2,1))


        IF (Neigh(INeigh,1).eq.I) THEN
           arete_neigh=1
        ELSE IF (Neigh(Ineigh,2).eq.I) THEN
           arete_neigh=2
        ELSE IF (Neigh(Ineigh,3).eq.I) THEN
           arete_neigh=3
        ENDIF



!! Computation of -{grad u}. [v]-{grad v}.[u]
        DO J=1,Nphi
           DO K=1,Order+1
              res=-sum((MATMUL(Jfinv,GradPhiIPhiJ(Arete,J,K,:)))*norm(Arete,:))
              L=PhiEdge(Arete,K)
!grad(phi_I,J)*phi_I,L/2
              A_pp(I,J,L)=A_pp(I,J,L)+res*sqrt(sum(V(Arete,:)**2))/(2.D0*rho(I))
!grad(phi_I,L)*phi_I,J/2
              A_pp(I,L,J)=A_pp(I,L,J)+res*sqrt(sum(V(Arete,:)**2))/(2.D0*rho(I))

              L=PhiEdge(Arete_neigh,2+Order-K)
!grad(phi_I,J)*phi_Ineigh,L/2 
!(attention, sur l'element voisin on tourne en sens inverse, d'ou le 2-Order-K
              A_pp(I,J,Arete*Nphi+L)=A_pp(I,J,Arete*Nphi+L)-res*sqrt(sum(V(Arete,:)**2))/(2D0*rho(I))
              res=-sum((MATMUL(Jfinv_voisin,GradPhiIPhiJ(Arete_neigh,J,K,:)))*norm(Arete,:))
              L=PhiEdge(Arete,2+Order-K)
!grad(phi_Ineigh,J)*phi_I,L/2 
!(attention, sur l'element voisin on tourne en sens inverse, d'ou le 2-Order-K
              A_pp(I,L,Arete*Nphi+J)=A_pp(I,L,Arete*Nphi+J)+res*sqrt(sum(V(Arete,:)**2))/(2.D0*rho(Ineigh))

           END DO
        END DO
!!alpha/hmin [u][v]
        DO J=1,Order+1
           DO K=1,Order+1
              res=PhiIPhiJ(J,K)*alpha*cmax/hmin
              L=PhiEdge(Arete,K)
              Q=PhiEdge(Arete,J)
!phi_I,Q*phi_I,L
              A_pp(I,Q,L)=A_pp(I,Q,L)+res*sqrt(sum(V(Arete,:)**2))
!phi_I,Q*phi_Ineigh,L
              L=PhiEdge(Arete_neigh,2+Order-K)
              A_pp(I,Q,Arete*Nphi+L)=A_pp(I,Q,Arete*Nphi+L)-res*sqrt(sum(V(Arete,:)**2))
           END DO
        END DO
     L1=L1+Nphi
!!! Dirichlet 
     elseif(Ineigh.eq.-2) then
!!! Computation of -{grad u}. [v]-{grad v}.[u]
        DO J=1,Nphi
           DO K=1,Order+1
              L=PhiEdge(Arete,K)
              res=-sum((MATMUL(Jfinv,GradPhiIPhiJ(Arete,J,K,:)))*norm(Arete,:))
!grad(phi_I,J)*phi_I,L/2
             A_pp(I,J,L)=A_pp(I,J,L)+res*sqrt(sum(V(Arete,:)**2))/rho(I)
!grad(phi_I,L)*phi_I,J/2
              A_pp(I,L,J)=A_pp(I,L,J)+res*sqrt(sum(V(Arete,:)**2))/rho(I)
           END DO
        END DO
!!alpha/hmin [u][v]
        DO J=1,Order+1
           DO K=1,Order+1
              res=PhiIPhiJ(J,K)*2.D0*alpha/(h1*rho(I))
              L=PhiEdge(Arete,K)
              Q=PhiEdge(Arete,J)
!phi_I,Q*phi_I,L
              A_pp(I,Q,L)=A_pp(I,Q,L)+res*sqrt(sum(V(Arete,:)**2))
           END DO
        END DO
     END IF
  END DO
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!    
!!!!!!!!! Computation of Minv*A
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (helmholtz.eq.0)then
DO I=1,Nflu+Nflusol
   A_pp(I,1:Nphi,1:4*Nphi)=Matmul(Minv,A_pp(I,1:Nphi,1:4*Nphi))/DFVec(I)*mu(I)
ENDDO
end if
END SUBROUTINE sub_defstiffmat_elasto_acoustic_pp
