SUBROUTINE sub_defstiffmat_elasto_acoustic_asymp_dir_uxux_courbe(epsilon,R0)
  Use m_mat
  Use m_gen
  Use m_mesh
  implicit none
  INTEGER :: I,INeigh,Node(3),J,K,L,L1,Q,Arete,Arete_neigh,Iloc,J1,J2
  INTEGER,allocatable :: PhiEdge(:,:)
  REAL*8 :: Jfinv(2,2),V(3,2),Test(3),cmax,V_voisin(3,2)
  real*8 :: Jfinv_voisin(2,2),normal(2)
  REAL*8,allocatable::GradPhiIGradPhiJ(:,:,:,:)
  REAL*8,allocatable::GradPhiIPhiJ(:,:,:,:)
  REAL*8,allocatable::PhiIPhiJ(:,:)
  REAL*8,allocatable :: Phi(:,:),Ainter(:,:,:)
  REAL*8,allocatable :: GradPhi1D(:,:),GradPhi2D(:,:,:)
  REAL*8 :: norm(3,2),h1,h2,hmin,epsilon,R0,theta1,theta2,d,pi,pt_courbe_x(Order+1),pt_courbe_y(Order+1),dtheta
  complex*16 :: res
pi=2.d0*dasin(1.d0)
  allocate(Phi(1+Order,NGL1D))
  allocate(PhiEdge(3,1+Order))
  allocate(GradPhi1D(Order+1,NGL1D))
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
! GradPhiIPhiJ(I,J)=phi_I*phi_J
DO I=1,1+Order
   DO J=1,1+Order
      PhiIPhiJ(I,J)=sum(Phi(I,:)*Phi(J,:)*wGL1D)
   END DO
END DO
!initialize h to avoid some problems :-)
h=sqrt((Coor(Tri(1,2),1)-Coor(Tri(1,1),1))**2+(Coor(Tri(1,2),2)-Coor(Tri(1,1),2))**2)

DO Iloc=1,Nsolflu
I=Nflu+Nflusol+Iloc
   Node = Tri(I,:)
   !!Computation of vectors V12,V23 and V31
   V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
   V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
   V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
   V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
   V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
   V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)
   h1=dsqrt(sum(V(1,:)**2))
   h1=min(h1,dsqrt(sum(V(2,:)**2)))
   h1=min(h1,dsqrt(sum(V(3,:)**2)))
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

     L1=Nphi
     DO Arete=1,3
!!! Searching Neighbor
     INeigh=Neigh(I,Arete)
!!! If the Neighbor is in the fluid part
     IF((Ineigh.gt.0).and.(Ineigh.le.Nflu+Nflusol)) then
        if(Helmholtz.eq.0) then
           Corres_solflu(Iloc)=Arete
        endif
     V_voisin(1,1)=Coor(Tri(Ineigh,3),1)-Coor(Tri(Ineigh,2),1)
     V_voisin(1,2)=Coor(Tri(Ineigh,3),2)-Coor(Tri(Ineigh,2),2)
     V_voisin(2,1)=Coor(Tri(Ineigh,1),1)-Coor(Tri(Ineigh,3),1)
     V_voisin(2,2)=Coor(Tri(Ineigh,1),2)-Coor(Tri(Ineigh,3),2)
     V_voisin(3,1)=Coor(Tri(Ineigh,2),1)-Coor(Tri(Ineigh,1),1)
     V_voisin(3,2)=Coor(Tri(Ineigh,2),2)-Coor(Tri(Ineigh,1),2)


!!! Compute hmin=min(h)
     h2=dsqrt(sum(V_voisin(1,:)**2))
     h2=min(h2,dsqrt(sum(V_voisin(2,:)**2)))
     h2=min(h2,dsqrt(sum(V_voisin(3,:)**2)))
     hmin=min(h1,h2)
     cmax=max(mu(I),mu(Ineigh))




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


SELECT CASE(Arete)
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
      pt_courbe_x(k)= R0*dcos(theta1+(k-1)*dtheta/real(order,8))
      pt_courbe_y(k)= R0*dsin(theta1+(k-1)*dtheta/real(order,8))
   ENDDO

!! Computation of p*u_x
        DO J=1,Order+1
           DO K=1,Order+1
              res=0D0
              DO L=1,NGL1D
!                 d=sum(pt_courbe_x*gradphi1D(:,L))**2
!                 d=d+sum(pt_courbe_y*gradphi1D(:,L))**2
                normal(1)=sum(pt_courbe_y*gradphi1D(:,L))
                normal(2)=-sum(pt_courbe_x*gradphi1D(:,L))
                d=sqrt(normal(1)**2+normal(2)**2)
                 res=res+Phi(J,L)*Phi(K,L)*wGL1D(L)*normal(1)**2/d
              ENDDO
              res=res*rho(Ineigh)
              res=res*omega**2*epsilon*(1D0-epsilon/2D0/0.01D0)
              L=PhiEdge(Arete,J)
              Q=PhiEdge(Arete,K)
!phi_I,L*phi_I_neigh,Q
If(helmholtz.eq.1) then
              A_uxux(Iloc,L,Q)=A_uxux(Iloc,L,Q)-res
else
!              A_uxp(Iloc,L,(Arete-1)*Nphi+Q)=A_uxp(Iloc,L,(Arete-1)*Nphi+Q)+res*sqrt(sum(V(Arete,:)**2))
endif
           END DO
        END DO
     END IF
  END DO
END DO
if (helmholtz.eq.0)then
DO I=1,Nsolflu
   A_uxp(I,1:Nphi,1:3*Nphi)=Matmul(Minv,A_uxp(I,1:Nphi,1:3*Nphi))/DFVec(Nflu+Nflusol+I)/rho(I+Nflu+Nflusol)
ENDDO
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!    
!!!!!!!!! Computation of Minv*A
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE sub_defstiffmat_elasto_acoustic_asymp_dir_uxux_courbe
