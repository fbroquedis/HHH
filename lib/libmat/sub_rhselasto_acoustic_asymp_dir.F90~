subroutine sub_rhselasto_acoustic_asymp_dir(F,theta0,epsilon)
  Use m_condbord    
  Use m_mat
  Use m_gen
  Use m_mesh
implicit none 
real*8,dimension(:),allocatable :: vectJ
complex*16, dimension(:), allocatable ::Valsourcebord
complex*16, dimension(:,:), allocatable ::Gradvalsourcebord
complex*16 :: GradPhiIG(2)
real*8 :: V(3,2),Jfinv(2,2),thetai(2),norm(3,2),h2,hmin,coor_x,coor_y,theta
real*8 :: coorx1,coorx2,coory1,coory2,cmax,epsilon
integer :: Node(3),I,int_j,J,l,NsommeNum,K,NCALC,Ineigh,PhiEdge(3,1+Order),Nsomme


real*8 :: demi_grand_axe_a,demi_petit_axe_b,h1,ksi,theta0,foc,pi

REAL*8 :: Phi(1+Order,NGL1D),coord_points(4)

complex*16 :: F((Nflu+Nflusol+2*Nsol+2*Nsolflu)*Nphi,1),res
REAL*8:: GradPhi1D(3,Nphi,2,NGL1D),d
	pi=2.d0*dasin(1.d0)
  SELECT CASE(Order)
  CASE(1)
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

pi=2.d0*dasin(1.d0)
	
	allocate(Valsourcebord(NGL1D))
	allocate(GradValsourcebord(NGL1D,2))
	
	
	
 	
	


  Do J=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu
     thetai=0.D0
     
     Node = Tri(J,:)
     
     !!Computation of vectors V12,V23 and V31
     V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
     V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
     V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
     V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
     V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
     V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)
     
      h1=sqrt(sum(V(1,:)**2))
     h1=min(h1,sqrt(sum(V(2,:)**2)))
     h1=min(h1,sqrt(sum(V(3,:)**2)))
     
     
     !!Computation of normal vectors n1,n2 and n3
     norm(1,1)=V(1,2)
     norm(1,2)=-V(1,1)
     norm(1,:)=norm(1,:)/sqrt(norm(1,1)**2+norm(1,2)**2)
     norm(2,1)=V(2,2)
     norm(2,2)=-V(2,1)
     norm(2,:)=norm(2,:)/sqrt(norm(2,1)**2+norm(2,2)**2)
     norm(3,1)=V(3,2)
     norm(3,2)=-V(3,1)
     norm(3,:)=norm(3,:)/sqrt(norm(3,1)**2+norm(3,2)**2)
     
     !! Computation of JFinv and DF
     DFVEC(J)=-V(3,1)*V(2,2)+V(3,2)*V(2,1)
!!! Attention : Jfinv est la transposée de JF^-1
     Jfinv(1,1)=-V(2,2)
     Jfinv(2,2)=V(3,1)
     Jfinv(1,2)=-V(3,2)
     Jfinv(2,1)=V(2,1)
     Jfinv=Jfinv/DFVEC(J)
     
     DO int_j=1,3
        INeigh=Neigh(J,int_j)	
   
        !! Si le triangle est sur l'interface    
      IF ((Ineigh.gt.Nflu).and.(Ineigh.le.Nflu+Nflusol)) then
        h2=sqrt((Coor(Tri(Ineigh,2),1)-Coor(Tri(Ineigh,1),1))**2+(Coor(Tri(Ineigh,2),2)-Coor(Tri(Ineigh,1),2))**2)
        h2=min(h2,sqrt((Coor(Tri(Ineigh,3),1)-Coor(Tri(Ineigh,1),1))**2+(Coor(Tri(Ineigh,3),2)-Coor(Tri(Ineigh,1),2))**2))
        h2=min(h2,sqrt((Coor(Tri(Ineigh,2),1)-Coor(Tri(Ineigh,3),1))**2+(Coor(Tri(Ineigh,2),2)-Coor(Tri(Ineigh,3),2))**2))
        hmin=min(h1,h2)
         SELECT CASE(int_j)
            CASE(1)
               coorx1=Coor(Node(2),1)
               coory1=Coor(Node(2),2)
               coorx2=Coor(Node(3),1)
               coory2=Coor(Node(3),2)
            CASE(2)
               coorx1=Coor(Node(3),1)
               coory1=Coor(Node(3),2)
               coorx2=Coor(Node(1),1)
               coory2=Coor(Node(1),2)
            CASE(3)
               coorx1=Coor(Node(1),1)
               coory1=Coor(Node(1),2)
               coorx2=Coor(Node(2),1)
               coory2=Coor(Node(2),2)
         END SELECT
           
		  
         d=sqrt(sum(V(int_j,:)**2)) 
           do l=1,NGL1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! On calcule p_incident
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              Valsourcebord(l) = exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!On calcule grad p_incident 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              GradValsourcebord(l,1) = +cmplx(0,omega/sqrt(mu(1)/rho(1)))*cos(theta0)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))
              GradValsourcebord(l,2) = +cmplx(0,omega/sqrt(mu(1)/rho(1)))*sin(theta0)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))

           enddo
I=(2*(J-1)-Nflu-Nflusol)*Nphi
           do K=1,Order+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Attention
!!!! On a besoin de p_incident (la dérivée en temps)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de  -i*omega*p_incident*n_x 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
F(I+PhiEdge(int_j,K),1) = F(I+PhiEdge(int_j,K),1)&
     &-d*sum(Phi(K,:)*Valsourcebord(:)*wGL1D)*norm(int_j,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de  -i*omega*p_incident*n_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
F(I+Nphi+PhiEdge(int_j,K),1) = F(I+Nphi+PhiEdge(int_j,K),1)&
     &-d*sum(Phi(K,:)*Valsourcebord(:)*wGL1D)*norm(int_j,2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de  -c*grad p_incident.n*nx 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
res=sum(Phi(K,:)*GradValsourcebord(:,1)*wGL1D)*norm(int_j,1)
res=res+sum(Phi(K,:)*GradValsourcebord(:,2)*wGL1D)*norm(int_j,2)
F(I+PhiEdge(int_j,K),1) = F(I+PhiEdge(int_j,K),1)&
     &-d*res*norm(int_j,1)*epsilon*(1D0-epsilon/2D0/0.01D0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de  -c*grad p_incident.n*n_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
F(I+Nphi+PhiEdge(int_j,K),1) = F(I+Nphi+PhiEdge(int_j,K),1)&
     &-d*res*norm(int_j,2)*epsilon*(1D0-epsilon/2D0/0.01D0)
              end do

           end IF

        enddo
     enddo
  deallocate(Valsourcebord)
end subroutine sub_rhselasto_acoustic_asymp_dir
