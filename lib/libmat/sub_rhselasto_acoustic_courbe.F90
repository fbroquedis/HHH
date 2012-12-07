subroutine sub_rhselasto_acoustic_courbe(F,theta0,R0)
  Use m_condbord    
  Use m_mat
  Use m_gen
  Use m_mesh
implicit none 
real*8,dimension(:),allocatable :: vectJ
complex*16, dimension(:), allocatable ::Valsourcebord
complex*16, dimension(:,:), allocatable ::Gradvalsourcebord
complex*16 :: GradPhiIG(2),res
real*8 :: V(3,2),Jfinv(2,2),thetai(2),norm(3,2),h2,hmin,coor_x,coor_y,theta
real*8 :: coorx1,coorx2,coory1,coory2,cmax,theta1,theta2,R0,normal(2)
integer :: Node(3),I,int_j,J,l,NsommeNum,K,NCALC,Ineigh,PhiEdge(3,1+Order),Nsomme,J1,J2


real*8 :: demi_grand_axe_a,demi_petit_axe_b,h1,ksi,theta0,foc,pi

REAL*8 :: Phi(1+Order,NGL1D),coord_points(4)

complex*16 :: F((Nflu+Nflusol+2*Nsol+2*Nsolflu)*Nphi,1)
REAL*8:: GradPhi1D(order+1,NGL1D),d,pt_courbe_x(Order+1),pt_courbe_y(Order+1),dtheta
	pi=2.d0*dasin(1.d0)
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

pi=2.d0*dasin(1.d0)
	
	allocate(Valsourcebord(NGL1D))
	allocate(GradValsourcebord(NGL1D,2))
	
	
	
 	
	

  Do J=Nflu+1,Nflu+Nflusol
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
      IF (Ineigh.gt.Nflu+Nflusol) then
        h2=sqrt((Coor(Tri(Ineigh,2),1)-Coor(Tri(Ineigh,1),1))**2+(Coor(Tri(Ineigh,2),2)-Coor(Tri(Ineigh,1),2))**2)
        h2=min(h2,sqrt((Coor(Tri(Ineigh,3),1)-Coor(Tri(Ineigh,1),1))**2+(Coor(Tri(Ineigh,3),2)-Coor(Tri(Ineigh,1),2))**2))
        h2=min(h2,sqrt((Coor(Tri(Ineigh,2),1)-Coor(Tri(Ineigh,3),1))**2+(Coor(Tri(Ineigh,2),2)-Coor(Tri(Ineigh,3),2))**2))
        hmin=min(h1,h2)
         SELECT CASE(int_j)
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
           do l=1,NGL1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!On calcule grad p_incident 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              GradValsourcebord(l,1) = +cmplx(0,omega/sqrt(mu(1)/rho(1)))*cos(theta0)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &sum(pt_courbe_x*phi(:,l))*Dcos(theta0)+&
                   &sum(pt_courbe_y*phi(:,l))*dsin(theta0)))))
              GradValsourcebord(l,2) = +cmplx(0,omega/dsqrt(mu(1)/rho(1)))*sin(theta0)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &sum(pt_courbe_x*phi(:,l))*dcos(theta0)+&
                   &sum(pt_courbe_y*phi(:,l))*dsin(theta0)))))
           enddo
	        do K=1,Order+1 
!!! Calcul du terme grad uinc*phiK
res=0D0
             DO L=1,NGL1D
!                 d=sum(pt_courbe_x*gradphi1D(:,L))**2
!                 d=d+sum(pt_courbe_y*gradphi1D(:,L))**2
                normal(1)=sum(pt_courbe_y*gradphi1D(:,L))
                normal(2)=-sum(pt_courbe_x*gradphi1D(:,L))
!                normal=normal/dsqrt(d)
                res=res+Phi(K,l)*wGL1D(l)*&
                     &(GradValsourcebord(l,1)*normal(1)+GradValsourcebord(l,2)*normal(2))
         end DO
               F((J-1)*Nphi+PhiEdge(int_j,K),1) = F((J-1)*Nphi+PhiEdge(int_j,K),1) -&
                    &res/rho(J)

            enddo
         end IF

      enddo
   enddo

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
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! On calcule p_incident
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           do l=1,NGL1D
              Valsourcebord(l) = exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &sum(pt_courbe_x*phi(:,l))*dcos(theta0)+&
                   &sum(pt_courbe_y*phi(:,l))*dsin(theta0)))))
           enddo
!!! Calcul du terme uinc*gradphiK
I=(2*(J-1)-Nflu-Nflusol)*Nphi
           do K=1,Order+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Attention
!!!! On a besoin de p_incident (la dérivée en temps)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de  -p_incident*n_x 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
res=0D0
             DO L=1,NGL1D
!                 d=sum(pt_courbe_x*gradphi1D(:,L))**2
!                 d=d+sum(pt_courbe_y*gradphi1D(:,L))**2
                normal(1)=sum(pt_courbe_y*gradphi1D(:,L))
                normal(2)=-sum(pt_courbe_x*gradphi1D(:,L))
!                normal=normal/dsqrt(d)
               res=res+Phi(K,l)*Valsourcebord(l)*wGL1D(l)*normal(1)
         end DO
F(I+PhiEdge(int_j,K),1) = F(I+PhiEdge(int_j,K),1)&
     &-res
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de  -p_incident*n_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
res=0D0
             DO L=1,NGL1D
!                 d=sum(pt_courbe_x*gradphi1D(:,L))**2
!                 d=d+sum(pt_courbe_y*gradphi1D(:,L))**2
                normal(1)=sum(pt_courbe_y*gradphi1D(:,L))
                normal(2)=-sum(pt_courbe_x*gradphi1D(:,L))
!                normal=normal/dsqrt(d)
               res=res+Phi(K,l)*Valsourcebord(l)*wGL1D(l)*normal(2)
         end DO
F(I+Nphi+PhiEdge(int_j,K),1) = F(I+Nphi+PhiEdge(int_j,K),1)&
     &-res
              end do
           end IF

        enddo
     enddo
  deallocate(Valsourcebord)

end subroutine sub_rhselasto_acoustic_courbe
