subroutine sub_rhshelmholtz_scalaire_vectoriel_exact(F,theta0)
  Use m_condbord    
  Use m_mat
  Use m_gen
  Use m_mesh
  implicit none 
  real*8,dimension(:),allocatable :: vectJ
  complex*16, dimension(:), allocatable ::Valsourcebord, Lapvalsourcebord
  complex*16, dimension(:,:), allocatable ::Gradvalsourcebord
  complex*16 :: GradPhiIG(2)
  real*8 :: V(3,2),Jfinv(2,2),thetai(2),norm(3,2),h2,hmin,coor_x,coor_y,theta
  real*8 :: coorx1,coorx2,coory1,coory2,cmax
  integer :: Node(3),I,int_j,J,l,NsommeNum,K,NCALC,Ineigh,PhiEdge(3,1+Order),Nsommen,Isol


  real*8 :: demi_grand_axe_a,demi_petit_axe_b,h1,ksi,theta0,foc,pi

  REAL*8 :: Phi(1+Order,NGL1D),coord_points(4)

  complex*16 :: F((Nflu+Nflusol+2*Nsol+2*Nsolflu)*Nphi,1)
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

  allocate(LapValsourcebord(NGL1D))
  allocate(Valsourcebord(NGL1D))
  allocate(GradValsourcebord(NGL1D,2))

  NCALC=0
  NbR0=0
  DO I=1,NTri
     DO int_j=1,3
        INeigh=Neigh(I,int_j)		    
        !! Si le triangle est sur l'interface fluide
        IF (Ineigh.eq.-3) then 
           NbR0=NbR0+1
           exit
        end IF
     enddo

  enddo


  !! Calcul du terme source f3 provenant de bord extérieur
  !! int dp/dn*q avec p onde plane

  !!Evaluation des mailles concernées

  allocate(TriBordR0(NbR0))
  NbR0=0
  DO I=1,Nflu
     DO int_j=1,3
        INeigh=Neigh(I,int_j)		

        !Dirichlet
        !! Si le triangle est sur l'interface
        IF (Ineigh.eq.-3) then
           NbR0=NbR0+1
           TriBordR0(NbR0)=I
           exit
        end IF
     enddo
  enddo

  !! Calcul de l'intégrale

  Do I=1,NbR0
     J = TriBordR0(I)
     theta=0.D0

     Node = Tri(J,:)

     !!Computation of vectors V12,V23 and V31
     V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
     V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
     V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
     V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
     V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
     V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)


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


     DO int_j=1,3
        INeigh=Neigh(J,int_j)	

        !! Si le triangle est sur l'interface extérieure du fluide
        IF (Ineigh.eq.-3) then

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!On calcule grad p
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
           do K=1,Order+1 
!!! Calcul du terme grad p*phiK
              do l=1,2
                 GradPhiIG(l)=sum(Phi(K,:)*GradValsourcebord(:,l)*wGL1D)
              end do

              F((J-1)*Nphi+PhiEdge(int_j,K),1) = F((J-1)*Nphi+PhiEdge(int_j,K),1) +&
                   &d*sum(GradPhiIG*norm(int_j,:))/rho(J)

           enddo
        end IF

     enddo
  enddo



  !!Calcul du terme source provenant de la première équation 
  !! int omega**2 rho_f * un + grad(p).n avec u=grad phi, où phi, et p ondes planes
  Isol=Nflu+Nflusol+1
  Do J=Nflu+1,Nflu+Nflusol
     theta=0.D0

     Node = Tri(J,:)

     !!Computation of vectors V12,V23 and V31
     V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
     V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
     V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
     V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
     V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
     V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)

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

     DO int_j=1,3
        INeigh=Neigh(J,int_j)	

        !! Si le triangle est sur l'interface fluide-solide
        IF (Ineigh.gt.Nflu+Nflusol) then

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!On calcule  w^2*rho_f*u 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              GradValsourcebord(l,1) = omega**2*rho(J)*cmplx(0,omega/sqrt(mu(Isol)/rho(Isol)))*cos(theta)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(Isol)/rho(Isol))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta)))))
              GradValsourcebord(l,2) = omega**2*rho(J)*cmplx(0,omega/sqrt(mu(Isol)/rho(Isol)))*sin(theta)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(Isol)/rho(Isol))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta)))))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!On calcule - grad p 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              GradValsourcebord(l,1) = GradValsourcebord(l,1) - cmplx(0,omega/sqrt(mu(1)/rho(1)))*cos(theta0)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))
              GradValsourcebord(l,2) =  GradValsourcebord(l,2) - cmplx(0,omega/sqrt(mu(1)/rho(1)))*sin(theta0)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))
           enddo
           do K=1,Order+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calcul du terme -(w^2*rho_f*u -grad p)*n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do l=1,2
                 GradPhiIG(l)=sum(Phi(K,:)*GradValsourcebord(:,l)*wGL1D)
              end do

              F((J-1)*Nphi+PhiEdge(int_j,K),1) = F((J-1)*Nphi+PhiEdge(int_j,K),1) -&
                   &d*sum(GradPhiIG*norm(int_j,:))/rho(J)
           end do
        end IF

     enddo
  enddo

  !!Calcul du terme source provenant de la seconde équation 
  !! int div(u)n + pn avec u=grad psi, où psi, et p ondes planes

  Do J=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu
     theta=0.D0

     Node = Tri(J,:)

     !!Computation of vectors V12,V23 and V31
     V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
     V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
     V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
     V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
     V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
     V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)


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


     DO int_j=1,3
        INeigh=Neigh(J,int_j)	

        !! Si le triangle est sur l'interface fluide-solide
        IF ((Ineigh.gt.Nflu).and.(Ineigh.le.Nflu+Nflusol)) then

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!On calcule  div(u)=lap(psi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              GradValsourcebord(l,1) = +cmplx(0,omega/sqrt(mu(Isol)/rho(Isol)))*cos(theta)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(Isol)/rho(Isol))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta)))))
              GradValsourcebord(l,2) = +cmplx(0,omega/sqrt(mu(Isol)/rho(Isol)))*sin(theta)&
                   &*exp(+cmplx(0,(omega/sqrt(mu(Isol)/rho(Isol))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta)))))

              LapValsourcebord(l) = cmplx(0,omega/sqrt(mu(Isol)/rho(Isol)))*cos(theta)*GradValsourcebord(l,1)
              LapValsourcebord(l) = LapValsourcebord(l)+ cmplx(0,omega/sqrt(mu(Isol)/rho(Isol)))*sin(theta)*GradValsourcebord(l,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! On calcule p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              Valsourcebord(l) = exp(+cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))
           enddo

           I=(2*(J-1)-Nflu-Nflusol)*Nphi
           do K=1,Order+1 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de   mu*div(u)*n_x*v_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              F(I+PhiEdge(int_j,K),1) = F(I+PhiEdge(int_j,K),1)&
                   &+d*sum(Phi(K,:)*LapValsourcebord(:)*wGL1D)*norm(int_j,1)*mu(J)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de   mu* div(u)*n_y*v_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              F(I+Nphi+PhiEdge(int_j,K),1) = F(I+Nphi+PhiEdge(int_j,K),1)&
                   &+d*sum(Phi(K,:)*LapValsourcebord(:)*wGL1D)*norm(int_j,2)*mu(J)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de  p*n_x*v_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              F(I+PhiEdge(int_j,K),1) = F(I+PhiEdge(int_j,K),1)&
                   &+d*sum(Phi(K,:)*Valsourcebord(:)*wGL1D)*norm(int_j,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Calcul de  p*n_y*v_y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              F(I+Nphi+PhiEdge(int_j,K),1) = F(I+Nphi+PhiEdge(int_j,K),1)&
                   &+d*sum(Phi(K,:)*Valsourcebord(:)*wGL1D)*norm(int_j,2)

           enddo
        end IF

     enddo
  enddo

  deallocate(LapValsourcebord)
  deallocate(GradValsourcebord)
  deallocate(Valsourcebord)
end subroutine sub_rhshelmholtz_scalaire_vectoriel_exact
