subroutine sub_rhsacoustique(F,Fcomp,theta0)
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
real*8 :: coorx1,coorx2,coory1,coory2
integer :: Node(3),I,int_j,J,l,NsommeNum,K,NCALC,Ineigh,PhiEdge(3,1+Order),Nsomme


real*8 :: demi_grand_axe_a,demi_petit_axe_b,h1,ksi,theta0,foc,pi

REAL*8 :: Phi(1+Order,NGL1D),coord_points(4)

!!!!complex*16,dimension(Ntri*Nphi,1) :: F
real*8 :: F(Nphi*Ntri,1),Fcomp(Nphi*Ntri,1)
REAL*8:: GradPhi1D(3,Nphi,2,NGL1D),d
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
write(6,*) 'toto'
write(6,*) 'ok dirich?'
	
	allocate(Valsourcebord(NGL1D))
	allocate(GradValsourcebord(NGL1D,2))
	
	
	
 NCALC=0
	NbR0=0
  DO I=1,NTri
     DO int_j=1,3
        INeigh=Neigh(I,int_j)		    
        !! Si le triangle est sur l'interface
      IF (Ineigh.gt.0) then 
         if (((type_media(I).eq.1).and.(type_media(INeigh).ne.1)).or.&
              &((type_media(I).ne.1).and.(type_media(INeigh).eq.1))) THEN
           NbR0=NbR0+1
        endif
     end IF
     enddo

  enddo

  write(6,*) 'nbr0',NbR0

  allocate(TriBordR0(NbR0))
  NbR0=0
  DO I=1,NTri
     DO int_j=1,3
        INeigh=Neigh(I,int_j)		
            
        !Dirichlet
        !! Si le triangle est sur l'interface
      IF (Ineigh.gt.0) then
         if (((type_media(I).eq.1).and.(type_media(INeigh).ne.1)).or.&
              &((type_media(I).ne.1).and.(type_media(INeigh).eq.1))) THEN
           NbR0=NbR0+1
            TriBordR0(NbR0)=I
         endif
      end IF
     enddo
  enddo
  
  
  
  
	
	
  !Dirichlet

  Do I=1,NbR0
     J = TriBordR0(I)
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
      IF (Ineigh.gt.0) then
         if (((type_media(J).eq.1).and.(type_media(INeigh).ne.1)).or.&
              &((type_media(J).ne.1).and.(type_media(INeigh).eq.1))) THEN
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
              Valsourcebord(l) = -exp(cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))
              GradValsourcebord(l,1) = -cmplx(0,omega)*cos(theta0)/sqrt(mu(1)/rho(1))&
                   &*exp(cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))
              GradValsourcebord(l,2) = -cmplx(0,omega)*sin(theta0)/sqrt(mu(1)/rho(1))&
                   &*exp(cmplx(0,(omega/sqrt(mu(1)/rho(1))*(&
                   &((coorx2-coorx1)*PtGl1D(l)+coorx1)*cos(theta0)+&
                   &((coory2-coory1)*PtGl1D(l)+coory1)*sin(theta0)))))
           enddo
	
!!! Calcul du terme uinc*gradphiK
           do K=1,Nphi
              do l=1,2
                 GradPhiIG(l) = sum(GradPhi1D(int_j,K,l,:)*Valsourcebord(:)*wGL1D)
              enddo
              If(TYPE_MEDIA(J).eq.1) then
                 F((J-1)*Nphi+K,1) = F((J-1)*Nphi+K,1)+d*sum((matmul(Jfinv,GradPhiIG))*norm(int_j,:))/(2D0*rho(J))
              else
                 F((J-1)*Nphi+K,1) = F((J-1)*Nphi+K,1)-d*sum((matmul(Jfinv,GradPhiIG))*norm(int_j,:))/(2D0*rho(J))
              end If
      end do
	        do K=1,Order+1 
!!! Calcul du terme grad uinc*phiK
            do l=1,2
               GradPhiIG(l)=sum(Phi(K,:)*GradValsourcebord(:,l)*wGL1D)
            end do
            If(Type_Media(J).eq.1) then
               F((J-1)*Nphi+PhiEdge(int_j,K),1) = F((J-1)*Nphi+PhiEdge(int_j,K),1) -&
                    &d*sum(GradPhiIG*norm(int_j,:))/(2D0*rho(1))
            else
               F((J-1)*Nphi+PhiEdge(int_j,K),1) = F((J-1)*Nphi+PhiEdge(int_j,K),1) +&
                    &d*sum(GradPhiIG*norm(int_j,:))/(2D0*rho(1))
           end If

If(Type_Media(J).eq.1) then
               F((J-1)*Nphi+PhiEdge(int_j,K),1) = F((J-1)*Nphi+PhiEdge(int_j,K),1) -&
               &alpha/hmin*d*sum(Phi(K,:)*Valsourcebord(:)*wGL1D)
else
               F((J-1)*Nphi+PhiEdge(int_j,K),1) = F((J-1)*Nphi+PhiEdge(int_j,K),1)+&
               &alpha/hmin*d*sum(Phi(K,:)*Valsourcebord(:)*wGL1D)
            end If

         enddo
      endif
   end IF
enddo
enddo
F((J-1)*Nphi+1:J*Nphi,1)=matmul(Minv,F((J-1)*Nphi+1:J*Nphi,1))/DFVEC(J)
Fcomp((J-1)*Nphi+1:J*Nphi,1)=matmul(Minv,Fcomp((J-1)*Nphi+1:J*Nphi,1))/DFVEC(J)
  deallocate(Valsourcebord)
  deallocate(TriBordR0)
end subroutine sub_rhsacoustique
