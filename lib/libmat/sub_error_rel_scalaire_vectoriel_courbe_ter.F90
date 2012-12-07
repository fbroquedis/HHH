subroutine sub_error_abs_scalaire_vectoriel_courbe3(error_p, error_ux, error_uy,R0)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K,L,J1,J2,indice,edge,node(3),test
  INTEGER,allocatable :: PhiEdge(:,:)
	real*8 :: error_p, error_ux, error_uy,pi,d,dtheta,pt_courbe_x(Nphi),pt_courbe_y(Nphi)
 REAL*8,allocatable::Phi2D(:,:),M_courbe(:,:)
  REAL*8,allocatable :: GradPhi2D(:,:,:)
  REAL*8 :: DF(2,2,NGL2D),DFVEC_bis(NGL2D),theta1,theta2,R0,vectmp1(2),vectmp2(2),DFINV(2,2,NGL2D)
 pi=2.d0*dasin(1.d0)
  allocate(GradPhi2D(Nphi,2,NGL2D))
  allocate(PhiEdge(3,1+Order))
  allocate(Phi2D(Nphi,NGL2D))
  allocate(M_courbe(Nphi,Nphi))
  SELECT CASE(Order)
  CASE(1)
     !Functions on Edge 23
     PhiEdge(1,1)=2
     PhiEdge(1,2)=3
     !Functions on Edge 31
     PhiEdge(2,1)=3
     PhiEdge(2,2)=1
     !Functions on Edge 12
     PhiEdge(3,1)=1
     PhiEdge(3,2)=2
     CALL GradPhiOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     
 

  CASE(2)
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
     CALL GradPhiOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     

  CASE(3)!Functions on Edge 23
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

     CALL GradPhiOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     
  END SELECT	
	!! Uold est la solution exacte et U la solution approchée
	error_p=0.d0
	error_ux=0.d0
	error_uy=0.d0

 Do I=1,Nflu+Nflusol
 Node = Tri(I,:)
pt_courbe_x(1:3)=Coor(Node(:),1)
pt_courbe_y(1:3)=Coor(Node(:),2)
indice=3
test=1
DO Edge=1,3
If (neigh(I,Edge).lt.0) then
test=0
end If
SELECT CASE(Edge)
CASE(1)
!! First Edge
J1=2
J2=3
indice=3+order-1
!! Nodes of the first edge are numbered from 3+order to 3+2*order-1
CASE(2)
!! Second Edge
J1=3
J2=1
!! Nodes of the second edge are numbered from 3+2*order-1 to 3+3*order-2
indice=3+2*order-2
CASE(3)
!! Third Edge
J1=1
J2=2
!! Nodes of the second edge are numbered from 4 to 3+order-1
indice=3
END SELECT
!If (neigh(I,Edge).gt.Nflu+nflusol) then
end DO
if (test.eq.1)then
       DO K=1,Nphi
          error_p = error_p + (P_analytic((I-1)*Nphi+K)-P_complexe((I-1)*Nphi+K))*&
               &conjg(P_analytic((I-1)*Nphi+K)-P_complexe((I-1)*Nphi+K))
       ENDDO
    End if
 end Do

	Do I=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu+Nsol
 Node = Tri(I,:)
pt_courbe_x(1:3)=Coor(Node(:),1)
pt_courbe_y(1:3)=Coor(Node(:),2)
indice=3
DO Edge=1,3
SELECT CASE(Edge)
CASE(1)
!! First Edge
J1=2
J2=3
indice=3+order-1
!! Nodes of the first edge are numbered from 3+order to 3+2*order-1
CASE(2)
!! Second Edge
J1=3
J2=1
!! Nodes of the second edge are numbered from 3+2*order-1 to 3+3*order-2
indice=3+2*order-2
CASE(3)
!! Third Edge
J1=1
J2=2
!! Nodes of the second edge are numbered from 4 to 3+order-1
indice=3
END SELECT
END DO

  	L=((I-1)-Nflu-Nflusol)*Nphi
!   If (neigh(I,Edge).le.Nflu+nflusol) then
   		DO K=1,Nphi
				error_ux = error_ux + (Ux_analytic(L+K)-Ux_complexe(L+K))*&
				&conjg(Ux_analytic(L+K)-Ux_complexe(L+K))
				error_uy = error_uy + (Uy_analytic(L+K)-Uy_complexe(L+K))*&
				&conjg(Uy_analytic(L+K)-Uy_complexe(L+K))
!end If
ENDDO
End Do





	error_p=sqrt(error_p)
        error_ux=sqrt(error_ux)
        error_uy=sqrt(error_uy)

      end subroutine sub_error_abs_scalaire_vectoriel_courbe3

subroutine sub_error_rel_scalaire_vectoriel_courbe3(error_p,error_ux,error_uy,R0)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K,L,J1,J2,indice,edge,node(3),test
  INTEGER,allocatable :: PhiEdge(:,:)
	real*8 :: error_p,error_ux,error_uy,error_p_num,error_p_den,error_ux_num,error_ux_den,error_uy_num,error_uy_den
	real*8::pi,d,dtheta,pt_courbe_x(Nphi),pt_courbe_y(Nphi)
 REAL*8,allocatable::Phi2D(:,:),M_courbe(:,:)
  REAL*8,allocatable :: GradPhi2D(:,:,:)
  REAL*8 :: DF(2,2,NGL2D),DFVEC_bis(NGL2D),theta1,theta2,R0,vectmp1(2),vectmp2(2),DFINV(2,2,NGL2D)
 pi=2.d0*dasin(1.d0)
 allocate(PhiEdge(3,1+Order))
  allocate(GradPhi2D(Nphi,2,NGL2D))
  allocate(Phi2D(Nphi,NGL2D))
  allocate(M_courbe(Nphi,Nphi))
  SELECT CASE(Order)
  CASE(1)
     !Functions on Edge 23
     PhiEdge(1,1)=2
     PhiEdge(1,2)=3
     !Functions on Edge 31
     PhiEdge(2,1)=3
     PhiEdge(2,2)=1
     !Functions on Edge 12
     PhiEdge(3,1)=1
     PhiEdge(3,2)=2
     CALL GradPhiOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     
 

  CASE(2)
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
     CALL GradPhiOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     

  CASE(3)
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
     CALL GradPhiOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     
  END SELECT	
	!! P_complexe est la solution exacte et P_analytic la solution approchée
	call sub_error_abs_scalaire_vectoriel_courbe2(error_p_num,error_ux_num,error_uy_num,R0)
	error_p_den=0.d0
	error_ux_den=0.d0
	error_uy_den=0.d0


 Do I=1,Nflu+Nflusol 
 Node = Tri(I,:)
pt_courbe_x(1:3)=Coor(Node(:),1)
pt_courbe_y(1:3)=Coor(Node(:),2)
indice=3
test=1
DO Edge=1,3
SELECT CASE(Edge)
CASE(1)
!! First Edge
J1=2
J2=3
indice=3+order-1
!! Nodes of the first edge are numbered from 3+order to 3+2*order-1
CASE(2)
!! Second Edge
J1=3
J2=1
!! Nodes of the second edge are numbered from 3+2*order-1 to 3+3*order-2
indice=3+2*order-2
CASE(3)
!! Third Edge
J1=1
J2=2
!! Nodes of the second edge are numbered from 4 to 3+order-1
indice=3
END SELECT
If (neigh(I,Edge).lt.0) then
test=0
end If
end DO
!If (neigh(I,Edge).gt.Nflu+nflusol) then
If (test.eq.1) then
DO K=1,Nphi
				error_p_den = error_p_den + (P_analytic((I-1)*Nphi+K))*&
				&conjg(P_analytic((I-1)*Nphi+K))
 ENDDO
end If
End Do

	error_p=error_p_num/sqrt(error_p_den)*100

	Do I=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu+Nsol
 Node = Tri(I,:)
pt_courbe_x(1:3)=Coor(Node(:),1)
pt_courbe_y(1:3)=Coor(Node(:),2)
indice=3
DO Edge=1,3
SELECT CASE(Edge)
CASE(1)
!! First Edge
J1=2
J2=3
indice=3+order-1
!! Nodes of the first edge are numbered from 3+order to 3+2*order-1
CASE(2)
!! Second Edge
J1=3
J2=1
!! Nodes of the second edge are numbered from 3+2*order-1 to 3+3*order-2
indice=3+2*order-2
CASE(3)
!! Third Edge
J1=1
J2=2
!! Nodes of the second edge are numbered from 4 to 3+order-1
indice=3
END SELECT
end DO
L=((I-1)-Nflu-Nflusol)*Nphi
!If (neigh(I,Edge).le.Nflu+nflusol) then
DO K=1,Nphi

				error_ux_den = error_ux_den + (Ux_analytic(L+K))*&
				&conjg(Ux_analytic(L+K))
				error_uy_den = error_uy_den + (Uy_analytic(L+K))*&
				&conjg(Uy_analytic(L+K))
 ENDDO
!End if
enddo

 error_ux=error_ux_num/sqrt(error_ux_den)*100
 error_uy=error_uy_num/sqrt(error_uy_den)*100	
 

end subroutine sub_error_rel_scalaire_vectoriel_courbe3
