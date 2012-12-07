subroutine sub_error_abs_scalaire_vectoriel_courbe(error_p, error_ux, error_uy,R0)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K,L,J1,J2,indice,edge,node(3)
	real*8 :: error_p, error_ux, error_uy,pi,d,dtheta,pt_courbe_x(Nphi),pt_courbe_y(Nphi)
 REAL*8,allocatable::Phi2D(:,:),M_courbe(:,:)
  REAL*8,allocatable :: GradPhi2D(:,:,:)
  REAL*8 :: DF(2,2,NGL2D),DFVEC_bis(NGL2D),theta1,theta2,R0,vectmp1(2),vectmp2(2),DFINV(2,2,NGL2D)
 pi=2.d0*dasin(1.d0)
  allocate(GradPhi2D(Nphi,2,NGL2D))
  allocate(Phi2D(Nphi,NGL2D))
  allocate(M_courbe(Nphi,Nphi))
  SELECT CASE(Order)
  CASE(1)

     CALL GradPhiOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     
 

  CASE(2)
     CALL GradPhiOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     

  CASE(3)
     CALL GradPhiOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     
  END SELECT	
	!! Uold est la solution exacte et U la solution approchée
	error_p=0.d0
	error_ux=0.d0
	error_uy=0.d0
 Do I=1,Nflu
    DO J=1,Nphi
       DO K=1,Nphi
          error_p = error_p + (P_analytic((I-1)*Nphi+J)-P_complexe((I-1)*Nphi+J))*&
               &conjg(P_analytic((I-1)*Nphi+K)-P_complexe((I-1)*Nphi+K))*DFVEC(I)*M(J,K)
       ENDDO
    End Do
 END Do
 Do I=Nflu+1,Nflu+Nflusol
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
IF(Neigh(I,Edge).le.Nflu+Nflusol) then
   DO K=2,order
indice=indice+1
!!!      x_k=x_1+(k-1)*(x_2-x_1)/order (x1 and x2 are the two vertices)
      pt_courbe_x(indice)= pt_courbe_x(J1)+(k-1)&
           &*(pt_courbe_x(J2)-pt_courbe_x(J1))/real(order,8)
      pt_courbe_y(indice)=pt_courbe_y(J1)+(k-1)&
           &*(pt_courbe_y(J2)-pt_courbe_y(J1))/real(order,8)
   ENDDO
else
call calc_theta(pt_courbe_x(J1), pt_courbe_y(J1),theta1)
call calc_theta(pt_courbe_x(J2), pt_courbe_y(J2),theta2)
dtheta=theta2-theta1
if (dtheta.gt.pi) then
dtheta=dtheta-2*pi
elseif (dtheta.lt.-pi) then
dtheta=dtheta+2*pi
end if
   DO K=2,order
indice=indice+1
      pt_courbe_x(indice)= R0*dcos(theta1+(k-1)*dtheta/real(order,8))
      pt_courbe_y(indice)= R0*dsin(theta1+(k-1)*dtheta/real(order,8))
   ENDDO
end IF
ENDDO
IF (ORDER.EQ.3) then
   indice=10
   pt_courbe_x(indice)=sum(pt_courbe_x(1:3))/3D0
   pt_courbe_y(indice)=sum(pt_courbe_y(1:3))/3D0
ENDIF
DF=0.D0
DO J=1,NGL2D
DF(1,1,J)=sum(GradPhi2D(:,1,J)*pt_courbe_x(:))
DF(1,2,J)=sum(GradPhi2D(:,2,J)*pt_courbe_x(:))
DF(2,1,J)=sum(GradPhi2D(:,1,J)*pt_courbe_y(:))
DF(2,2,J)=sum(GradPhi2D(:,2,J)*pt_courbe_y(:))
DFVEC_BIS(J)=abs(DF(1,1,J)*DF(2,2,J)-DF(1,2,J)*DF(2,1,J))
!!! Attention, DFINV est la transposée de JFinv
DFINV(1,1,J)=DF(2,2,J)
DFINV(2,2,J)=DF(1,1,J)
DFINV(1,2,J)=-DF(2,1,J)
DFINV(2,1,J)=-DF(1,2,J)
DFINV(:,:,J)=DFINV(:,:,J)/(DF(1,1,J)*DF(2,2,J)-DF(1,2,J)*DF(2,1,J))
END DO
 DO J=1,NPhi
     DO K=1,NPhi
        M_courbe(J,K)=sum(Phi2D(J,:)*Phi2D(K,:)*wGL2D*DFVEC_BIS)
     END DO
  END DO
    DO J=1,Nphi
       DO K=1,Nphi
          error_p = error_p + (P_analytic((I-1)*Nphi+J)-P_complexe((I-1)*Nphi+J))*&
               &conjg(P_analytic((I-1)*Nphi+K)-P_complexe((I-1)*Nphi+K))*M_courbe(J,K)
       ENDDO
    End Do
 end Do

	Do I=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu
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
IF(Neigh(I,Edge).gt.Nflu+Nflusol) then
   DO K=2,order
indice=indice+1
!!!      x_k=x_1+(k-1)*(x_2-x_1)/order (x1 and x2 are the two vertices)
      pt_courbe_x(indice)= pt_courbe_x(J1)+(k-1)&
           &*(pt_courbe_x(J2)-pt_courbe_x(J1))/real(order,8)
      pt_courbe_y(indice)=pt_courbe_y(J1)+(k-1)&
           &*(pt_courbe_y(J2)-pt_courbe_y(J1))/real(order,8)
   ENDDO
else
call calc_theta(pt_courbe_x(J1), pt_courbe_y(J1),theta1)
call calc_theta(pt_courbe_x(J2), pt_courbe_y(J2),theta2)
dtheta=theta2-theta1
if (dtheta.gt.pi) then
dtheta=dtheta-2*pi
elseif (dtheta.lt.-pi) then
dtheta=dtheta+2*pi
end if
   DO K=2,order
indice=indice+1
      pt_courbe_x(indice)= R0*dcos(theta1+(k-1)*dtheta/real(order,8))
      pt_courbe_y(indice)= R0*dsin(theta1+(k-1)*dtheta/real(order,8))
   ENDDO
end IF
ENDDO
IF (ORDER.EQ.3) then
   indice=10
   pt_courbe_x(indice)=sum(pt_courbe_x(1:3))/3D0
   pt_courbe_y(indice)=sum(pt_courbe_y(1:3))/3D0
ENDIF
DF=0.D0
DO J=1,NGL2D
DF(1,1,J)=sum(GradPhi2D(:,1,J)*pt_courbe_x(:))
DF(1,2,J)=sum(GradPhi2D(:,2,J)*pt_courbe_x(:))
DF(2,1,J)=sum(GradPhi2D(:,1,J)*pt_courbe_y(:))
DF(2,2,J)=sum(GradPhi2D(:,2,J)*pt_courbe_y(:))
DFVEC_BIS(J)=abs(DF(1,1,J)*DF(2,2,J)-DF(1,2,J)*DF(2,1,J))
!!! Attention, DFINV est la transposée de JFinv
DFINV(1,1,J)=DF(2,2,J)
DFINV(2,2,J)=DF(1,1,J)
DFINV(1,2,J)=-DF(2,1,J)
DFINV(2,1,J)=-DF(1,2,J)
DFINV(:,:,J)=DFINV(:,:,J)/(DF(1,1,J)*DF(2,2,J)-DF(1,2,J)*DF(2,1,J))
END DO
 DO J=1,NPhi
     DO K=1,NPhi
        M_courbe(J,K)=sum(Phi2D(J,:)*Phi2D(K,:)*wGL2D*DFVEC_BIS)
     END DO
  END DO

  	L=((I-1)-Nflu-Nflusol)*Nphi
   		DO J=1,Nphi
			DO K=1,Nphi
				error_ux = error_ux + (Ux_analytic(L+J)-Ux_complexe(L+J))*&
				&conjg(Ux_analytic(L+K)-Ux_complexe(L+K))*M_courbe(J,K)
				error_uy = error_uy + (Uy_analytic(L+J)-Uy_complexe(L+J))*&
				&conjg(Uy_analytic(L+K)-Uy_complexe(L+K))*M_courbe(J,K)
			END DO
		ENDDO
End Do


	Do I=Nflu+Nflusol+Nsolflu+1,Nflu+Nflusol+Nsolflu+Nsol
  	L=((I-1)-Nflu-Nflusol)*Nphi
   		DO J=1,Nphi
			DO K=1,Nphi
				error_ux = error_ux + (Ux_analytic(L+J)-Ux_complexe(L+J))*&
				&conjg(Ux_analytic(L+K)-Ux_complexe(L+K))*DFVEC(I)*M(J,K)
				error_uy = error_uy + (Uy_analytic(L+J)-Uy_complexe(L+J))*&
				&conjg(Uy_analytic(L+K)-Uy_complexe(L+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
End Do


	error_p=sqrt(error_p)
        error_ux=sqrt(error_ux)
        error_uy=sqrt(error_uy)

      end subroutine sub_error_abs_scalaire_vectoriel_courbe

subroutine sub_error_rel_scalaire_vectoriel_courbe(error_p,error_ux,error_uy,R0)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K,L,J1,J2,indice,edge,node(3)
	real*8 :: error_p,error_ux,error_uy,error_p_num,error_p_den,error_ux_num,error_ux_den,error_uy_num,error_uy_den
	real*8::pi,d,dtheta,pt_courbe_x(Nphi),pt_courbe_y(Nphi)
 REAL*8,allocatable::Phi2D(:,:),M_courbe(:,:)
  REAL*8,allocatable :: GradPhi2D(:,:,:)
  REAL*8 :: DF(2,2,NGL2D),DFVEC_bis(NGL2D),theta1,theta2,R0,vectmp1(2),vectmp2(2),DFINV(2,2,NGL2D)
 pi=2.d0*dasin(1.d0)
  allocate(GradPhi2D(Nphi,2,NGL2D))
  allocate(Phi2D(Nphi,NGL2D))
  allocate(M_courbe(Nphi,Nphi))
  SELECT CASE(Order)
  CASE(1)

     CALL GradPhiOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     
 

  CASE(2)
     CALL GradPhiOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     

  CASE(3)
     CALL GradPhiOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,GradPhi2D)     
     CALL Phi2DOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi2D)     
  END SELECT	
	!! P_complexe est la solution exacte et P_analytic la solution approchée
	call sub_error_abs_scalaire_vectoriel_courbe(error_p_num,error_ux_num,error_uy_num,R0)
	error_p_den=0.d0
	error_ux_den=0.d0
	error_uy_den=0.d0
	Do I=1,Nflu 
		DO J=1,Nphi
			DO K=1,Nphi
				error_p_den = error_p_den + (P_analytic((I-1)*Nphi+J))*&
				&conjg(P_analytic((I-1)*Nphi+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
	End Do

	Do I=Nflu+1,Nflu+Nflusol 
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
IF(Neigh(I,Edge).le.Nflu+Nflusol) then
   DO K=2,order
indice=indice+1
!!!      x_k=x_1+(k-1)*(x_2-x_1)/order (x1 and x2 are the two vertices)
      pt_courbe_x(indice)= pt_courbe_x(J1)+(k-1)&
           &*(pt_courbe_x(J2)-pt_courbe_x(J1))/real(order,8)
      pt_courbe_y(indice)=pt_courbe_y(J1)+(k-1)&
           &*(pt_courbe_y(J2)-pt_courbe_y(J1))/real(order,8)
   ENDDO
else
call calc_theta(pt_courbe_x(J1), pt_courbe_y(J1),theta1)
call calc_theta(pt_courbe_x(J2), pt_courbe_y(J2),theta2)
dtheta=theta2-theta1
if (dtheta.gt.pi) then
dtheta=dtheta-2*pi
elseif (dtheta.lt.-pi) then
dtheta=dtheta+2*pi
end if
   DO K=2,order
indice=indice+1
      pt_courbe_x(indice)= R0*dcos(theta1+(k-1)*dtheta/real(order,8))
      pt_courbe_y(indice)= R0*dsin(theta1+(k-1)*dtheta/real(order,8))
   ENDDO
end IF
ENDDO
IF (ORDER.EQ.3) then
   indice=10
   pt_courbe_x(indice)=sum(pt_courbe_x(1:3))/3D0
   pt_courbe_y(indice)=sum(pt_courbe_y(1:3))/3D0
ENDIF
DF=0.D0
DO J=1,NGL2D
DF(1,1,J)=sum(GradPhi2D(:,1,J)*pt_courbe_x(:))
DF(1,2,J)=sum(GradPhi2D(:,2,J)*pt_courbe_x(:))
DF(2,1,J)=sum(GradPhi2D(:,1,J)*pt_courbe_y(:))
DF(2,2,J)=sum(GradPhi2D(:,2,J)*pt_courbe_y(:))
DFVEC_BIS(J)=abs(DF(1,1,J)*DF(2,2,J)-DF(1,2,J)*DF(2,1,J))
!!! Attention, DFINV est la transposée de JFinv
DFINV(1,1,J)=DF(2,2,J)
DFINV(2,2,J)=DF(1,1,J)
DFINV(1,2,J)=-DF(2,1,J)
DFINV(2,1,J)=-DF(1,2,J)
DFINV(:,:,J)=DFINV(:,:,J)/(DF(1,1,J)*DF(2,2,J)-DF(1,2,J)*DF(2,1,J))
END DO
 DO J=1,NPhi
     DO K=1,NPhi
        M_courbe(J,K)=sum(Phi2D(J,:)*Phi2D(K,:)*wGL2D*DFVEC_BIS)
     END DO
  END DO
		DO J=1,Nphi
			DO K=1,Nphi
				error_p_den = error_p_den + (P_analytic((I-1)*Nphi+J))*&
				&conjg(P_analytic((I-1)*Nphi+K))*M_courbe(J,K)
			END DO
		ENDDO
	End Do

	error_p=error_p_num/sqrt(error_p_den)*100

	Do I=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu
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
IF(Neigh(I,Edge).gt.Nflu+Nflusol) then
   DO K=2,order
indice=indice+1
!!!      x_k=x_1+(k-1)*(x_2-x_1)/order (x1 and x2 are the two vertices)
      pt_courbe_x(indice)= pt_courbe_x(J1)+(k-1)&
           &*(pt_courbe_x(J2)-pt_courbe_x(J1))/real(order,8)
      pt_courbe_y(indice)=pt_courbe_y(J1)+(k-1)&
           &*(pt_courbe_y(J2)-pt_courbe_y(J1))/real(order,8)
   ENDDO
else
call calc_theta(pt_courbe_x(J1), pt_courbe_y(J1),theta1)
call calc_theta(pt_courbe_x(J2), pt_courbe_y(J2),theta2)
dtheta=theta2-theta1
if (dtheta.gt.pi) then
dtheta=dtheta-2*pi
elseif (dtheta.lt.-pi) then
dtheta=dtheta+2*pi
end if
   DO K=2,order
indice=indice+1
      pt_courbe_x(indice)= R0*dcos(theta1+(k-1)*dtheta/real(order,8))
      pt_courbe_y(indice)= R0*dsin(theta1+(k-1)*dtheta/real(order,8))
   ENDDO
ENDIF
ENDDO
IF (ORDER.EQ.3) then
   indice=10
   pt_courbe_x(indice)=sum(pt_courbe_x(1:3))/3D0
   pt_courbe_y(indice)=sum(pt_courbe_y(1:3))/3D0
ENDIF
DF=0.D0
DO J=1,NGL2D
DF(1,1,J)=sum(GradPhi2D(:,1,J)*pt_courbe_x(:))
DF(1,2,J)=sum(GradPhi2D(:,2,J)*pt_courbe_x(:))
DF(2,1,J)=sum(GradPhi2D(:,1,J)*pt_courbe_y(:))
DF(2,2,J)=sum(GradPhi2D(:,2,J)*pt_courbe_y(:))
DFVEC_BIS(J)=abs(DF(1,1,J)*DF(2,2,J)-DF(1,2,J)*DF(2,1,J))
!!! Attention, DFINV est la transposée de JFinv
DFINV(1,1,J)=DF(2,2,J)
DFINV(2,2,J)=DF(1,1,J)
DFINV(1,2,J)=-DF(2,1,J)
DFINV(2,1,J)=-DF(1,2,J)
DFINV(:,:,J)=DFINV(:,:,J)/(DF(1,1,J)*DF(2,2,J)-DF(1,2,J)*DF(2,1,J))
END DO
 DO J=1,NPhi
     DO K=1,NPhi
        M_courbe(J,K)=sum(Phi2D(J,:)*Phi2D(K,:)*wGL2D*DFVEC_BIS)
     END DO
  END DO

  	L=((I-1)-Nflu-Nflusol)*Nphi

		DO J=1,Nphi
			DO K=1,Nphi

				error_ux_den = error_ux_den + (Ux_analytic(L+J))*&
				&conjg(Ux_analytic(L+K))*M_courbe(J,K)
				error_uy_den = error_uy_den + (Uy_analytic(L+J))*&
				&conjg(Uy_analytic(L+K))*M_courbe(J,K)
			END DO
		ENDDO
	End Do        

	Do I=Nflu+Nflusol+Nsolflu+1,Nflu+Nflusol+Nsolflu+Nsol
  	L=((I-1)-Nflu-Nflusol)*Nphi

		DO J=1,Nphi
			DO K=1,Nphi
				error_ux_den = error_ux_den + (Ux_analytic(L+J))*&
				&conjg(Ux_analytic(L+K))*DFVEC(I)*M(J,K)
				error_uy_den = error_uy_den + (Uy_analytic(L+J))*&
				&conjg(Uy_analytic(L+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
	End Do        

 error_ux=error_ux_num/sqrt(error_ux_den)*100
 error_uy=error_uy_num/sqrt(error_uy_den)*100	
 

end subroutine sub_error_rel_scalaire_vectoriel_courbe
