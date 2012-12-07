subroutine sub_error_abs_scalaire_vectoriel(error_p, error_ux, error_uy)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K,L
	real*8 :: error_p, error_ux, error_uy
	
	!! Uold est la solution exacte et U la solution approchée
	error_p=0.d0
	error_ux=0.d0
	error_uy=0.d0
	Do I=1,Nflu+Nflusol

   		DO J=1,Nphi
			DO K=1,Nphi
				error_p = error_p + (P_analytic((I-1)*Nphi+J)-P_complexe((I-1)*Nphi+J))*&
				&conjg(P_analytic((I-1)*Nphi+K)-P_complexe((I-1)*Nphi+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
	End Do

	Do I=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu+Nsol
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

end subroutine sub_error_abs_scalaire_vectoriel

subroutine sub_error_rel_scalaire_vectoriel(error_p,error_ux,error_uy)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K,L
	real*8 :: error_p,error_ux,error_uy,error_p_num,error_p_den,error_ux_num,error_ux_den,error_uy_num,error_uy_den
	
	!! P_complexe est la solution exacte et P_analytic la solution approchée
	call sub_error_abs_scalaire_vectoriel(error_p_num,error_ux_num,error_uy_num)
	error_p_den=0.d0
	error_ux_den=0.d0
	error_uy_den=0.d0
	Do I=1,Nflu+Nflusol  
		DO J=1,Nphi
			DO K=1,Nphi
				error_p_den = error_p_den + (P_analytic((I-1)*Nphi+J))*&
				&conjg(P_analytic((I-1)*Nphi+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
	End Do
	error_p=error_p_num/sqrt(error_p_den)*100

	Do I=Nflu+Nflusol+1,Nflu+Nflusol+Nsolflu+Nsol
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
 

end subroutine sub_error_rel_scalaire_vectoriel
