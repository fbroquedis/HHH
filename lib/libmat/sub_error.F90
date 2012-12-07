subroutine sub_error_abs(error)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K
	real*8 :: error
	

	!! Uold est la solution exacte et U la solution approchée
	error=0.d0
	Do I=1,Ntri
	
     
		DO J=1,Nphi
			DO K=1,Nphi
				error = error + (P_analytic((I-1)*Nphi+J)-P_complexe((I-1)*Nphi+J))*&
				&conjg(P_analytic((I-1)*Nphi+K)-P_complexe((I-1)*Nphi+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
	End Do
        error=sqrt(error)
end subroutine
subroutine sub_error_abs2(error)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K
	real*8 :: error
	

	!! P_complexe est la solution exacte et P_analytic la solution approchée
	error=0.d0
	Do I=1,Ntri
	
     
		DO J=1,Nphi
			DO K=1,Nphi
				error = error + (P_analytic((I-1)*Nphi+J)-P_complexe((I-1)*Nphi+J))*&
				&(P_analytic((I-1)*Nphi+K)-P_complexe((I-1)*Nphi+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
	End Do
        error=sqrt(error)
end subroutine

subroutine sub_error_rel(error)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K
	real*8 :: error,error_num,error_den
	
	!! P_complexe est la solution exacte et P_analytic la solution approchée
	call sub_error_abs(error_num)
	error_den=0.d0
	Do I=1,Ntri     
		DO J=1,Nphi
			DO K=1,Nphi
				error_den = error_den + (P_analytic((I-1)*Nphi+J))*&
				&conjg(P_analytic((I-1)*Nphi+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
	End Do
	error=error_num/sqrt(error_den)*100
	
end subroutine
subroutine sub_error_rel2(error)
	Use m_mat
	Use m_gen
	Use m_mesh
	implicit none
	
	INTEGER :: I,J,K
	real*8 :: error,error_num,error_den
	
	!! P_complexe est la solution exacte et P_analytic la solution approchée
	call sub_error_abs2(error_num)
	error_den=0.d0
	Do I=1,Ntri     
		DO J=1,Nphi
			DO K=1,Nphi
				error_den = error_den + (P_analytic((I-1)*Nphi+J))*&
				&(P_analytic((I-1)*Nphi+K))*DFVEC(I)*M(J,K)
			END DO
		ENDDO
	End Do

	error=(1.D0-error_num/sqrt(error_den))*100
	
end subroutine
