SUBROUTINE sub_defmassmat
  Use m_mat
  Use m_gen
  implicit none
  INTEGER :: I,J,Lwork,INFO
  REAL*8 :: Phi(Nphi,NGL2D),Diag(Nphi,Nphi),Eig(Nphi),Work(340)
  allocate(M(Nphi,Nphi))
  SELECT CASE(ORDER)
  CASE(1)
     CALL Phi2DOrder1(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi)     
  CASE(2)
     CALL Phi2DOrder2(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi)     
  CASE(3)
     CALL Phi2DOrder3(PtGL2D(:,1),PtGL2D(:,2),NGL2D,Phi)     
  END SELECT
  DO I=1,NPhi
     DO J=1,NPhi
        M(I,J)=sum(Phi(I,:)*Phi(J,:)*wGL2D)
     END DO
  END DO
  allocate(Minv(Nphi,Nphi))
  Minv=M
  CALL DPOTRF('U',Nphi,Minv,Nphi,INFO)
  CALL DPOTRI('U',Nphi,Minv,Nphi,INFO)
  DO I=1,Nphi
     DO J=I+1,Nphi
       Minv(J,I)=Minv(I,J)
     end DO
  end DO
END SUBROUTINE sub_defmassmat
