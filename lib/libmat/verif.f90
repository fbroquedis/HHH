SUBROUTINE verif
  Use m_mat
  Use m_gen
  Use m_mesh
  implicit none
  INTEGER :: I,INeigh,J,K,L1,L2,Face,Face2,TEST,L
  INTEGER,allocatable :: PhiEdge(:,:)
  allocate(PhiEdge(3,1+Order))
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

  END SELECT

     DO I=1,Ntri
        DO J=1,Nphi
           DO K=1,Nphi
              IF(A(I,J,K).ne.A(I,K,J)) then
                 write(6,*) 'error'
                 write(6,*) I,J,K
                 write(6,*) A(I,J,K)
                 write(6,*) A(I,K,J)
                 stop
              END IF
           END DO
        END DO
        L1=Nphi
     Do Face=1,3
        INeigh=Neigh(I,Face)
        IF(Ineigh.ne.-1) then
           L2=Nphi
        DO K=1,3
            IF (Neigh(INeigh,K).eq.I) THEN
               Face2=K
              EXIT 
           ELSE
              IF (Neigh(INeigh,K).ne.-1) THEN
                 L2=L2+Nphi
              ENDIF
           ENDIF
        ENDDO
        DO J=1,Nphi
           DO K=1,Nphi
              Test=1
              DO L=1,1+Order
                 IF((J.eq.PhiEdge(Face,L)).or.(K.eq.PhiEdge(Face2,L))) then
                    test=0
                 ENDIF
              ENDDO
              If (Test.eq.1) then
                 If(A(I,J,L1+K).ne.0) then
                    write(6,*) 'Error',A(I,j,L1+K)
                    write(6,*) I,J,K,L1
                 ENDIF
              ENDIF
              IF(abs(A(I,J,L1+K)-A(INeigh,K,L2+J)).ge.1E6*spacing(A(I,J,L1+K))) then
                 write(6,*) 'errorr'
                 write(6,*) I,Ineigh,J,K,L1,L2
                 write(6,*) A(I,J,L1+K)
                 write(6,*) A(INeigh,K,L2+J)
                 stop
              END IF
              
           END DO
        END DO
     L1=L1+Nphi
  END IF
  END DO
END DO
END SUBROUTINE verif
