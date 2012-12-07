SUBROUTINE sub_refine
  Use m_mesh
  implicit none
  INTEGER :: I,J,K,Node(6),Test
  INTEGER :: Iraff,Inoraff,Iboundraff,Iboundnoraff,Ineigh
  INTEGER,allocatable :: Testraff(:),Testraffinter(:),Corres(:),TestK(:)
  integer,dimension(:,:),allocatable,save :: TriInter,NeighInter
  real*8,dimension (:,:),allocatable,save :: CoorInter
  real*8 :: V12(2),V13(2),V23(2),x,y,h1
  h=0.D0
  hraff=0.D0
  Nraff=0
  allocate(Testraff(Ntri))
  DO I=1,Ntri
     Node(1:3) = Tri(I,:)
     V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
     V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
     V13(1)=Coor(Node(3),1)-Coor(Node(1),1)
     V13(2)=Coor(Node(3),2)-Coor(Node(1),2)
     V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
     V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
     h1=sqrt(sum(V12**2))
     h1=min(h1,sqrt(sum(V23**2)))
     h1=min(h1,sqrt(sum(V13**2)))
     x=Coor(Node(1),1)+1/3.D0*V12(1)+1/3.D0*V13(1)
     y=Coor(Node(1),2)+1/3.D0*V12(2)+1/3.D0*V13(2)
     IF (((x-x1reff)*(x-x2reff).le.0.d0).and.((y-y1reff)*(y-y2reff).le.0.d0)) then
        Testraff(I)=3
        Nraff=Nraff+1
        IF(hraff.eq.0) then 
           hraff=h1/2.D0
        ELSE
           hraff=min(hraff,h1/2.D0)
        ENDIF
     ELSE
        Testraff(I)=1
       IF(h.eq.0) then 
           h=h1
        ELSE
           h=min(h,h1)
        ENDIF
     END IF
  ENDDO
  write(6,*) Nraff
  allocate(TriInter(Ntri+3*Nraff,3))
  allocate(NeighInter(Ntri+3*Nraff,6))
  allocate(CoorInter(NPoints+3*Nraff,2))
  allocate(TestraffInter(Ntri))
  TestraffInter=Testraff
  deallocate(TestRaff)
  allocate(Testraff(Ntri+3*Nraff))
  Testraff(1:Ntri)=TestraffInter
  deallocate(TestRaffInter)
  Iraff=0
  Inoraff=0
  Iboundraff=0
  Iboundnoraff=0
  TriInter=0
  TriInter(1:Ntri,:)=Tri
  NeighInter=0
  NeighInter(1:Ntri,1:3)=Neigh
  CoorInter=0.D0
  CoorInter(1:NPoints,:)=Coor
  DO I=1,Ntri
     IF(TestRaff(I).le.2) then
        test=0
        DO J=1,3
           IF (Neigh(I,J).ne.-1) then
              IF (TestRaff(Neigh(I,J)).ge.3) then
                 test=1
              END IF
           END IF
        END DO
        IF(test.eq.1) then
           Iboundnoraff=Iboundnoraff+1
           TestRaff(I)=2
        ELSE
           Inoraff=Inoraff+1
        ENDIF
     ELSE
        TestRaff(Ntri+1:Ntri+3)=3
        Node(1:3) = Tri(I,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V13(1)=Coor(Node(3),1)-Coor(Node(1),1)
        V13(2)=Coor(Node(3),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        test=0
!        Edge 1 (23)
        x=Coor(Node(2),1)+1/2.D0*V23(1)
        y=Coor(Node(2),2)+1/2.D0*V23(2)
        Ineigh=Neigh(I,1)
           IF (Ineigh.ne.-1) then
              IF (TestRaff(Ineigh).le.2) then
                 test=test+1
                 Npoints=Npoints+1
                 CoorInter(Npoints,1)=x
                 CoorInter(Npoints,2)=y
                 Node(4)=Npoints
                 NeighInter(Ntri+2,1)=Ineigh
                 NeighInter(Ntri+3,1)=Ineigh
                 NeighInter(Ntri+2,4)=2
                 NeighInter(Ntri+3,4)=1
                 TestRaff(Ntri+2)=4
                 TestRaff(Ntri+3)=4
                 DO J=1,3
                    IF (Neigh(INeigh,J).eq.I) then
                       NeighInter(INeigh,J)=Ntri+3
                       NeighInter(INeigh,J+3)=Ntri+2
                       EXIT
                    END IF
                 END DO
              ELSE
                 IF(I.lt.Ineigh) then
                    Npoints=Npoints+1
                    CoorInter(Npoints,1)=x
                    CoorInter(Npoints,2)=y
                    Node(4)=Npoints                    
                    DO J=1,3
                       IF (Neigh(INeigh,J).eq.I) then
                          NeighInter(INeigh,J)=Ntri+3
                          NeighInter(INeigh,J+3)=Ntri+2
                          EXIT
                       END IF
                    END DO
                 ELSE
                    NeighInter(Ntri+2,1)=NeighInter(I,1)
                    NeighInter(Ntri+3,1)=NeighInter(I,4)
                    DO J=1,3
                       IF (Neigh(INeigh,J).eq.I) then
                          NeighInter(NeighInter(I,1),J)=Ntri+2
                          NeighInter(NeighInter(I,4),J)=Ntri+3
                          SELECT CASE(J)
                          CASE(1)
                             Node(4)=TriInter(NeighInter(Ntri+2,1),2)
                             IF (Node(4).ne.TriInter(NeighInter(Ntri+3,1),3)) then
                                write(6,*) 'error1'
                                stop
                             END IF
                          CASE(2)
                             Node(4)=TriInter(NeighInter(Ntri+2,1),3)
                             IF (Node(4).ne.TriInter(NeighInter(Ntri+3,1),1)) then
                                write(6,*) 'error4'
                                stop
                             END IF
                          CASE(3)
                             Node(4)=TriInter(NeighInter(Ntri+2,1),1)
                             IF (Node(4).ne.TriInter(NeighInter(Ntri+3,1),2)) then
                                write(6,*) 'error7'
                                stop
                             END IF
                          END SELECT
                          EXIT
                       END IF
                    END DO
                 END IF
              END IF
           ELSE
              Npoints=Npoints+1
              CoorInter(Npoints,1)=x
              CoorInter(Npoints,2)=y
              Node(4)=Npoints
              NeighInter(Ntri+2,1)=-1
              NeighInter(Ntri+3,1)=-1
           END IF
!        Edge 2 (31)
        Ineigh=Neigh(I,2)
        x=Coor(Node(1),1)+1/2.D0*V13(1)
        y=Coor(Node(1),2)+1/2.D0*V13(2)
           IF (Ineigh.ne.-1) then
              IF (TestRaff(Neigh(I,2)).le.2) then
                 test=test+1
                 Npoints=Npoints+1
                 CoorInter(Npoints,1)=x
                 CoorInter(Npoints,2)=y
                 Node(5)=Npoints
                 NeighInter(Ntri+1,2)=Ineigh
                 NeighInter(Ntri+3,2)=Ineigh
                 NeighInter(Ntri+1,5)=1
                 NeighInter(Ntri+3,5)=2
                 TestRaff(Ntri+1)=4
                 TestRaff(Ntri+3)=4
                 DO J=1,3
                    IF (Neigh(INeigh,J).eq.I) then
                       NeighInter(INeigh,J)=Ntri+1
                       NeighInter(INeigh,J+3)=Ntri+3
                       EXIT
                    END IF
                 END DO
              ELSE
                 IF(I.lt.Ineigh) then
                 Npoints=Npoints+1
                 CoorInter(Npoints,1)=x
                 CoorInter(Npoints,2)=y
                 Node(5)=Npoints                    
                 DO J=1,3
                    IF (Neigh(INeigh,J).eq.I) then
                       NeighInter(INeigh,J)=Ntri+1
                       NeighInter(INeigh,J+3)=Ntri+3
                    EXIT
                    END IF
                 END DO
                 ELSE
                    NeighInter(Ntri+3,2)=NeighInter(I,2)
                    NeighInter(Ntri+1,2)=NeighInter(I,5)
                    DO J=1,3
                       IF (Neigh(INeigh,J).eq.I) then
                          NeighInter(NeighInter(I,2),J)=Ntri+3
                          NeighInter(NeighInter(I,5),J)=Ntri+1
                          SELECT CASE(J)
                          CASE(1)
                             Node(5)=TriInter(NeighInter(Ntri+3,2),2)
                             IF (Node(5).ne.TriInter(NeighInter(Ntri+1,2),3)) then
                                write(6,*) 'error11'
                                stop
                             END IF
                             IF (CoorInter(Node(5),1).ne.x) then
                                write(6,*) 'error12'
                                stop
                             END IF
                             IF (CoorInter(Node(5),2).ne.y) then
                                write(6,*) 'error13'
                                stop
                             END IF
                          CASE(2)
                             Node(5)=TriInter(NeighInter(Ntri+3,2),3)
                             IF (Node(5).ne.TriInter(NeighInter(Ntri+1,2),1)) then
                                write(6,*) 'error14'
                                stop
                             END IF
                          CASE(3)
                             Node(5)=TriInter(NeighInter(Ntri+3,2),1)
                             IF (Node(5).ne.TriInter(NeighInter(Ntri+1,2),2)) then
                                write(6,*) 'error17',NeighInter(Ntri+3,2)
                                stop
                             END IF
                             IF (CoorInter(Node(5),1).ne.x) then
                                write(6,*) 'error18'
                                stop
                             END IF
                             IF (CoorInter(Node(5),2).ne.y) then
                                write(6,*) 'error19'
                                stop
                             END IF               
                          END SELECT
                          EXIT
                       END IF
                    END DO
                 END IF                 
              END IF
           ELSE
              Npoints=Npoints+1
              CoorInter(Npoints,1)=x
              CoorInter(Npoints,2)=y
              Node(5)=Npoints
              NeighInter(Ntri+1,2)=-1
              NeighInter(Ntri+3,2)=-1
           END IF
!        Edge 3 (12)
           x=Coor(Node(1),1)+1/2.D0*V12(1)
           y=Coor(Node(1),2)+1/2.D0*V12(2)
           Ineigh=Neigh(I,3)
           IF (Ineigh.ne.-1) then
              IF (TestRaff(Ineigh).le.2) then
                 test=test+1
                 Npoints=Npoints+1
                 CoorInter(Npoints,1)=x
                 CoorInter(Npoints,2)=y
                 Node(6)=Npoints
                 NeighInter(Ntri+1,3)=Ineigh
                 NeighInter(Ntri+2,3)=Ineigh
                 NeighInter(Ntri+1,6)=2
                 NeighInter(Ntri+2,6)=1
                 TestRaff(Ntri+1)=4
                 TestRaff(Ntri+2)=4
                 DO J=1,3
                    IF (Neigh(INeigh,J).eq.I) then
                       NeighInter(INeigh,J)=Ntri+2
                       NeighInter(INeigh,J+3)=Ntri+1
                       EXIT
                    END IF
                 END DO
              ELSE
                 IF(I.lt.Ineigh) then
                    Npoints=Npoints+1
                    CoorInter(Npoints,1)=x
                    CoorInter(Npoints,2)=y
                    Node(6)=Npoints                    
                    DO J=1,3
                       IF (Neigh(INeigh,J).eq.I) then
                          NeighInter(INeigh,J)=Ntri+2
                          NeighInter(INeigh,J+3)=Ntri+1
                          EXIT
                       END IF
                    END DO
                 ELSE
                    NeighInter(Ntri+1,3)=NeighInter(I,3)
                    NeighInter(Ntri+2,3)=NeighInter(I,6)
                 DO J=1,3
                       IF (Neigh(INeigh,J).eq.I) then
                          NeighInter(NeighInter(I,3),J)=Ntri+1
                          NeighInter(NeighInter(I,6),J)=Ntri+2
                          SELECT CASE(J)
                          CASE(1)
                             Node(6)=TriInter(NeighInter(Ntri+1,3),2)
                          CASE(2)
                             Node(6)=TriInter(NeighInter(Ntri+1,3),3)
                             IF (Node(6).ne.TriInter(NeighInter(Ntri+2,3),1)) then
                                write(6,*) 'error24'
                                stop
                             END IF
                             IF (CoorInter(Node(6),1).ne.x) then
                                write(6,*) 'error25'
                                stop
                             END IF
                             IF (CoorInter(Node(6),2).ne.y) then
                                write(6,*) 'error26'
                                stop
                             END IF
                          CASE(3)
                             Node(6)=TriInter(NeighInter(Ntri+1,3),1)
                             IF (Node(6).ne.TriInter(NeighInter(Ntri+2,3),2)) then
                                write(6,*) 'error27'
                                stop
                             END IF
                          END SELECT
                          EXIT
                       END IF
                    END DO
                 END IF
              END IF
           ELSE
              Npoints=Npoints+1
              CoorInter(Npoints,1)=x
              CoorInter(Npoints,2)=y
              Node(6)=Npoints
              NeighInter(Ntri+1,3)=-1
              NeighInter(Ntri+2,3)=-1
           END IF
        IF(test.eq.1) then
           Iboundraff=Iboundraff+2
           Iraff=Iraff+2
        ELSEIF (test.ge.2) then
           Iboundraff=Iboundraff+3
           Iraff=Iraff+1
        ELSE
           Iraff=Iraff+4
        ENDIF
        TriInter(I,1)=Node(4)
        TriInter(I,2)=Node(5)
        TriInter(I,3)=Node(6)
        NeighInter(I,1)=Ntri+1
        NeighInter(I,2)=Ntri+2
        NeighInter(I,3)=Ntri+3
        NeighInter(I,4:6)=0
        TriInter(Ntri+1,1)=Node(1)
        TriInter(Ntri+1,2)=Node(6)
        TriInter(Ntri+1,3)=Node(5)
        NeighInter(Ntri+1,1)=I
        TriInter(Ntri+2,1)=Node(6)
        TriInter(Ntri+2,2)=Node(2)
        TriInter(Ntri+2,3)=Node(4)
        NeighInter(Ntri+2,2)=I
        TriInter(Ntri+3,1)=Node(5)
        TriInter(Ntri+3,2)=Node(4)
        TriInter(Ntri+3,3)=Node(3)
        NeighInter(Ntri+3,3)=I
        Ntri=Ntri+3
     END IF
  END DO
deallocate(Tri)
deallocate(Neigh)
deallocate(Coor)
allocate(Tri(Ntri,3))
allocate(Neigh(Ntri,6))
allocate(Coor(NPoints,2))
allocate(Corres(Ntri))
allocate(TestK(Ntri))
Tri=0
Neigh=0
Coor=Coorinter(1:Npoints,:)
deallocate(CoorInter)
!!!! renumerotation
Nraff=Iraff
Nnoraff=Inoraff
Nboundraff=Iboundraff
Nboundnoraff=Iboundnoraff
IF(Nraff+Nnoraff+Nboundraff+Nboundnoraff.ne.Ntri) then
   write(6,*)'error',Nraff+Nnoraff+Nboundraff+Nboundnoraff,Ntri
ENDIF

Iraff=0
Inoraff=0
Iboundraff=0
Iboundnoraff=0
Corres=0
TestK=0
DO I=1,Ntri
   SELECT CASE(TestRaff(I))
   CASE(1)
      Inoraff=Inoraff+1
      K=Inoraff
   CASE(2)
      Iboundnoraff=Iboundnoraff+1
      K=Nnoraff+Iboundnoraff
   CASE(3)
      Iraff=Iraff+1
      K=Nnoraff+Nboundnoraff+Nboundraff+Iraff
   CASE(4)
      Iboundraff=Iboundraff+1
      K=Nnoraff+Nboundnoraff+Iboundraff      
   END SELECT
   IF (TestK(K).ne.0)then
      write(6,*) 'error',TestK(K),K,I,TestRaff(I)
      stop
   ENDIF
   IF(K.gt.Ntri) then
      write(6,*) 'error',K,Ntri
      stop
ENDIF
      Corres(I)=K
      TestK(K)=I
   Tri(K,:)=TriInter(I,:)
ENDDO
DO I=1,NTri
   K=Corres(I)
   IF ((K.le.Nnoraff).or.(k.gt.Nnoraff+Nboundnoraff+Nboundraff)) then
      DO J=1,3
         Ineigh=NeighInter(I,J)
         IF(Ineigh.gt.0) then
            Neigh(K,J)=Corres(Ineigh)
         ELSE
            Neigh(K,J)=Ineigh
         ENDIF
      ENDDO
      Neigh(K,4:6)=0
   ENDIF
   IF ((K.gt.Nnoraff).and.(K.le.Nnoraff+Nboundnoraff)) then
      K=Corres(I)
      DO J=1,6
         Ineigh=NeighInter(I,J)
         IF(Ineigh.gt.0) then
            Neigh(K,J)=Corres(Ineigh)
         ELSE
            Neigh(K,J)=Ineigh
         ENDIF
      ENDDO
   ENDIF
   IF ((K.gt.Nnoraff+Nboundnoraff).and.(K.le.Nnoraff+Nboundnoraff+Nboundraff)) then
      K=Corres(I)
      DO J=1,3
         Ineigh=NeighInter(I,J)
         IF(Ineigh.gt.0) then
            Neigh(K,J)=Corres(Ineigh)
         ELSE
            Neigh(K,J)=Ineigh
         ENDIF
      ENDDO
      DO J=4,6
         Ineigh=NeighInter(I,J)
         Neigh(K,J)=Ineigh
      ENDDO
   ENDIF
ENDDO
deallocate(TriInter)
deallocate(NeighInter)
END SUBROUTINE sub_refine
