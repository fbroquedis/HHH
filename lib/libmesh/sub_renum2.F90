SUBROUTINE sub_renum2
  Use m_mesh
  implicit none
  INTEGER :: I,J,K,Test
  INTEGER :: Iraff,Inoraff,Iboundraff,Iboundnoraff,Ineigh
  INTEGER,allocatable :: Testraffinter(:),Corres(:),TestK(:)
  INTEGER,dimension(:),allocatable,save :: Testraff
  INTEGER,dimension(:,:),allocatable,save :: Flusol
  integer,dimension(:,:),allocatable,save :: TriInter,NeighInter
  Nflu=0
  Nsol=0     
  Nflusol=0
  Nsolflu=0
  allocate(flusol(2,Ntri))

  DO I=1,Ntri
      if (TYPE_MEDIA(I).eq.1) then
         test=0      
         DO J=1,3
            IF (Neigh(J,I).gt.0) then
               if (Type_Media(Neigh(J,I)).ne.1) then
                  test=1
               end if
            end IF
         END DO
         if (test.eq.0) then 
            Nflu=Nflu+1
            Flusol(1,I)=1
            Flusol(2,I)=Nflu
         else
            Nflusol=Nflusol+1
            Flusol(1,I)=2
            Flusol(2,I)=Nflusol
         end if
      else
         test=0      
         DO J=1,3
            IF (Neigh(J,I).gt.0) then
               if (Type_Media(Neigh(J,I)).eq.1) then
                  test=1
               end if
            end IF
         END DO
         if (test.eq.0) then 
            Nsol=Nsol+1
            Flusol(1,I)=4
            Flusol(2,I)=Nsol
         else
            Nsolflu=Nsolflu+1
            Flusol(1,I)=3
            Flusol(2,I)=Nsolflu
         end if
      end if
  ENDDO
  IF(Nflu+Nflusol+Nsolflu+Nsol.ne.Ntri) then
     write(6,*)'error',Nflu+Nflusol+Nsolflu+Nsol,Ntri
stop
  ENDIF
allocate(Triinter(Ntri,3))
 DO I=1,NTri
    SELECT CASE(Flusol(1,I))
    CASE(1)
       Triinter(Flusol(2,I),:)=Tri(I,:)
    CASE(2)
       Triinter(Nflu+Flusol(2,I),:)=Tri(I,:)
    CASE(3)
       Triinter(Nflu+Nflusol+Flusol(2,I),:)=Tri(I,:)
   CASE(4)
       Triinter(Nflu+Nflusol+Nsolflu+Flusol(2,I),:)=Tri(I,:)
    end SELECT
 END DO
Tri=Triinter
Triinter=0
 DO I=1,NTri
    SELECT CASE(Flusol(1,I))
    CASE(1)
       Triinter(Flusol(2,I),:)=Neigh(:,I)
    CASE(2)
       Triinter(Nflu+Flusol(2,I),:)=Neigh(:,I)
    CASE(3)
       Triinter(Nflu+Nflusol+Flusol(2,I),:)=Neigh(:,I)
   CASE(4)
       Triinter(Nflu+Nflusol+Nsolflu+Flusol(2,I),:)=Neigh(:,I)
    end SELECT
 END DO
Neigh(1,1:Ntri)=Triinter(1:Ntri,1)
Neigh(2,1:Ntri)=Triinter(1:Ntri,2)
Neigh(3,1:Ntri)=Triinter(1:Ntri,3)
deallocate(Triinter)
TYPE_MEDIA(1:Nflu+Nflusol)=1
TYPE_MEDIA(Nflu+Nflusol+1:Nflu+Nflusol+Nsolflu+Nsol)=2.
  DO I=1,Ntri
     DO J=1,3
        IF (Neigh(J,I).gt.0) then        
           SELECT CASE(Flusol(1,Neigh(J,I)))
        CASE(1)
           Neigh(J,I)=Flusol(2,Neigh(J,I))
        CASE(2)
           Neigh(J,I)=Nflu+Flusol(2,Neigh(J,I))
        CASE(3)
           Neigh(J,I)=Nflu+Nflusol+Flusol(2,Neigh(J,I))
       CASE(4)
           Neigh(J,I)=Nflu+Nflusol+Nsolflu+Flusol(2,Neigh(J,I))
        end SELECT
        end IF
     END DO
  ENDDO
END SUBROUTINE sub_renum2

SUBROUTINE my_sub_renum2
  Use m_mesh
  implicit none

  Nflu=Ntri
  Nsol=0     
  Nflusol=0
  Nsolflu=0
  
END SUBROUTINE my_sub_renum2
  
  

