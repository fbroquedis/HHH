SUBROUTINE sub_renum
  Use m_mesh
  implicit none
  INTEGER :: I,J,K,Node(6),Test,compf,comp,prem
  INTEGER :: Iraff,Inoraff,Iboundraff,Iboundnoraff,Ineigh
  INTEGER,allocatable :: Testraff(:),Testraffinter(:),Corres(:),TestK(:),Testnoderaff(:),Stock(:)
  integer,dimension(:,:),allocatable,save :: TriInter,NeighInter
  real*8,dimension (:),allocatable,save :: Type_MediaInter
  real*8,dimension (:,:),allocatable,save :: CoorInter
  real*8 :: V12(2),V13(2),V23(2),x,y,h1
  Nraff=0
  Nnoraff=0
  allocate(Testraff(Ntri))
  allocate(Testnoderaff(NPoints))
  allocate(Stock(Ntri))    
  hraff=0
  CFLraff=0
  h=0
  CFL=0
  Testnoderaff=1
  compf=0
  prem=1    

  DO I=1,Ntri
      if(Type_Media(i)>c_lim) then
        Testraff(I)=1
        Nraff=Nraff+1 
        Stock(Nraff)=i
     ELSE
        Testraff(I)=0
        Nnoraff=Nnoraff+1
        Testnoderaff(Node(1:3))=0
 
     END IF
  ENDDO
  IF(Nraff+Nnoraff.ne.Ntri) then
     write(6,*)'error',Nraff+Nnoraff,Ntri
  ENDIF
  do i=1,dec
      comp=0    
      do j=prem,Nraff
         do k=1,3
            if (testraff(Neigh(Stock(j),k)).eq.0) then
               testraff(Neigh(Stock(j),k))=1
               comp=comp+1
               Stock(Nraff+comp)=Neigh(Stock(j),k)
            endif
         enddo
      enddo
         prem=Nraff+1
         Nraff=Nraff+comp
         compf=compf+comp
  enddo

Nnoraff=Nnoraff-compf
deallocate(Stock)
  IF(Nraff+Nnoraff.ne.Ntri) then
     write(6,*)'error',Nraff+Nnoraff,Ntri
  ENDIF

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
     IF (TESTRAFF(I).eq.1) THEN
        IF(hraff.eq.0) then 
           hraff=h1
           CFLraff=h1/TYPE_MEDIA(I)
        ELSE
           hraff=min(hraff,h1)
           CFLraff=min(CFLraff,h1/TYPE_MEDIA(I))
        ENDIF
     ELSE
        IF(h.eq.0) then 
           h=h1
           CFL=h/TYPE_MEDIA(I)
        ELSE
           h=min(h,h1)
           CFL=min(CFL,h1/TYPE_MEDIA(I))
        ENDIF
        DO J=1,3
           IF (NEIGH(I,J).gt.0) THEN
              IF (TESTRAFF(NEIGH(I,J)).eq.1) THEN
                 IF(hraff.eq.0) then 
                    hraff=h1
                    CFLraff=h1/TYPE_MEDIA(I)
                 ELSE
                    hraff=min(hraff,h1)
                    CFLraff=min(CFLraff,h1/TYPE_MEDIA(I))
                 ENDIF
              END IF
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  
  allocate(TriInter(Ntri,3))
  allocate(Type_MediaInter(Ntri))    
  allocate(NeighInter(Ntri,3))
  allocate(Corres(Ntri))
  Iraff=Nnoraff
  Inoraff=0
  DO I=1,Ntri
     IF(testraff(I).eq.0) then
        Inoraff=Inoraff+1
        TriInter(Inoraff,:)=Tri(I,:)
        Type_MediaInter(Inoraff)=TYPE_MEDIA(I)
        NeighInter(Inoraff,:)=Neigh(I,:)
        Corres(I)=Inoraff
     else
        Iraff=Iraff+1
        TriInter(Iraff,:)=Tri(I,:)
        Type_MediaInter(Iraff)=TYPE_MEDIA(I)
        NeighInter(Iraff,:)=Neigh(I,:)
        Corres(I)=Iraff
     ENDIF
  ENDDO
  Tri=TriInter
  TYPE_MEDIA=Type_MediaInter
  deallocate(TriInter)
  deallocate(Type_MediaInter)
  DO I=1,Ntri
     DO J=1,3
        IF (NeighInter(I,J).gt.0) then
           NeighInter(I,J)=Corres(NeighInter(I,J))
        ENDIF
     ENDDO
  ENDDO
  Neigh=NeighInter
  deallocate(NeighInter)
write(6,*) h,hraff,CFL,CFLraff

END SUBROUTINE sub_renum
