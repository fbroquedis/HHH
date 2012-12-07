SUBROUTINE sub_source_elasto_acoustic_localisation
      Use m_mesh
      Use m_gen
      Use m_source
      implicit none
      real*8,save :: coord_s(2),a(3,2),coord_s_ref(2)
      integer :: I,node(3),J
REAL*8 :: Jfinv(2,2),V(3,2)
!      x_source=0
 !     y_source=0.9

      coord_s(1)=x0
      coord_s(2)=y0
      DO I=1,Ntri
         a=Coor(Tri(I,:),:)
         call sub_transform(a,coord_s,coord_s_ref)
         IF ((coord_s_ref(1).ge.0.d0).and.(coord_s_ref(2).ge.0.d0).and.(coord_s_ref(1)+coord_s_ref(2).le.1D0)) then
            write(6,*)'coucou'
            Isource=I
if(Isource.le.Nflu+Nflusol) then
      allocate(ValPhisource(Nphi,1))

	    SELECT CASE(ORDER)
            CASE(1)
               CALL Phi2DOrder1(coord_s_ref(1),coord_s_ref(2),1,ValPhisource)     
	    CASE(2)
               CALL Phi2DOrder2(coord_s_ref(1),coord_s_ref(2),1,ValPhisource)     
            CASE(3)
               CALL Phi2DOrder3(coord_s_ref(1),coord_s_ref(2),1,ValPhisource)     

            END SELECT
else
      allocate(ValPhisource(Nphi,2))

	    SELECT CASE(ORDER)
            CASE(1)
               CALL GradPhiOrder1(coord_s_ref(1),coord_s_ref(2),1,ValPhisource)     
	    CASE(2)
               CALL GradPhiOrder2(coord_s_ref(1),coord_s_ref(2),1,ValPhisource)     
            CASE(3)
               CALL GradPhiOrder3(coord_s_ref(1),coord_s_ref(2),1,ValPhisource)     

            END SELECT
Node = Tri(Isource,:)
   V(1,1)=Coor(Node(3),1)-Coor(Node(2),1)
   V(1,2)=Coor(Node(3),2)-Coor(Node(2),2)
   V(2,1)=Coor(Node(1),1)-Coor(Node(3),1)
   V(2,2)=Coor(Node(1),2)-Coor(Node(3),2)
   V(3,1)=Coor(Node(2),1)-Coor(Node(1),1)
   V(3,2)=Coor(Node(2),2)-Coor(Node(1),2)
Jfinv(1,1)=-V(2,2)
     Jfinv(2,2)=V(3,1)
     Jfinv(1,2)=-V(3,2)
     Jfinv(2,1)=V(2,1)
     Jfinv=(Jfinv/(-V(3,1)*V(2,2)+V(3,2)*V(2,1)))
do j=1,Nphi
Valphisource(j,:)=matmul(Jfinv,Valphisource(j,:))
enddo
end if
            !write(6,*)'valphisource',ValPhisource
            return
         END IF
      ENDDO
    END SUBROUTINE sub_source_elasto_acoustic_localisation
