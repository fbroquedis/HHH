subroutine sub_centre
  Use m_mesh
  Use m_gen
  
  implicit none
  
  integer                                :: i,node(3)
  real*8                                 :: x1,x2,x3,xm,y1,y2,y3,ym
  real*8,dimension(:,:),allocatable,save :: cgrav_ref
  real*8,save                            :: a(3,2)
  
  allocate (cgrav(ntri,2))
  allocate (cgrav_ref(ntri,2))
  allocate (ValPhicentre(ntri,nphi))

  do i=1,ntri
     node=tri(i,:)
     x1=coor(node(1),1)
     x2=coor(node(2),1)
     x3=coor(node(3),1)
     xm=(x1+x2+x3)/3
     cgrav(i,1)=xm
     
     y1=coor(node(1),2)
     y2=coor(node(2),2)
     y3=coor(node(3),2)
     ym=(y1+y2+y3)/3
     cgrav(i,2)=ym
     
     a=coor(tri(i,:),:)
     call sub_transform(a,cgrav(i,:),cgrav_ref(i,:))
     
     IF ((cgrav_ref(i,1).ge.0.d0).and.(cgrav_ref(i,2).ge.0.d0).and.(cgrav_ref(i,1)+cgrav_ref(i,2).le.1D0)) then
        SELECT CASE(ORDER)
        CASE(1)
           CALL Phi2DOrder1(cgrav_ref(i,1),cgrav_ref(i,2),1,ValPhicentre(i,:))     
        CASE(2)
           CALL Phi2DOrder2(cgrav_ref(i,1),cgrav_ref(i,2),1,ValPhicentre(i,:))     
        CASE(3)
           CALL Phi2DOrder3(cgrav_ref(i,1),cgrav_ref(i,2),1,ValPhicentre(i,:))     
        END SELECT
     END IF
  enddo

end subroutine sub_centre







