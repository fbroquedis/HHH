SUBROUTINE sub_U0
  Use m_mesh
  Use m_gen
  Use m_mat
  Use m_source
  implicit none
  real*8 :: V12(2),V23(2),V31(2),x,y
  integer :: Node(3),i,j,k

DO I=1,Ntri
   DO J=1,3
      x=Coor(Tri(I,J),1)
      y=Coor(Tri(I,J),2)
!      IF ((x-x0)**2+(y-y0)**2.le. 2*r**2) THEN
         Uold(Nphi*(I-1)+J)=exp(-((x-x0)**2+(y-y0)**2)/r**2)
!      END IF
   ENDDO
END DO
IF (Order.ge.2) then
DO I=1,Ntri
   Node = Tri(I,:)
   V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
   V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
   V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
   V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
   V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
   V31(2)=Coor(Node(1),2)-Coor(Node(3),2)
   DO K=1,Order-1
      x=Coor(Tri(I,1),1)+K*V12(1)/Order
      y=Coor(Tri(I,1),2)+K*V12(2)/Order
      !      IF ((x-x0)**2+(y-y0)**2.le. 2*r**2) THEN
      Uold(Nphi*(I-1)+3+K)=exp(-((x-x0)**2+(y-y0)**2)/r**2)
      !      END IF
       x=Coor(Tri(I,2),1)+K*V23(1)/Order
      y=Coor(Tri(I,2),2)+K*V23(2)/Order
      !      IF ((x-x0)**2+(y-y0)**2.le. 2*r**2) THEN
      Uold(Nphi*(I-1)+3+Order-1+K)=exp(-((x-x0)**2+(y-y0)**2)/r**2)
      !      END IF
      x=Coor(Tri(I,3),1)+K*V31(1)/Order
      y=Coor(Tri(I,3),2)+K*V31(2)/Order
      !      IF ((x-x0)**2+(y-y0)**2.le. 2*r**2) THEN
      Uold(Nphi*(I-1)+3+2*(Order-1)+K)=exp(-((x-x0)**2+(y-y0)**2)/r**2)
      !      END IF
   ENDDO
   IF (Order.eq.3) then
      x=Coor(Tri(I,1),1)+V12(1)/3.-V31(1)/3.
      y=Coor(Tri(I,1),2)+V12(2)/3.-V31(2)/3.
      !      IF ((x-x0)**2+(y-y0)**2.le. 2*r**2) THEN
      Uold(Nphi*(I-1)+10)=exp(-((x-x0)**2+(y-y0)**2)/r**2)
   ENDIF
END DO
ENDIF
DO I=1,Ntri
U(Nphi*(I-1)+1:Nphi*I)=Uold(Nphi*(I-1)+1:Nphi*I)
END DO
Uold=U
!!$DO I=1,Ntri
!!$   J=Corres(I,1)
!!$   UU(Corres(I,2:Nphi+1))=MATMUL(A(I,1:Nphi,1:J),U(Corres(I,2:J+1)))
!!$ENDDO
END SUBROUTINE sub_U0
