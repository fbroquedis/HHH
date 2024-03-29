SUBROUTINE sub_writeunstruct_p
  Use m_mesh
  Use m_output
  Use m_condbord
  Use m_gen
  Use m_mat
  implicit none
  INTEGER :: I,J,K
  character*100 :: nam
  real*8 :: Ecrit(Nphi)
  real*8 :: DF,V12(2),V23(2),V31(2),x,y,N(2,2),Node2(2),JJ
  integer :: Node(3),Ncells,Nnodes
  real*8 ::Ecritbis,h5
  If (Order.eq.3) then
     Ncells=9*(Nflu+Nflusol)
     Nnodes=Nphi*(Nflu+Nflusol)
     else
        Ncells=(Nflu+Nflusol)
        Nnodes=3*(Nflu+Nflusol)
     Endif
if(Nbinst.eq.0) then
  open(20,FILE="FILM/P_real.vtu")
else
  open(20,FILE="FILM/P_imag.vtu")
end if
write(20,FMT='(A21)')'<?xml version="1.0"?>' 
   write(20,*) '<VTKFile type="UnstructuredGrid">' 
   write(20,*) '<UnstructuredGrid>' 
   write(20,*) '<Piece NumberOfPoints="',Nnodes,'" NumberOfCells="',Ncells,'">'
   write(20,*) '<Points>'
write(20,*) '<DataArray type="Float64" Name="Position" NumberOfComponents="3" format="ascii">'
DO I=1,Nflu+Nflusol
   DO J=1,3
      write(20,*) Coor(Tri(I,J),1),Coor(Tri(I,J),2), 0
   END DO
   IF (Order.eq.3) then
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
         write(20,*) x,y, 0
      ENDDO
      DO K=1,Order-1
         x=Coor(Tri(I,2),1)+K*V23(1)/Order
         y=Coor(Tri(I,2),2)+K*V23(2)/Order
         write(20,*) x,y, 0
      ENDDO
      DO K=1,Order-1
         x=Coor(Tri(I,3),1)+K*V31(1)/Order
         y=Coor(Tri(I,3),2)+K*V31(2)/Order
         write(20,*) x,y, 0
      ENDDO
      x=Coor(Tri(I,1),1)+V12(1)/3.-V31(1)/3.
      y=Coor(Tri(I,1),2)+V12(2)/3.-V31(2)/3.
         write(20,*) x,y, 0
   ENDIF
ENDDO
   write(20,*) '</DataArray>'
   write(20,*) '</Points>'
   write(20,*) '<Cells>'
   write(20,*) '<DataArray type="Int32" Name="connectivity"   format="ascii">'
IF (Order.eq.3) then
   DO I=1,Nflu+Nflusol
      write(20,*) Nphi*(I-1), Nphi*(I-1)+3,Nphi*(I-1)+8
      write(20,*) Nphi*(I-1)+3, Nphi*(I-1)+9,Nphi*(I-1)+8
      write(20,*) Nphi*(I-1)+3,Nphi*(I-1)+4,Nphi*(I-1)+9
      write(20,*) Nphi*(I-1)+4,Nphi*(I-1)+5,Nphi*(I-1)+9
      write(20,*) Nphi*(I-1)+4,Nphi*(I-1)+1,Nphi*(I-1)+5
      write(20,*) Nphi*(I-1)+8,Nphi*(I-1)+9,Nphi*(I-1)+7
      write(20,*) Nphi*(I-1)+9,Nphi*(I-1)+6,Nphi*(I-1)+7
      write(20,*) Nphi*(I-1)+9,Nphi*(I-1)+5,Nphi*(I-1)+6
      write(20,*) Nphi*(I-1)+7,Nphi*(I-1)+6,Nphi*(I-1)+2
   END DO
ELSE
   DO I=1,Nflu+Nflusol
      write(20,*) Nphi*(I-1), Nphi*(I-1)+1,Nphi*(I-1)+2
   END DO
END IF
   write(20,*)'</DataArray>'
   write(20,*)'<DataArray type="Int32" Name="offsets"  format="ascii">'
IF (Order.eq.3) then
   DO I=1,Nflu+Nflusol
      write(20,*) 27*(I-1)+3, 27*(I-1)+6, 27*(I-1)+9, 27*(I-1)+12, 27*(I-1)+15, 27*(I-1)+18, 27*(I-1)+21, 27*(I-1)+24, 27*(I-1)+27
   END DO
else
   DO I=1,Nflu+Nflusol
     write(20,*) 3*I
   END DO
ENDIF
   write(20,*)'</DataArray>'
   write(20,*)' <DataArray type="UInt8"  Name="types"  format="ascii">'
IF (Order.eq.3) then
   DO I=1,Nflu+Nflusol
      write(20,*) 5,5,5,5,5,5,5,5,5
   END DO
ELSE
   DO I=1,Nflu+Nflusol
      write(20,*) 5
   END DO
ENDIF
   write(20,*)'</DataArray>'
   write(20,*) '</Cells>'
   write(20,*) '<PointData Scalars="U">'
   write(20,*) '<DataArray type="Float64" Name="U" format="ascii">'
IF (Order.eq.3) then
   DO I=1,(Nflu+Nflusol)
      DO J=1,Nphi
      write(20,*) P((I-1)*Nphi+J)
   END DO
   END DO
else
   DO I=1,(Nflu+Nflusol)
      DO J=1,3
      write(20,*) P((I-1)*Nphi+J)
   END DO
   END DO
endif
   write(20,*) '</DataArray>'
   write(20,*) '</PointData>'
  write(20,*) '</Piece>'
   write(20,*) '</UnstructuredGrid>' 
   write(20,*) '</VTKFile>' 
 close(20)
     
END SUBROUTINE sub_writeunstruct_p
