SUBROUTINE sub_writemeshCLA_rampant
  Use m_condbord
!  Use m_cla
  Use m_mesh
  Use m_mat
  implicit none
  INTEGER :: I,J,K

  open(20,FILE="FILM/mesh.vtu")
  write(20,FMT='(A21)')'<?xml version="1.0"?>' 
  write(20,*) '<VTKFile type="UnstructuredGrid">' 
  write(20,*) '<UnstructuredGrid>' 
  write(20,*) '<Piece NumberOfPoints="',3*Ntri,'" NumberOfCells="',Ntri,'">'
  write(20,*) '<Points>'
  write(20,*) '<DataArray type="Float64" Name="Position" NumberOfComponents="3" format="ascii">'
  DO K=1,Narete_du_bord
     I=tri_bord(K)
     DO J=1,3
        write(20,*) Coor(Tri(I,J),1),Coor(Tri(I,J),2), 0
     END DO
  ENDDO
  DO K=1,Ntri-Narete_du_bord
     I=tri_non_bord(K)
     DO J=1,3
        write(20,*) Coor(Tri(I,J),1),Coor(Tri(I,J),2), 0
     END DO
  ENDDO
  write(20,*) '</DataArray>'
  write(20,*) '</Points>'
  write(20,*) '<Cells>'
  write(20,*) '<DataArray type="Int32" Name="connectivity"   format="ascii">'
   DO K=1,Narete_du_bord
     I=tri_bord(K)
     write(20,*) 3*(I-1), 3*(I-1)+1, 3*(I-1)+2 
  END DO
   DO K=1,Ntri-Narete_du_bord
     I=tri_non_bord(K)
     write(20,*) 3*(I-1), 3*(I-1)+1, 3*(I-1)+2 
  END DO
  write(20,*)'</DataArray>'
  write(20,*)'<DataArray type="Int32" Name="offsets"  format="ascii">'
  DO I=1,NTri
     write(20,*) 3*I
  END DO
  write(20,*)'</DataArray>'
  write(20,*)' <DataArray type="UInt8"  Name="types"  format="ascii">'
  DO I=1,NTri
     write(20,*) 5
  END DO
  write(20,*)'</DataArray>'
  write(20,*) '</Cells>'
  write(20,*) '<PointData Scalars="U">'
  write(20,*) '<DataArray type="Float64" Name="U" format="ascii">'
   DO K=1,Narete_du_bord
     write(20,*) 1
     write(20,*) 1
     write(20,*) 1
  END DO
   DO K=1,Ntri-Narete_du_bord
     write(20,*) 2
     write(20,*) 2
     write(20,*) 2
  END DO
  write(20,*) '</DataArray>'
  write(20,*) '</PointData>'
  write(20,*) '</Piece>'
  write(20,*) '</UnstructuredGrid>' 
  write(20,*) '</VTKFile>' 
  close(20)
END SUBROUTINE sub_writemeshCLA_rampant
