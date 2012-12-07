SUBROUTINE sub_writeregular_vero
  Use m_mesh
  Use m_output
  Use m_gen
  Use m_mat
  implicit none
  INTEGER :: I,J
  character*100 :: nam
  real*8 :: Ecrit(Nphi)
  write(nam,*) NbInst2
  if (Nbinst2.le.9) then
     nam="00"//trim(adjustl(nam))
  elseif (Nbinst2.le.99) then
     nam="0"//trim(adjustl(nam))
  end if

  open(20,FILE="FILM_donut_3_1_raff_source_0_1_3/evanescent/beta_1_fois_courb/test."//trim(nam)//".vti")
!  write(20,FMT='(A21)')'<?xml version="1.0"?>' 
!   write(20,*) '<VTKFile type="ImageData">' 
!   write(20,*) '<ImageData WholeExtent="',1,' ',nbx,' ',1,' ', nby,' 0 0"'
!write(20,*) 'Origin="0 0 0" Spacing="',stepx,' ', stepy,' ',' 0">' 
!   write(20,*) '<Piece Extent="',1,' ',nbx,' ',1,' ', nby,'0 0">'
!   write(20,*) '<PointData Scalars="U">'
!   write(20,*) '<DataArray type="Float64" Name="U" format="ascii">'
   DO J=1,nbx*nby
      I=Triinterp(J)
      IF (I.gt.0) then
      Ecrit=U((I-1)*Nphi+1:I*Nphi)
      write(20,*) sum(Ecrit*ValPhiGrid(J,:))
   else
     write(20,*) 0.D0
  END IF
   END DO
!   write(20,*) '</DataArray>'
!   write(20,*) '</PointData>'
!   write(20,*) '</Piece>'
!   write(20,*) '</ImageData>' 
!   write(20,*) '</VTKFile>' 
close(20)
END SUBROUTINE sub_writeregular_vero
