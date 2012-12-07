subroutine sub_write_val_centre
  use m_mesh
  use m_output
  use m_gen
  use m_mat
  
  implicit none
  real*8         :: Ecrit(Nphi)
  real*8         :: Stock
  integer        :: i
   character*100 :: nam

  write(nam,*) NbInst
  if (Nbinst.le.9) then
     nam="00"//trim(adjustl(nam))
  elseif (Nbinst.le.99) then
     nam="0"//trim(adjustl(nam))
  end if

  open(20,File="FILM/Valcentre."//trim(nam)//"@.H",access='direct',recl=4)

  do i=1,ntri
     Ecrit=U((i-1)*Nphi+1:i*Nphi)
     Stock=sum(Ecrit*Valphicentre(i,:))
     
     write(20,rec=i) real(stock)
  enddo

  close(20)

end subroutine sub_write_val_centre
