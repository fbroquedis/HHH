SUBROUTINE sub_writeximage
  Use m_mesh
  Use m_output
  Use m_gen
  Use m_mat
  implicit none
  INTEGER :: I,J,K
  character*100 :: nam
  real*8 :: Ecrit(Nphi)
      real*8,dimension (:,:),allocatable ::Stock
allocate(Stock(nbx,nby))
  write(nam,*) NbInst
  if (Nbinst.le.9) then
     nam="00"//trim(adjustl(nam))
  elseif (Nbinst.le.99) then
     nam="0"//trim(adjustl(nam))
  end if


   DO J=1,nbx
      DO K=1,nby
      I=Triinterp((K-1)*nbx+J)
      IF (I.gt.0) then
      Ecrit=U((I-1)*Nphi+1:I*Nphi)
      Stock(J,K)=sum(Ecrit*ValPhiGrid((K-1)*nbx+J,:))
   else
     Stock(J,K)=0.D0
  END IF
      ENDDO
   END DO
write(6,*) nbx,nby
  open(20,FILE="FILM/test."//trim(nam)//"@.H",access='direct',recl=4*nbx*nby)
write(20,rec=1) real(Stock(1:nbx,1:nby))
close(20)
deallocate(Stock)
END SUBROUTINE sub_writeximage
