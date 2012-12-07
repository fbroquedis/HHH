subroutine sub_write_sismos_acoustic_elasto(nt)

      Use m_mat
      Use m_output
      Use m_gen
      Use m_source
      Use m_mesh
      implicit none
      character*100 :: nam
      real*8 :: Ecrit(Nphi)
      real*8,dimension(:,:),allocatable,save::StockP(:)
      real*8,dimension(:,:),allocatable,save::StockUx(:)
      real*8,dimension(:,:),allocatable,save::StockUy(:)
      integer :: nt,j,i,iflu,isol


      allocate(StockP(nbrcvflu))
      allocate(StockUx(nbrcvsol))
      allocate(StockUy(nbrcvsol))
	
      nbsis=nbsis+1
iflu=0
isol=0
      do i=1,nbrcv
         j=ircv(i)
!write(6,*) i,j
         if (j.le.nflu+nflusol) then
            iflu=iflu+1
            Ecrit=P((j-1)*Nphi+1:j*Nphi)        
            StockP(iflu)=sum(Ecrit*Valphircv(i,:))
      else
         isol=isol+1
         j=j-nflu-nflusol
            Ecrit=Ux((j-1)*Nphi+1:j*Nphi)        
            StockUx(isol)=sum(Ecrit*Valphircv(i,:))
            Ecrit=Uy((j-1)*Nphi+1:j*Nphi)        
            StockUy(isol)=sum(Ecrit*Valphircv(i,:))
         end if
      enddo

         write(23,*) StockP
         write(24,*) StockUx
         write(25,*) StockUy

         

deallocate(StockP)
deallocate(StockUx)
deallocate(StockUy)
end subroutine sub_write_sismos_acoustic_elasto
      
