subroutine sub_write_sismos(nt)

      Use m_mat
      Use m_output
      Use m_gen
      Use m_source

      implicit none
      character*100 :: nam
      real*8 :: Ecrit(Nphi)
      real*8,dimension(:,:),allocatable,save::Stock(:)
      integer :: nt,j,i


      allocate(Stock(nbrcv))
	
      nbsis=nbsis+1

      do i=1,nbrcv
         j=ircv(i)
!write(6,*) i,j
         Ecrit=U((j-1)*Nphi+1:j*Nphi)

 Stock(i)=sum(Ecrit*Valphircv(i,:))

      enddo

         write(21,*) Stock

         

deallocate(Stock)
end subroutine sub_write_sismos
      
