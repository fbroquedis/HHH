module m_source

  implicit none
  real*8,save ::x0,y0,z0,r,x_rcv,y_rcv,fpeak
  real*8,dimension (:,:),allocatable,save :: ValPhisource,ValPhircv
  real*8,dimension (:),allocatable,save :: Valsource,Valsourcederiv2
  integer,save ::Isource
  integer,dimension(:),allocatable,save  :: Ircv
  integer, save :: nbrcv,nbrcvflu,nbrcvsol
  real*8,dimension(:,:),allocatable,save :: coord_rcv

end module m_source
