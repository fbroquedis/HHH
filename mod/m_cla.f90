module m_cla

  implicit none

  integer,save :: NbCla
  integer,dimension(:),allocatable,save :: TriCla,Tri_bord
  real*8,dimension (:,:,:),allocatable,save :: BCla,Bcla2
  !real*8,save :: anglealpha

end module m_cla

