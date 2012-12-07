module m_condbord

  implicit none

  integer,save :: CondBord,NbCla,isCla,Narete_du_bord,Ntri_non_bord
  integer,dimension(:),allocatable,save :: TriCla,CondPoint,Tri_bord,test_bord,tri_non_bord
  real*8,dimension (:,:,:),allocatable,save :: BCla,Bcla2,Bcla3,Bcla4

  integer,save :: Nbinc
  integer,dimension(:),allocatable,save :: Triinc,TriNoinc
  real*8,dimension (:,:,:),allocatable,save :: Binc
  real*8,dimension (:),allocatable,save :: courbure
end module m_condbord

