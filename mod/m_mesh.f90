module m_mesh

  implicit none
  integer,save ::NPoints,Ntri,Nneigh,Nnoraff,Nraff,Nboundraff,Nboundnoraff
  integer,save :: Nflu,Nsol,Nflusol,Nsolflu
  integer,dimension(:,:),allocatable,save :: Tri,Neigh,Neigh_bord,Arete_du_bord
  real*8,dimension (:,:),allocatable,save :: Coor
  real*8,dimension (:),allocatable,save :: mu,rho
  integer,dimension (:),allocatable,save :: Type_Media
  real*8,dimension (:,:,:),allocatable,save :: Cij
  real*8, save :: c_lim
  integer, save :: dec
  real*8,save ::  h,hraff,dt,x1reff,y1reff,x2reff,y2reff,CFL,CFLraff
  real*8,save ::  xgrid1,xgrid2,stepx,ygrid1,ygrid2,stepy
  integer,save ::nbx,nby
  real*8,dimension (:,:),allocatable,save :: ValphiGrid
  integer,dimension(:),allocatable,save :: Triinterp

!! tout ce qui relatif au rampant
  integer :: nombre_noeuds_DF
  real*8, dimension(:), allocatable :: pas_espace_DF,coeff_der_frac,W2_new_deriv
  integer, dimension(:), allocatable :: liste_sommets_bord
  real*8, dimension(:,:), allocatable :: liste_noeuds_DF
  real*8 :: pas_espace_DF_max
complex*16,save::alpha,alpha2
end module m_mesh

