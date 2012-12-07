module m_mat

  implicit none
!!$  complex*16,dimension(:,:,:),allocatable,save:: A
!!$  complex*16,dimension(:,:,:),allocatable,save:: A_pp,A_uxux,A_uxuy,A_uyux,A_uyuy,A_pux,A_puy,A_uxp,A_uyp
  real*8,dimension(:,:,:),allocatable,save:: A
  real*8,dimension(:,:,:),allocatable,save:: A_pp,A_uxux,A_uxuy,A_uyux,A_uyuy,A_pux,A_puy,A_uxp,A_uyp
  complex*16, dimension(:), allocatable,save :: P_complexe,P_analytic, Ux_analytic,Ux_complexe,Uy_analytic,Uy_complexe
!  real*8,dimension(:,:,:),allocatable,save:: A
  real*8,dimension(:,:,:),allocatable,save    :: Dcla,B,G,S1,Dclabis,Gbis,Ebis
  real*8,dimension(:,:,:),allocatable,save    :: BBcla,Bclainv,BBclabis,Bbis
  real*8,dimension(:,:),allocatable,save   :: Minv,sqrtMinv,sqrtM,M,phiglob,phiglobold
  real*8,dimension(:),allocatable,save   :: Ux,Uy,phi1old,phi1,Xiold,Xi,UU,DFVEC,DFVEC_arete,U_inc,P
  real*8,dimension(:),allocatable,save   :: AP,AUx,AUy,Pold,Pinter,Uxold,Uxinter,Uyold,Uyinter
  real*8,dimension(:),allocatable,save   :: psiold,psi,psiinter,W,Wold,TriBordR0,var_aux
  real*8,dimension(:,:),allocatable,save  ::u_bord !variable auxiliaire pour le schéma avec dérivée fractionnaire en temps avec tout stocker
  
  real*8,dimension(:,:),allocatable,save,target :: Us
  real*8,dimension(:),allocatable,save          :: AU
  real*8,dimension(:),pointer                   :: U,Uold
    

  integer,dimension(:,:),allocatable :: Corres_arete,Corres_arete_u,Corres_arete_rampant
  integer,dimension(:),allocatable,save :: Corres_flusol,Corres_solflu
  integer,save :: NGL2D,NGL1D,NbR0
  real*8,allocatable,save :: wGL2D(:)
  real*8,allocatable,save :: PtGL2D(:,:)
  real*8,allocatable,save :: wGL1D(:)
  real*8,allocatable,save :: PtGL1D(:)
end module m_mat
