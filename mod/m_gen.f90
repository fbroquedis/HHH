module m_gen

  implicit none

  INTEGER,save :: Order,Nphi,reff,hang,nangle,helmholtz
  real*8,dimension(:),allocatable,save :: anglealpha
  real*8,save ::anglealpha1,omega
  real*8,dimension(:,:),allocatable,save :: cgrav, ValPhicentre,coeffl

  real*8,save :: k1,k2,rayon_obstacle,omega_incident
end module m_gen