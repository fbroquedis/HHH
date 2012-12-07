SUBROUTINE sub_val_source

Use m_mesh
Use m_source

implicit none

real               :: lam,pi,t,tpeak,ots
integer            :: k,nnt!,fpeak


nnt=600005
allocate(Valsource(nnt))
allocate(Valsourcederiv2(nnt))
 pi=4.*atan(1.)

!fpeak=5
fpeak=1e6
 tpeak=0.
!tpeak=0
 Valsource=0.
 !fpeak=10
 lam=pi*pi*fpeak*fpeak
 ots=-1./fpeak
   do k=1,nnt
      t=ots+(k-1)*dt
      Valsource(k)=2.*lam*(2.*lam*(t-tpeak)*(t-tpeak)-1)*&
      exp(-lam*(t-tpeak)*(t-tpeak))
      Valsource(k)=-(2.*lam*(t-tpeak))*&
      exp(-lam*(t-tpeak)*(t-tpeak))
      
!      Valsourcederiv2(k)=4.*lam**2*exp(-lam*(t-tpeak)*(t-tpeak))*&
!      ((2*lam*(t-tpeak)**2-3)**2-6)     

!      open(20,File="source@.H",access='direct',recl=4)
!write(20,rec=k) real(valsource(k))
!    close(20)

   enddo
END SUBROUTINE sub_val_source
