SUBROUTINE  sub_defgrid_modif
  Use m_mesh 
  Use m_gen
  implicit none
  integer :: i,ixmin,ixmax,iymin,iymax,ix,iy,nb
  real*8 :: xmin,xmax,ymin,ymax,point(2),a(3,2),pointbis(2),det
  nbx=(xgrid2-xgrid1)/stepx+1
  nby=(ygrid2-ygrid1)/stepy+1
  allocate(Triinterp(nbx*nby))
  allocate(ValphiGrid(nbx*nby,Nphi))
  Triinterp=0
  DO I=1,Ntri
     a=Coor(Tri(I,:),:)
     IF (((a(1,1)**2/9.D0+a(1,2)**2/9.D0).le.1.01).and.((a(2,1)**2/9.D0+a(2,2)**2/9.D0).le.1.01)&
			.and.((a(3,1)**2/9.D0+a(3,2)**2/9.D0).le.1.01)) then
!write(6,*) I
     xmin=min(a(1,1),a(2,1),a(3,1))
     xmax=max(a(1,1),a(2,1),a(3,1))
     ymin=min(a(1,2),a(2,2),a(3,2))
     ymax=max(a(1,2),a(2,2),a(3,2))
     xmin=xmin*(1-1E-5)
     xmax=xmax*(1+1E-5)
     ymin=ymin*(1-1E-5)
     ymax=ymax*(1+1E-5)
     ixmin=max(INT((xmin-xgrid1)/stepx+1)-1,1)
     ixmax=min(INT((xmax-xgrid1)/stepx+1)+2,nbx)
     iymin=max(INT((ymin-ygrid1)/stepy+1)-1,1)
     iymax=min(INT((ymax-ygrid1)/stepy+1)+2,nby)
     IF ((ixmax.ge.1).and.(iymax.ge.1)) then
        DO iy=iymin,iymax
           DO ix=ixmin,ixmax
              nb=(iy-1)*nbx+ix
              point(1)=xgrid1+(ix-1)*stepx
              point(2)=ygrid1+(iy-1)*stepy
!              write(6,*) point
!!$              write(6,*) a(1,:)
!!$              write(6,*) a(2,:)
!!$              write(6,*) a(3,:)
!              write(6,*) xmin,xmax,ymin,ymax
!              write(6,*) ixmin,ixmax,iymin,iymax
!!$              write(6,*) xgrid1,(ix-1)*stepx
!              write(6,*) 'bedew'
              point(:)=point(:)-a(1,:)
              pointbis(1)=(a(3,2)-a(1,2))*point(1)+(a(1,1)-a(3,1))*point(2)
              pointbis(2)=(a(1,2)-a(2,2))*point(1)+(a(2,1)-a(1,1))*point(2)
              det=(a(3,2)-a(1,2))*(a(2,1)-a(1,1))-(a(3,1)-a(1,1))*(a(2,2)-a(1,2))
              point=pointbis/det
!     write(6,*) point
              IF ((point(1).ge.0.d0).and.(point(2).ge.0.d0).and.(point(1)+point(2).le.1D0)) then
                 Triinterp(nb)=i
                 SELECT CASE(ORDER)
                 CASE(1)
                    CALL Phi2DOrder1(point(1),point(2),1,ValPhigrid(nb,:))     
                 CASE(2)
                    CALL Phi2DOrder2(point(1),point(2),1,ValPhigrid(nb,:))     
                 CASE(3)
                    CALL Phi2DOrder3(point(1),point(2),1,ValPhigrid(nb,:))     
                 END SELECT
              END IF
           END DO
        ENDDO
     ENDIF
     ENDIF
  END DO
END SUBROUTINE sub_defgrid_modif
