SUBROUTINE  sub_transform(a,point,point1)

      implicit none
      real*8 :: point(2),a(3,2),pointbis(2),det,point1(2)

      point(:)=point(:)-a(1,:)
      pointbis(1)=(a(3,2)-a(1,2))*point(1)+(a(1,1)-a(3,1))*point(2)
      pointbis(2)=(a(1,2)-a(2,2))*point(1)+(a(2,1)-a(1,1))*point(2)
      det=(a(3,2)-a(1,2))*(a(2,1)-a(1,1))-(a(3,1)-a(1,1))*(a(2,2)-a(1,2))
      point1=pointbis/det
      point(:)=point(:)+a(1,:)

END SUBROUTINE sub_transform
