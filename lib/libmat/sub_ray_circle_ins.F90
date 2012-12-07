SUBROUTINE sub_ray_circle_ins(a,b,c,VAL)

  implicit none
  REAL*8 ::a(2),b(2),c(2),ab(2),ac(2),bc(2),VAL
  REAL*8 ::long_ab, long_ac, long_bc, aire,vol
  INTEGER ::I

  
  do I=1,2 
     ab(I)=b(I)-a(I)
     ac(I)=c(I)-a(I)
     bc(I)=c(I)-b(I)
  enddo

  long_ab=sqrt(sum(ab**2))
  long_ac=sqrt(sum(ac**2))
  long_bc=sqrt(sum(bc**2))

  aire=1/2.*(long_ab+long_ac+long_bc)
  aire=aire*(aire-long_ab)*(aire-long_ac)*(aire-long_bc)
  aire=sqrt(aire)

  VAL=2*aire/(long_ab+long_ac+long_bc)
END SUBROUTINE sub_ray_circle_ins
