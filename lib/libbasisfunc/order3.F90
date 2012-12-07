!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Order 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Phi2DOrder3(X,Y,DIM,VAL)
!!! Phi1=9*(1-x-y)(2/3-x-y)(1/3-x-y)/2
!!! Phi2=9x(x-1/3)(x-2/3)/2
!!! Phi3=9y(y-1/3)(y-2/3)/2
!!! Phi4=27x(1-x-y)(2/3-x-y)/2
!!! Phi5=27x(x-1/3)(1-x-y)/2
!!! Phi6=27x(x-1/3)y/2
!!! Phi7=27xy(y-1/3)/2
!!! Phi8=27y(y-1/3)(1-x-y)/2
!!! Phi9=27y(1-x-y)(2/3-x-y)/2
!!! Phi10=27xy(1-x-y)
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X,Y,T
  real*8, dimension(DIM):: X1,X2,Y1,Y2,T1,T2
  real*8, dimension(10,DIM):: VAL
  T=1D0-X-Y
  T1=2/3D0-X-Y
  T2=1/3D0-X-Y
  X1=X-1/3D0
  X2=X-2/3D0
  Y1=Y-1/3D0
  Y2=Y-2/3D0
  VAL(1,:)=T*T1*T2
  VAL(2,:)=X*X1*X2
  VAL(3,:)=Y*Y1*Y2
  VAL(1:3,:)=9/2D0*VAL(1:3,:)
  VAL(4,:)=X*T*T1
  VAL(5,:)=X*X1*T
  VAL(6,:)=X*Y*X1
  VAL(7,:)=X*Y*Y1
  VAL(8,:)=Y*Y1*T
  VAL(9,:)=Y*T*T1
  VAL(4:9,:)=27/2D0*VAL(4:9,:)
  VAL(10,:)=27D0*X*Y*T
END SUBROUTINE Phi2DOrder3

SUBROUTINE Phi1DOrder3(X,DIM,VAL)
!!! Phi1=9/2(1-x)(2/3-x)(1/3-x)
!!! Phi2=27/2(1-x)(2/3-x)x
!!! Phi3=27/2x(x-1/3)(1-x)
!!! Phi4=9/2x(x-1/3)(x-2/3)
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X,X1,X2,X3
  real*8, dimension(4,DIM):: VAL
  X1=1/3D0-X
  X2=2/3D0-X
  X3=1D0-X
  VAL(1,:)=9/2D0*X1*X2*X3
  VAL(2,:)=27/2D0*X*X2*X3
  VAL(3,:)=-27/2D0*X*X1*X3
  VAL(4,:)=9/2D0*X*X1*X2
END SUBROUTINE Phi1DOrder3

SUBROUTINE GradPhi1DOrder3(X,DIM,VAL)
!!! Phi1=9/2(1-x)(2/3-x)(1/3-x)
!!! Phi2=27/2(1-x)(2/3-x)x
!!! Phi3=27/2x(x-1/3)(1-x)
!!! Phi4=9/2x(x-1/3)(x-2/3)
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X,X1,X2,X3
  real*8, dimension(4,DIM):: VAL
  X1=1/3D0-X
  X2=2/3D0-X
  X3=1D0-X
  VAL(1,:)=-9/2D0*(X1*X2+X2*X3+X1*X3)
  VAL(2,:)=27/2D0*(X2*X3-X*X2-X*X3)
  VAL(3,:)=-27/2D0*(X1*X3-X*X1-X*X3)
  VAL(4,:)=9/2D0*(X1*X2-X*X1-X*X2)
END SUBROUTINE GradPhi1DOrder3

SUBROUTINE GradPhiOrder3(X,Y,DIM,VAL)
!!! GradPhi1x=-9/2*[(2/3-x-y)(1/3-x-y)+(1-x-y)(1/3-x-y)+(1-x-y)(2/3-x-y)]
!!! GradPhi1y=-9/2*[(2/3-x-y)(1/3-x-y)+(1-x-y)(1/3-x-y)+(1-x-y)(2/3-x-y)]
!!! GradPhi2x=9/2*[(x-1/3)(x-2/3)+x(x-1/3)+x(x-2/3)]
!!! GradPhi2y=0
!!! GradPhi3x=0
!!! GradPhi3y=9*[(y-1/3)(y-2/3)+y(y-1/3)+y(y-2/3)]
!!! GradPhi4x=27/2*[(1-x-y)(2/3-x-y)-x(2/3-x-y)-x(1-x-y)]
!!! GradPhi4y=-27/2*x*[(2/3-x-y)+(1-x-y)]
!!! GradPhi5x=27/2*[(x-1/3)(1-x-y)+x(1-x-y)-x(x-1/3)]
!!! GradPhi5y=-27/2x(x-1/3)
!!! GradPhi6x=27/2y(2x-1/3)
!!! GradPhi6y=27/2x(x-1/3)
!!! GradPhi7x=27/2y(y-1/3)
!!! GradPhi7y=27/2x(2y-1/3)
!!! GradPhi8x=-27/2y(y-1/3)
!!! GradPhi8y=27/2[(y-1/3)(1-x-y)+y(1-x-y)-y(y-1/3)]
!!! GradPhi9x=-27/2y[(2/3-x-y)+(1-x-y)]
!!! GradPhi9y=27/2[(1-x-y)(2/3-x-y)-y(2/3-x-y)-y(1-x-y)]
!!! GradPhi10x=27y(1-2x-y)
!!! GradPhi10y=27x(1-x-2y)
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X,Y,T
  real*8, dimension(DIM):: X1,X2,Y1,Y2,T1,T2
  real*8, dimension(10,2,DIM):: VAL
  T=1D0-X-Y
  T1=2/3D0-X-Y
  T2=1/3D0-X-Y
  X1=X-1/3D0
  X2=X-2/3D0
  Y1=Y-1/3D0
  Y2=Y-2/3D0
  VAL(1,1,:)=-(T2*T1+T*T1+T*T2)
  VAL(1,2,:)=-(T2*T1+T*T1+T*T2)
  VAL(2,1,:)=X1*X2+X*X1+X*X2
  VAL(2,2,:)=0D0
  VAL(3,1,:)=0D0
  VAL(3,2,:)=Y1*Y2+Y*Y1+Y*Y2
  VAL(1:3,:,:)=9/2D0*VAL(1:3,:,:)
  VAL(4,1,:)=T*T1-X*T1-X*T
  VAL(4,2,:)=-X*(T1+T)
  VAL(5,1,:)=X1*T+X*T-X*X1
  VAL(5,2,:)=-X*X1
  VAL(6,1,:)=Y*(X+X1)
  VAL(6,2,:)=X*X1
  VAL(7,1,:)=Y*Y1
  VAL(7,2,:)=X*(Y+Y1)
  VAL(8,1,:)=-Y*Y1
  VAL(8,2,:)=Y1*T+Y*T-Y*Y1
  VAL(9,1,:)=-Y*(T1+T)
  VAL(9,2,:)=T*T1-Y*T1-Y*T
  VAL(4:9,:,:)=27/2D0*VAL(4:9,:,:)
  VAL(10,1,:)=27D0*Y*(T-X)
  VAL(10,2,:)=27D0*X*(T-Y)
END SUBROUTINE GradPhiOrder3

SUBROUTINE DerivPhi1DOrder3(X,DIM,VAL)
!!! DerivPhi1=(-11+36x-27x**2)/2
!!! DerivPhi2=9(1-5x+(9/2)x**2)
!!! DerivPhi3=9(4x-0.5-4.5x**2)
!!! DerivPhi4=(27x**2-18x+2)/2
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X
  real*8, dimension(4,DIM):: VAL
  !real*8 :: X  

  VAL(1,:)=(-11+36*X-27*X**2)/2
  VAL(2,:)=9*(1-5*X+9/2D0*X**2)
  VAL(3,:)=9*(4*X-1/2D0-9/2D0*X**2)
  VAL(4,:)=(27*X**2-18*X+2)/2
END SUBROUTINE DerivPhi1DOrder3

