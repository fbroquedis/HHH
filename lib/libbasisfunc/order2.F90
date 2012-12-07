!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Order 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Phi2DOrder2(X,Y,DIM,VAL)
!!! Phi1=2(1-x-y)(0.5-x-y)
!!! Phi2=2x(x-0.5)
!!! Phi3=2y(y-0.5)
!!! Phi4=4x(1-x-y)
!!! Phi5=4xy
!!! Phi6=4y(1-x-y)
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X,Y,T
  real*8, dimension(DIM):: X1,Y1,T1
  real*8, dimension(6,DIM):: VAL
  T=1D0-X-Y
  T1=0.5D0-X-Y
  X1=X-0.5D0
  Y1=Y-0.5D0  
  VAL(1,:)=2D0*T*T1
  VAL(2,:)=2D0*X*X1
  VAL(3,:)=2D0*Y*Y1
  VAL(4,:)=4D0*X*T
  VAL(5,:)=4D0*X*Y
  VAL(6,:)=4D0*Y*T
END SUBROUTINE Phi2DOrder2

SUBROUTINE Phi1DOrder2(X,DIM,VAL)
!!! Phi1=2(1-x)(1/2-x)
!!! Phi2=4x(1-x)
!!! Phi3=2x(x-1/2)
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X,X1,X2
  real*8, dimension(3,DIM):: VAL
  X1=X-0.5D0
  X2=1D0-X
  VAL(1,:)=-2D0*X2*X1
  VAL(2,:)=4D0*X*X2
  VAL(3,:)=2D0*X*X1
END SUBROUTINE Phi1DOrder2

SUBROUTINE GradPhi1DOrder2(X,DIM,VAL)
!!! Phi1=2(1-x)(1/2-x)
!!! Phi2=4x(1-x)
!!! Phi3=2x(x-1/2)
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X,X1,X2
  real*8, dimension(3,DIM):: VAL
  X1=X-0.5D0
  X2=1D0-X
  VAL(1,:)=-2D0*(X2-X1)
  VAL(2,:)=4D0*(X2-X)
  VAL(3,:)=2D0*(X1+X)
END SUBROUTINE gradPhi1DOrder2

SUBROUTINE GradPhiOrder2(X,Y,DIM,VAL)
!!! GradPhi1x=2*(2x+2y-1.5)
!!! GradPhi1y=2*(2x+2y-1.5)
!!! GradPhi2x=4x-1
!!! GradPhi2y=0
!!! GradPhi3x=0
!!! GradPhi3y=4y-1
!!! GradPhi4x=4(1-2x-y)
!!! GradPhi4y=-4x
!!! GradPhi5x=4y
!!! GradPhi5y=4x
!!! GradPhi6x=-4y
!!! GradPhi6y=4(1-x-2y)
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X
  real*8, dimension(DIM) ::Y
  real*8, dimension(6,2,DIM):: VAL
  VAL(1,1,:)=2D0*(2D0*X+2D0*Y-1.5D0)
  VAL(1,2,:)=2D0*(2D0*X+2D0*Y-1.5D0)
  VAL(2,1,:)=4D0*X-1D0
  VAL(2,2,:)=0.D0
  VAL(3,1,:)=0.D0
  VAL(3,2,:)=4D0*Y-1D0
  VAL(4,1,:)=4D0*(1D0-2D0*X-Y)
  VAL(4,2,:)=-4D0*X
  VAL(5,1,:)=4D0*Y
  VAL(5,2,:)=4D0*X
  VAL(6,1,:)=-4D0*Y
  VAL(6,2,:)=4D0*(1D0-X-2D0*Y)

END SUBROUTINE GradPhiOrder2


SUBROUTINE DerivPhi1DOrder2(X,DIM,VAL)
!!! DerivPhi1=4x-3
!!! DerivPhi2=4-8x
!!! DerivPhi3=4x-1
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X
  real*8, dimension(3,DIM):: VAL

  VAL(1,:)=4D0*X-3D0
  VAL(2,:)=4D0-8D0*X
  VAL(3,:)=4D0*X-1D0
END SUBROUTINE DerivPhi1DOrder2
