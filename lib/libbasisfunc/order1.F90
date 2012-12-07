!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Order 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Phi2DOrder1(X,Y,DIM,VAL)
!!! Phi1=1-x-y
!!! Phi2=x
!!! Phi3=y
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X,Y
  real*8, dimension(3,DIM):: VAL
  VAL(1,:)=(1.-X-Y)
  VAL(2,:)=X
  VAL(3,:)=Y
END SUBROUTINE Phi2DOrder1

SUBROUTINE Phi1DOrder1(X,DIM,VAL)
!!! Phi1=1-x
!!! Phi2=x
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X
  real*8, dimension(2,DIM):: VAL
  
  VAL(1,:)=1-X
  VAL(2,:)=X
END SUBROUTINE Phi1DOrder1

SUBROUTINE gradPhi1DOrder1(X,DIM,VAL)
!!! Phi1=1-x
!!! Phi2=x
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X
  real*8, dimension(2,DIM):: VAL
  
  VAL(1,:)=-1
  VAL(2,:)=1
END SUBROUTINE gradPhi1DOrder1

SUBROUTINE GradPhiOrder1(X,Y,DIM,VAL)
!!! GradPhi1x=-1
!!! GradPhi1y=-1
!!! GradPhi2x=1
!!! GradPhi2y=0
!!! GradPhi3x=0
!!! GradPhi3y=1
  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X
  real*8, dimension(DIM) ::Y
  real*8, dimension(3,2,DIM):: VAL
  VAL(1,1,:)=-1
  VAL(1,2,:)=-1
  VAL(2,1,:)=1
  VAL(2,2,:)=0
  VAL(3,1,:)=0
  VAL(3,2,:)=1

END SUBROUTINE GradPhiOrder1

SUBROUTINE DerivPhi1DOrder1(X,DIM,VAL)
!!!DerivPhi1D1=-1
!!!DerivPhi1D2=1

  IMPLICIT NONE 
  Integer :: DIM
  real*8, dimension(DIM):: X
  real*8, dimension(2,DIM):: VAL
!  real*8 :: X

  VAL(1,:)=-1
  VAL(2,:)=1
  
END SUBROUTINE DerivPhi1DOrder1



