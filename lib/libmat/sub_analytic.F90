subroutine calc_theta(x,y,angle)
	implicit none
	real*8 :: x,y,angle,pi,theta1,r
	
	pi=2.d0*asin(1.d0)
	if(x.lt.0) then
		
		angle=atan(y/x) + pi
		
	else if (x.eq.0) then
		if (y.gt.0) then
			angle=pi/2.
		else
			angle=3.*pi/2.
			
		end if
	else
		if(y.ge.0) then
			angle=atan(y/x)
		else
			angle=atan(y/x) + 2.*pi
		endif
	end if
	
	r=sqrt(x*x+y*y)
	theta1=2*atan((y/r)/(1+(x/r)))
	
! 	if (angle.ne.theta1) then
! 		write(6,*) 'calc theta pas bon',angle,theta1
! 	end if
! 	angle=theta1
end

subroutine analytic_Neumann(k,R0,R1,Nsomme,PhiEdge,type_cla,alphab)

  Use m_mat
  Use m_mesh
  Use m_gen
  implicit none
  real*8 :: somme,k,R0,R1,coor_x,coor_y,r,theta,V12(2),V23(2),V31(2),gamma13,gamma23
  integer :: i,l, Nsomme, NCALC,int_j,KK,type_cla
  real*8,dimension(:),allocatable :: J0,Y0,J1,Y1,Jr,Yr
  complex*16,dimension(:),allocatable :: Hankel0,Hankel1,coeffA,coeffB

  complex*16 :: tmp,c,d,dd,alphab,betab
  INTEGER :: PhiEdge(3,Order+1),Node(3)

  allocate(J0(Nsomme+2))
  allocate(Y0(Nsomme+2))
  allocate(Hankel0(Nsomme+2))
  allocate(J1(Nsomme+2))
  allocate(Y1(Nsomme+2))
  allocate(Hankel1(Nsomme+2))

  allocate(Jr(Nsomme+2))
  allocate(Yr(Nsomme+2))

  allocate(coeffA(Nsomme+1))
  allocate(coeffB(Nsomme+1))

  NCALC=0

  call RJBESL(k*R0,0.D0,Nsomme+2,J0,NCALC)	
  call RYBESL(k*R0,0.D0,Nsomme+2,Y0,NCALC)

  Hankel0 = cmplx(J0,Y0)



  call RJBESL(k*R1,0.D0,Nsomme+2,J1,NCALC)	
  call RYBESL(k*R1,0.D0,Nsomme+2,Y1,NCALC)

  Hankel1 = cmplx(J1,Y1)

  SELECT CASE(type_cla)
  CASE(1) ! Sommerfeld
     c=-(0,1)
  CASE(2) !Courbure
     c=1/(2*R1*k*R0)-(0,1)
  CASE(3) !OSRC
     c=-(0,1)*(k*R0+(0,1)*0.4*(k*R0)**(1.D0/3.D0)/R1**(2.D0/3.D0))/(k*R0)
  CASE(4) !propa+evanescents
     !   alphab=1D0
     dd=(alphab*(1/(2*R1)-(0,1)*k*R0)-(k*R0)**2)/(k*(-1/(2*R1)-(0,1)*k*R0+alphab)); 
     d=k*(1/(2*R1)-(0,1)*k*R0-alphab); 	
     !   dd=(alphab*(-1/(2*R1)-(0,1)*k*R0)+(k*R0)**2)/(k*(1/(2*R1)-(0,1)*k*R0-alphab)); 
     !   d=k*(-1/(2*R1)-(0,1)*k*R0+alphab); 	
  CASE(5)  !propa+rampants
     gamma13=2.6789385347d0
     gamma23=1.3541179394d0
     alphab=(6/R1)**(1/3D0)*gamma23/gamma13*(1/2D0-(0,1)*sqrt(3D0)/2D0)*(k*R0)**(2/3D0) 

     d=k*R0*(-1/(2*R1)-(0,1)*k*R0+alphab)
     dd=(alphab*(1/(2D0*R1)-(0,1)*k*R0)-(k*R0)**2)/d; 

  CASE(6) !Complete
     gamma13=2.6789385347d0
     gamma23=1.3541179394d0
     alphab=(6/R1)**(1D0/3D0)*gamma23/gamma13*(1/2-(0,1)*sqrt(3D0)/2D0)*(k*R0)**(2D0/3D0)
     betab=1D0 
  END SELECT






  do i=1,Nsomme+1
     if (type_cla.ge.4) then
        c=(dd+(i-1)*(i-1)/(R1)**2/d)
        if(type_cla.eq.6) then
           c=((3/R1-betab-alphab)/(4*R1**2)-(0,1)*k*R0/(2*R1**2)+(1/R1-betab-alphab)*(i-1)*(i-1)/R1**2)
           c=c/(1/(2*R1)*(3D0/R1-betab-alphab)+(0,1)*k*R0*(1/R1-betab-alphab)-(k*R0)**2+betab*alphab+(i-1)*(i-1)/R1**2);
           c=c+(0,1)*k*R0-1/(2*R1)
           c=-c/(k*R0)
        End if
     end if
     tmp=(-Hankel1(i+1)+((i-1)/(k*R1)+c)*Hankel1(i))&
          &/(-conjg(Hankel1(i+1))+((i-1)/(k*R1)+c)*conjg(Hankel1(i)))

     coeffA(i)=(-J0(i+1)+((i-1)/(k*R0))*J0(i))/&
          &(-Hankel0(i+1)+((i-1)/(k*R0))*Hankel0(i)-tmp*conjg(-Hankel0(i+1)+((i-1)/(k*R0))*Hankel0(i)))
     coeffB(i)=-coeffA(i)*tmp
  enddo

  !deallocate(J)
  !deallocate(Y)
  deallocate(Hankel0)
  deallocate(Hankel1)

  DO I=1,Ntri
     DO int_j=1,3
        coor_x = Coor(Tri(I,int_j),1)
        coor_y = Coor(Tri(I,int_j),2)
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)
        Uold((I-1)*Nphi+int_j)=coeffA(1)*cmplx(Jr(1),Yr(1))+coeffB(1)*cmplx(Jr(1),-Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+int_j)=Uold((I-1)*Nphi+int_j)+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l))+coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
        enddo
     END DO

     SELECT CASE(Order)

     CASE(2)
        Node = Tri(I,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)

        coor_x=Coor(Tri(I,1),1)+V12(1)/Order
        coor_y=Coor(Tri(I,1),2)+V12(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)


        Uold((I-1)*Nphi+PhiEdge(3,2))=coeffA(1)*cmplx(Jr(1),Yr(1))+coeffB(1)*cmplx(Jr(1),-Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(3,2))=Uold((I-1)*Nphi+PhiEdge(3,2))+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l))+coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
        enddo

        coor_x=Coor(Tri(I,2),1)+V23(1)/Order
        coor_y=Coor(Tri(I,2),2)+V23(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)
        Uold((I-1)*Nphi+PhiEdge(1,2))=coeffA(1)*cmplx(Jr(1),Yr(1))+coeffB(1)*cmplx(Jr(1),-Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(1,2))=Uold((I-1)*Nphi+PhiEdge(1,2))+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l))+coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
        enddo

        coor_x=Coor(Tri(I,3),1)+V31(1)/Order
        coor_y=Coor(Tri(I,3),2)+V31(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)
        Uold((I-1)*Nphi+PhiEdge(2,2))=coeffA(1)*cmplx(Jr(1),Yr(1))+coeffB(1)*cmplx(Jr(1),-Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(2,2))=Uold((I-1)*Nphi+PhiEdge(2,2))+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l))+coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
        enddo

     CASE(3)
        Node = Tri(I,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)
        DO KK=1,Order-1
           coor_x=Coor(Tri(I,1),1)+KK*V12(1)/Order
           coor_y=Coor(Tri(I,1),2)+KK*V12(2)/Order

        ENDDO
        DO KK=1,Order-1
           coor_x=Coor(Tri(I,2),1)+KK*V23(1)/Order
           coor_y=Coor(Tri(I,2),2)+KK*V23(2)/Order

        ENDDO
        DO KK=1,Order-1
           coor_x=Coor(Tri(I,3),1)+KK*V31(1)/Order
           coor_y=Coor(Tri(I,3),2)+KK*V31(2)/Order

        ENDDO
        coor_x=Coor(Tri(I,1),1)+V12(1)/3.-V31(1)/3.
        coor_y=Coor(Tri(I,1),2)+V12(2)/3.-V31(2)/3.


     end SELECT

  ENDDO

  deallocate(coeffA)
  deallocate(coeffB)
end subroutine analytic_Neumann
subroutine analytic_Neumann3(k,R0,R1,Nsomme,PhiEdge)

  Use m_mat
  Use m_mesh
  Use m_gen
  implicit none
  real*8 :: somme,k,R0,R1,coor_x,coor_y,r,theta,V12(2),V23(2),V31(2),gamma13,gamma23
  integer :: i,l, Nsomme, NCALC,int_j,KK
  real*8,dimension(:),allocatable :: J0,Y0,J1,Y1,Jr,Yr
  complex*16,dimension(:),allocatable :: Hankel0,Hankel1,coeffA,coeffB

  complex*16 :: c,d,dd
  complex*16 :: tmp,alphab
  INTEGER :: PhiEdge(3,Order+1),Node(3)

  allocate(J0(Nsomme+2))
  allocate(Y0(Nsomme+2))
  allocate(Hankel0(Nsomme+2))
  allocate(J1(Nsomme+2))
  allocate(Y1(Nsomme+2))
  allocate(Hankel1(Nsomme+2))

  allocate(Jr(Nsomme+2))
  allocate(Yr(Nsomme+2))

  allocate(coeffA(Nsomme+1))
  allocate(coeffB(Nsomme+1))

  NCALC=0

  call RJBESL(k*R0,0.D0,Nsomme+2,J0,NCALC)	
  call RYBESL(k*R0,0.D0,Nsomme+2,Y0,NCALC)

  Hankel0 = cmplx(J0,Y0)



  call RJBESL(k*R1,0.D0,Nsomme+2,J1,NCALC)	
  call RYBESL(k*R1,0.D0,Nsomme+2,Y1,NCALC)

  Hankel1 = cmplx(J1,Y1)	
  c=-(0,1)
  c=1/(2*R1*k*R0)-(0,1)
  ! c=-(0,1)*(k*R0+(0,1)*0.4*(k*R0)**(1.D0/3.D0)/R1**(2.D0/3.D0))/(k*R0)
  alphab=1D0
  !dd=(alphab*(1/(2*R1)-(0,1)*k*R0)-(k*R0)**2)/(k*R0*(-1/(2*R1)-(0,1)*k*R0+alphab)); 
  !d=R1**2*k*R0*(-1/(2*R1)-(0,1)*k*R0+alphab); 


  !		gamma13=2.6789385347d0
  !		gamma23=1.3541179394d0
  !         alphab=(6/R1)**(1/3D0)*gamma23/gamma13*(1/2D0-(0,1)*sqrt(3D0)/2D0)*(k*R0)**(2/3D0) 


  do i=1,Nsomme+1
     ! c=(alphab*(1/(2D0*R1)-(0,1)*k*R0)-(k*R0)**2+(i-1)*(i-1)/(R1)**2)/(k*R0*(-1/(2*R1)-(0,1)*k*R0+alphab)); 
     !c=(dd+(i-1)*(i-1)/d)
     tmp=(-Hankel1(i+1)+((i-1)/(R1)+c)*Hankel1(i))&
          &/(-conjg(Hankel1(i+1))+((i-1)/(R1)+c)*conjg(Hankel1(i)))

     coeffA(i)=(-J0(i+1)+((i-1)/(k*R0))*J0(i))/&
          &(-Hankel0(i+1)+((i-1)/(k*R0))*Hankel0(i)-tmp*conjg(-Hankel0(i+1)+((i-1)/(k*R0))*Hankel0(i)))
     coeffB(i)=-coeffA(i)*tmp
  enddo

  !deallocate(J)
  !deallocate(Y)
  deallocate(Hankel0)
  deallocate(Hankel1)

  DO I=1,NTri
     DO int_j=1,3
        coor_x = Coor(Tri(I,int_j),1)
        coor_y = Coor(Tri(I,int_j),2)
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)
        Uold((I-1)*Nphi+int_j)=coeffA(1)*cmplx(Jr(1),Yr(1))+coeffB(1)*cmplx(Jr(1),-Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+int_j)=Uold((I-1)*Nphi+int_j)+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l))+coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
        enddo
     END DO

     SELECT CASE(Order)

     CASE(2)
        Node = Tri(I,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)

        coor_x=Coor(Tri(I,1),1)+V12(1)/Order
        coor_y=Coor(Tri(I,1),2)+V12(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)


        Uold((I-1)*Nphi+PhiEdge(3,2))=coeffA(1)*cmplx(Jr(1),Yr(1))+coeffB(1)*cmplx(Jr(1),-Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(3,2))=Uold((I-1)*Nphi+PhiEdge(3,2))+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l))+coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
        enddo

        coor_x=Coor(Tri(I,2),1)+V23(1)/Order
        coor_y=Coor(Tri(I,2),2)+V23(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)
        Uold((I-1)*Nphi+PhiEdge(1,2))=coeffA(1)*cmplx(Jr(1),Yr(1))+coeffB(1)*cmplx(Jr(1),-Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(1,2))=Uold((I-1)*Nphi+PhiEdge(1,2))+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l))+coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
        enddo

        coor_x=Coor(Tri(I,3),1)+V31(1)/Order
        coor_y=Coor(Tri(I,3),2)+V31(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)
        Uold((I-1)*Nphi+PhiEdge(2,2))=coeffA(1)*cmplx(Jr(1),Yr(1))+coeffB(1)*cmplx(Jr(1),-Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(2,2))=Uold((I-1)*Nphi+PhiEdge(2,2))+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l))+coeffB(l)*cmplx(Jr(l),-Yr(l)))*dcos((l-1)*theta)
        enddo

     CASE(3)
        Node = Tri(I,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)
        DO KK=1,Order-1
           coor_x=Coor(Tri(I,1),1)+KK*V12(1)/Order
           coor_y=Coor(Tri(I,1),2)+KK*V12(2)/Order

        ENDDO
        DO KK=1,Order-1
           coor_x=Coor(Tri(I,2),1)+KK*V23(1)/Order
           coor_y=Coor(Tri(I,2),2)+KK*V23(2)/Order

        ENDDO
        DO KK=1,Order-1
           coor_x=Coor(Tri(I,3),1)+KK*V31(1)/Order
           coor_y=Coor(Tri(I,3),2)+KK*V31(2)/Order

        ENDDO
        coor_x=Coor(Tri(I,1),1)+V12(1)/3.-V31(1)/3.
        coor_y=Coor(Tri(I,1),2)+V12(2)/3.-V31(2)/3.


     end SELECT

  ENDDO

  deallocate(coeffA)
  deallocate(coeffB)
end subroutine analytic_Neumann3

subroutine sol_exacte_neumann(k,R0,Nsomme,PhiEdge,lhank)
  Use m_mat
  Use m_mesh
  Use m_gen
  implicit none
  real*8 :: somme,k,R0,coor_x,coor_y,r,theta,V12(2),V23(2),V31(2)
  integer :: i,l, Nsomme, NCALC,int_j,KK,lhank
  real*8,dimension(:),allocatable :: J0,Y0,Jr,Yr
  complex*16,dimension(:),allocatable :: Hankel0,coeffA,coeffB


  INTEGER :: PhiEdge(3,Order+1),Node(3)

  allocate(J0(Nsomme+2))
  allocate(Y0(Nsomme+2))
  allocate(Hankel0(Nsomme+2))


  allocate(Jr(Nsomme+2))
  allocate(Yr(Nsomme+2))

  allocate(coeffA(Nsomme+1))


  NCALC=0

  call RJBESL(k*R0,0.D0,Nsomme+2,J0,NCALC)	
  call RYBESL(k*R0,0.D0,Nsomme+2,Y0,NCALC)

  Hankel0 = cmplx(J0,Y0)

  do i=1,Nsomme+1
     coeffA(i)=(-J0(i+1)+((i-1)/(k*R0))*J0(i))/&
          &(-Hankel0(i+1)+((i-1)/(k*R0))*Hankel0(i))
  enddo
  !deallocate(J)
  !deallocate(Y)
  deallocate(Hankel0)


  DO I=1,NTri
     DO int_j=1,3
        coor_x = Coor(Tri(I,int_j),1)
        coor_y = Coor(Tri(I,int_j),2)
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)

        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)

        call calc_theta(coor_x,coor_y,theta)

        if(lhank.eq.1) then
           Uold((I-1)*Nphi+int_j)=coeffA(1)*cmplx(Jr(1),Yr(1))
        else
           !      do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+int_j)=2*cmplx(0,-1)**(lhank-1)*&
                &(coeffA(lhank)*cmplx(Jr(lhank),Yr(lhank)))*dcos((lhank-1)*theta)
           !	  enddo
        end if
     END DO

     SELECT CASE(Order)

     CASE(2)
        Node = Tri(I,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)

        coor_x=Coor(Tri(I,1),1)+V12(1)/Order
        coor_y=Coor(Tri(I,1),2)+V12(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)

        Uold((I-1)*Nphi+PhiEdge(3,2))=coeffA(1)*cmplx(Jr(1),Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(3,2))=Uold((I-1)*Nphi+PhiEdge(3,2))+2*cmplx(0,-1)**(l-1)*&
                &coeffA(l)*cmplx(Jr(l),Yr(l))*dcos((l-1)*theta)
        enddo

        coor_x=Coor(Tri(I,2),1)+V23(1)/Order
        coor_y=Coor(Tri(I,2),2)+V23(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)

        Uold((I-1)*Nphi+PhiEdge(1,2))=coeffA(1)*cmplx(Jr(1),Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(1,2))=Uold((I-1)*Nphi+PhiEdge(1,2))+2*cmplx(0,-1)**(l-1)*&
                &coeffA(l)*cmplx(Jr(l),Yr(l))*dcos((l-1)*theta)
        enddo

        coor_x=Coor(Tri(I,3),1)+V31(1)/Order
        coor_y=Coor(Tri(I,3),2)+V31(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)

        Uold((I-1)*Nphi+PhiEdge(2,2))=coeffA(1)*cmplx(Jr(1),Yr(1))

        do l=2,Nsomme+1     		
           Uold((I-1)*Nphi+PhiEdge(2,2))=Uold((I-1)*Nphi+PhiEdge(2,2))+2*cmplx(0,-1)**(l-1)*&
                &coeffA(l)*cmplx(Jr(l),Yr(l))*dcos((l-1)*theta)
        enddo
     end SELECT
  ENDDO

  deallocate(J0)
  deallocate(Y0)



  deallocate(Jr)
  deallocate(Yr)

  deallocate(coeffA)
end subroutine sol_exacte_neumann
subroutine sol_exacte_neumann2(k,R0,Nsomme,PhiEdge)
  Use m_mat
  Use m_mesh
  Use m_gen
  implicit none
  real*8 :: somme,k,R0,coor_x,coor_y,r,theta,V12(2),V23(2),V31(2)
  integer :: i,l, Nsomme, NCALC,int_j,KK,lhank
  real*8,dimension(:),allocatable :: J0,Y0,Jr,Yr
  complex*16,dimension(:),allocatable :: Hankel0,coeffA,coeffB


  INTEGER :: PhiEdge(3,Order+1),Node(3)

  allocate(J0(Nsomme+2))
  allocate(Y0(Nsomme+2))
  allocate(Hankel0(Nsomme+2))


  allocate(Jr(Nsomme+2))
  allocate(Yr(Nsomme+2))

  allocate(coeffA(Nsomme+1))


  NCALC=0

  call RJBESL(k*R0,0.D0,Nsomme+2,J0,NCALC)	
  call RYBESL(k*R0,0.D0,Nsomme+2,Y0,NCALC)

  Hankel0 = cmplx(J0,Y0)

  do i=1,Nsomme+1
     coeffA(i)=(-J0(i+1)+((i-1)/(k*R0))*J0(i))/&
          &(-Hankel0(i+1)+((i-1)/(k*R0))*Hankel0(i))
  enddo
  !deallocate(J)
  !deallocate(Y)
  deallocate(Hankel0)


  DO I=1,NTri
     DO int_j=1,3
        coor_x = Coor(Tri(I,int_j),1)
        coor_y = Coor(Tri(I,int_j),2)
        r=dsqrt(coor_x*coor_x+coor_y*coor_y)

        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)

        call calc_theta(coor_x,coor_y,theta)

        !      if(lhank.eq.1) then
        P_analytic((I-1)*Nphi+int_j)=coeffA(1)*cmplx(Jr(1),Yr(1))
        !      else
        do l=2,Nsomme+1     		
           P_analytic((I-1)*Nphi+int_j)=P_analytic((I-1)*Nphi+int_j)+2*cmplx(0,-1)**(l-1)*&
                &(coeffA(l)*cmplx(Jr(l),Yr(l)))*dcos((l-1)*theta)
        enddo
        !   end if
     END DO

     SELECT CASE(Order)

     CASE(2)
        Node = Tri(I,:)
        V12(1)=Coor(Node(2),1)-Coor(Node(1),1)
        V12(2)=Coor(Node(2),2)-Coor(Node(1),2)
        V23(1)=Coor(Node(3),1)-Coor(Node(2),1)
        V23(2)=Coor(Node(3),2)-Coor(Node(2),2)
        V31(1)=Coor(Node(1),1)-Coor(Node(3),1)
        V31(2)=Coor(Node(1),2)-Coor(Node(3),2)

        coor_x=Coor(Tri(I,1),1)+V12(1)/Order
        coor_y=Coor(Tri(I,1),2)+V12(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)

        P_analytic((I-1)*Nphi+PhiEdge(3,2))=coeffA(1)*cmplx(Jr(1),Yr(1))

        do l=2,Nsomme+1     		
           P_analytic((I-1)*Nphi+PhiEdge(3,2))=P_analytic((I-1)*Nphi+PhiEdge(3,2))+2*cmplx(0,-1)**(l-1)*&
                &coeffA(l)*cmplx(Jr(l),Yr(l))*dcos((l-1)*theta)
        enddo

        coor_x=Coor(Tri(I,2),1)+V23(1)/Order
        coor_y=Coor(Tri(I,2),2)+V23(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)

        P_analytic((I-1)*Nphi+PhiEdge(1,2))=coeffA(1)*cmplx(Jr(1),Yr(1))

        do l=2,Nsomme+1     		
           P_analytic((I-1)*Nphi+PhiEdge(1,2))=P_analytic((I-1)*Nphi+PhiEdge(1,2))+2*cmplx(0,-1)**(l-1)*&
                &coeffA(l)*cmplx(Jr(l),Yr(l))*dcos((l-1)*theta)
        enddo

        coor_x=Coor(Tri(I,3),1)+V31(1)/Order
        coor_y=Coor(Tri(I,3),2)+V31(2)/Order

        r=dsqrt(coor_x*coor_x+coor_y*coor_y)
        call RJBESL(k*r,0.D0,Nsomme+1,Jr,NCALC)	
        call RYBESL(k*r,0.D0,Nsomme+1,Yr,NCALC)
        call calc_theta(coor_x,coor_y,theta)

        P_analytic((I-1)*Nphi+PhiEdge(2,2))=coeffA(1)*cmplx(Jr(1),Yr(1))

        do l=2,Nsomme+1     		
           P_analytic((I-1)*Nphi+PhiEdge(2,2))=P_analytic((I-1)*Nphi+PhiEdge(2,2))+2*cmplx(0,-1)**(l-1)*&
                &coeffA(l)*cmplx(Jr(l),Yr(l))*dcos((l-1)*theta)
        enddo
     end SELECT
  ENDDO

  deallocate(J0)
  deallocate(Y0)



  deallocate(Jr)
  deallocate(Yr)

  deallocate(coeffA)
end subroutine sol_exacte_neumann2

