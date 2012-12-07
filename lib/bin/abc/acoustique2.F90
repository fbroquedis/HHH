program acoustique
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
! En acoustique, pense à déclarer A réel (dans mod/m_mat.f90)
! Attention aussi a la courbure
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
  Use m_condbord    
  Use m_mat
  Use m_gen
  Use m_mesh
  Use m_output
  Use m_source
  Use m_inter


  implicit none 
  character*100 :: nommaillage
 ! character*100 :: nam
  integer :: I,J,toto,K,N,INFO,PhiEdge(3,2),II,taille,III,tmp,nmax,forme_obstacle
!  real*8,dimension (:,:,:),allocatable :: Araff,Anoraff,Bcla2inverse
  integer, dimension (:,:),allocatable ::IPIV,IPIV2
  real*8, dimension(:,:), allocatable,save :: F,Fcomp
  integer :: ntfinal,recep,obstacle,Ineigh,Arete,l,int_j,periodic
integer, parameter ::Nphi1=1,Nphi2=6,Nphi3=10
  real*8 :: pi,x1,y1,x,y,x2,y2,x3,y3,demi_grand_axe_a,demi_petit_axe_b,rayon_cercle
  real*8 :: courbure1,courbure2,aire_arete,theta0
  real*8 :: x_centre,y_centre,d1,d2,d3,dtmp,theta_snap
  real*8,dimension (:),allocatable :: Ftest,Ftest1
!integer,dimension (:),allocatable :: PhiEdge
!  integer, dimension(:),allocatable :: num_rcv
REAL*8,allocatable :: Phi(:,:)
	REAL*8,allocatable :: GradPhi1D(:,:,:,:)	
!!!Pour un milieu périodique :::
periodic=0
helmholtz=0
  pi=acos(-1.0)

!allocate(PhiEdge(3,2))
     PhiEdge(1,1)=2
     PhiEdge(1,2)=3
     PhiEdge(2,1)=3
     PhiEdge(2,2)=1
     PhiEdge(3,1)=1
     PhiEdge(3,2)=2

!!!!!!!!!!!!!!!!!!!!!!
!!
!!   TEMPS FINAL
!!
!!!!!!!!!!!!!!!!!!!!!!
 ! ntfinal=5000
  ntfinal=6000

xgrid1=-1.5D0
xgrid2=1.5D0
ygrid1=-1.5D0
ygrid2=1.5D0
stepx=0.005
stepy=0.005

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   LECTURE DES DONNEES : MAILLAGE, ORDRE DES ELEMENTS, ...
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        nom du maillage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(6,*) 'Name of mesh files?'
  read(5,*) nommaillage
  J=len_trim(nommaillage)
  write(6,*) trim(nommaillage)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      ordre des éléments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(6,*) 'Order of the elements'
  read(5,*) Order
!! Nphi : Nb de ddl par element. 
!! alpha : parametre de penalisation
  SELECT CASE(Order)
  CASE(1)
     Nphi=3
     alpha=3
  CASE(2)
     Nphi=6
     alpha=5
  CASE(3)
      Nphi=10
     alpha=11
     alpha2=6
  CASE DEFAULT
     write(6,*) 'Not implemented'
     stop
  END SELECT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! type de conditions de bord
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(6,*) 'Conditions autres que Neumann sur les bords?'
  read(5,*) CondBord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     type de domaine : 
!!     cercle ou donut ? 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! write(6,*) 'Rayon du cercle ?'
 ! read(5,*) rayon
  write(6,*) 'Y-a-t''il un obstacle ?'
  read(5,*) obstacle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! source ou onde incidente
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (obstacle==1) then
        write(6,*) 'Omega ?'
	read(5,*) omega_incident
        write(6,*)omega_incident
        write(6,*) 'theta incident?'
        read(5,*) theta0
        write(6,*) theta0
    else
	write(6,*) 'Center of Initial Condition'
	read(5,*) x0,y0
	write(6,*) x0,y0	
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    Lecture du maillage
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(6,*) 'Reading file ',  nommaillage(1:J)//".1.node"
  open(20,FILE=trim(nommaillage)//".1.node")
  read(20,*) NPoints
  write(6,*) NPoints, ' nodes'


  allocate(Coor(Npoints,2))
!!! On lit les coordonnées de chaque sommet et eventuellement la condition de bord
  IF (CondBord==0) then
     iscla=0
     DO I=1,NPoints
        READ(20,*) toto,Coor(I,1),Coor(I,2)
     ENDDO
  else
     allocate(CondPoint(Npoints))
     DO I=1,NPoints
        READ(20,*) toto,Coor(I,1),Coor(I,2),CondPoint(I)
     ENDDO
     IF (maxval(CondPoint(1:Npoints))==3) then
        iscla=1
     else
        iscla=0
     end IF
  END IF
  close(20)

  write(6,*) 'Reading file ',  nommaillage(1:J)//".1.ele"
  open(20,FILE=trim(nommaillage)//".1.ele")
  read(20,*) NTri
  write(6,*) NTri, ' triangles'
  if(periodic.eq.1) then
     Ntri=2*Ntri
  end if
  allocate(Tri(NTri,3))
  allocate(Type_Media(NTri))
!!! On lit les trois sommets de chaque triangle et sa vitesse
  DO I=1,NTri
     READ(20,*) toto,Tri(I,1),Tri(I,2),Tri(I,3),Type_Media(I)
  ENDDO
  close(20)
  write(6,*) 'Reading file ',  nommaillage(1:J)//".1.neigh"
  open(20,FILE=trim(nommaillage)//".1.neigh")
  read(20,*) NNeigh
  if (NNeigh.ne.NTri) then
     write(6,*) 'Error : NNeigh is different from NTri'
     stop
  endif
     allocate(Neigh(NTri,3))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!       CONSTRUCTION DES VOISINS A L'INTERIEUR
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(6,*)'debut construction voisin'
!On lit les voisins de chaque triangles (le voisin i est opposé au sommet i)
  DO I=1,NTri
    READ(20,*) toto,Neigh(I,1),Neigh(I,2),Neigh(I,3)
!! How to read Neigh : triangle Neigh(I,1) is a neighbourg of triangle I. They are connected by edge Tri(I,2),Tri(I,3).
  ! a-t-on autre chose que Neumann sur le bord
    IF (Condbord==1) then
       DO J=1,3
	!le triangle I a-t-il un voisin connecté par l'arête J
         IF (Neigh(I,J)==-1) then
             If (CondPoint(Tri(I,PhiEdge(J,1))).gt.0) then           
     Neigh(I,J)=-CondPoint(Tri(I,PhiEdge(J,1)))
             elseif (CondPoint(Tri(I,PhiEdge(J,2))).gt.0) then
                Neigh(I,J)=-CondPoint(Tri(I,PhiEdge(J,2)))
             END If
          END IF
       END DO
    END IF
  END DO
call sub_renum3
call my_sub_param_phys
!!!!!!!!!!!!!!!!!!!!!
!!
!!    SISMOGRAMMES
!!
!!!!!!!!!!!!!!!!!!!!!
     write(6,*) 'Combien de sismogrammes voulez-vous?'
     read(5,*) nbrcv
     write(6,*) nbrcv
! Cas où on donne toutes les coordonnées
     if(nbrcv.ge.1) then
        allocate(coord_rcv(nbrcv,2))
	recep=1
	DO i=1,nbrcv
          write(6,*)'coordonnées du récepteur ', i ,' :'
          read(5,*)x_rcv,y_rcv
          write(6,*) x_rcv,y_rcv
          coord_rcv(i,1)=x_rcv
          coord_rcv(i,2)=y_rcv
	END DO
     endif


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!           CONSTRUCTION DES ARETES DU BORD
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!On libère un peu de mémoire
    IF (Condbord==1) then
  deallocate(Condpoint)
  endif
  close(20)

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    APPEL DES MATRICES
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Call sub_allocate
!!!Définition des formules d'interpolation 2D et 3D (noeuds et poids sur l'élément de référence
  CALL sub_defGL2D
  CALL sub_defGL1D

!!!On construit la matrice de masse et celle de raideur
  CALL sub_defmassmat


  Call sub_defstiffmat
 
  if (iscla==1) then
    ! On calcule les matrices associées aux CLA
    call sub_CLA
   !  call sub_defmatBBCLA_bis
  endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    CALCUL DU dt EN FONCTION DE LA CFL
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  write(6,*) CFL
  
  SELECT CASE(Order)
  CASE(1)
    dt=cfl*0.125
  CASE(2)
    dt=cfl*.156
  CASE(3)
!    dt=cfl*0.05
!    dt=cfl*0.048
    dt=cfl*0.048*.5
 !  dt=0.066375
  END SELECT
!dt=3.00000000000000006E-003
write(6,*) cfl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    Source ou onde incidente
!! Pour l'instant on ne regarde que le deuxième cas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (obstacle==1) then
	!! onde incidente
!        write(6,*) 'hello'
	call sub_rcv_localisation

!        write(6,*) 'hello 2'
    else
	!! source
	call sub_source_localisation
	call sub_val_source
	call sub_rcv_localisation
	call sub_centre
ValPhisource=Matmul(Minv,ValPhisource)/DFVec(Isource)
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ON COMPTE LES INSTANTS POUR NE PAS TOUT REPRESENTER
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NBINST=0
NBINST2=0
nbsis=0
write(6,*) CFL,dt
write(6,*) 'dt',dt




!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     ON STOCKE LES SISMOS, ...
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!fichier qui stocke le sismo
 open(21,File="FILM/sismos_.dat")
 open(22,File="FILM/energ_cinetik.dat")
!write(6,*) 'test7'
if (obstacle==1)then
allocate(F(Ntri*Nphi,1))
allocate(Fcomp(Ntri*Nphi,1))
F=0.D0
Fcomp=0.D0
write(6,*) 'toto1'
 call sub_rhsacoustique(F,Fcomp,omega_incident,theta0)
write(6,*) 'toto2'
else
allocate(F(Nphi,1))
endif
!fichier qui stocke les valeurs de phi1
!open(22,File="FILM/phi1.dat")

!saving the mesh in FILM/mesh.vtu
! if (iscla==1) then
! call sub_writemeshCLA
!else
 call sub_writemesh
!endif

 call sub_write_pas_de_temps(dt)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    ALLOCATIONS
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! U
if (obstacle==1) then
allocate(U_inc(Nphi*Ntri))
U_inc=0.D0
else
F=0.D0
endif
Uold=0.D0
taille=size(Uold)
theta0 = 0.d0

!!les vecteurs AU, ...
allocate(AUold(Nphi*Ntri))
allocate(AUinter(Nphi*Ntri))
allocate(MUnew(taille))
allocate(MUold(taille))
allocate(MUinter(Nphi*Ntri))
allocate(CUold(taille))


!!vecteur contenant U et W
allocate(Ftest(Nphi))
allocate(Ftest1(Nphi))


write(6,*) 'COUCOU C MOI'
!write(6,*) Dclabis(:,:,:)
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    LA BOUCLE EN TEMPS
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(6,*) Narete_du_bord

AUold=0.D0
AUinter=0.D0
MUnew=0.D0
MUold=0.D0
MUinter=0.D0
Cuold=0.D0

  if (iscla==1) then
allocate(IPIV(Nbcla,Nphi))
! On calcule et on factorise I+dt/2*Bcla par la factorisation LU
     Bcla2=0.D0      
     DO I=1,NbCla
        DO J=1,Nphi
          Bcla2(I,J,J)=1.D0
        END DO
        !! Attention, la j'ai mis la courbure
!        Bcla2(I,:,:)=Bcla2(I,:,:)+(dt/2.+dt**2/4d0/2D0)*Bcla(I,:,:)
        Bcla2(I,:,:)=Bcla2(I,:,:)+dt/2.*Bcla(I,:,:)/sqrt(rho(TriCla(I))*mu(TriCla(I)))
        call DGETRF( Nphi, Nphi,Bcla2(I,:,:), Nphi, IPIV(I,:),INFO )
     ENDDO
  end if

!! U=0.D0
!! AU=0.D0
 !OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(J)
    DO I=1,Ntri
      AU((I-1)*Nphi+1:I*Nphi)=0.D0
      U((I-1)*Nphi+1:I*Nphi)=0.D0
    END DO
!OMP END PARALLEL DO


!!! La boucle en temps (enfin)
DO N=0,ntfinal
!     write(6,*) N

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!
     !!       TRIANGLE A L'INTERIEUR
     !!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$$$     SELECT CASE(ORDER)
!$$$     CASE(1)
!$$$  !$OMP PARALLEL DO    DEFAULT(SHARED) PRIVATE(I,J)    
!$$$     DO I=1,Ntri
!$$$        AU((I-1)*Nphi1+1:I*Nphi1)=-MATMUL(A(I,1:Nphi1,1:Nphi1),U((I-1)*Nphi1+1:I*Nphi1))
!$$$
!$$$        DO J=1,3
!$$$           IF (Neigh(I,J).gt.0) then
!$$$              AU((I-1)*Nphi1+1:I*Nphi1)=AU((I-1)*Nphi1+1:I*Nphi1)-MATMUL(A(I,1:Nphi1,J*Nphi1+1:(J+1)*Nphi1)&
!$$$                   &,U((Neigh(I,J)-1)*Nphi1+1:Neigh(I,J)*Nphi1))
!$$$           end IF
!$$$        END DO
!$$$!!!!!!!!!!!!!
!$$$      END DO  
!$$$ !$OMP END PARALLEL DO
!$$$
!$$$CASE(2)
!$$$  !$OMP PARALLEL DO    DEFAULT(SHARED) PRIVATE(I,J)    
!$$$     DO I=1,Ntri
!$$$        AU((I-1)*Nphi2+1:I*Nphi2)=-MATMUL(A(I,1:Nphi2,1:Nphi2),U((I-1)*Nphi2+1:I*Nphi2))
!$$$
!$$$        DO J=1,3
!$$$           IF (Neigh(I,J).gt.0) then
!$$$              AU((I-1)*Nphi2+1:I*Nphi2)=AU((I-1)*Nphi2+1:I*Nphi2)-MATMUL(A(I,1:Nphi2,J*Nphi2+1:(J+1)*Nphi2)&
!$$$                   &,U((Neigh(I,J)-1)*Nphi2+1:Neigh(I,J)*Nphi2))
!$$$           end IF
!$$$        END DO
!$$$!!!!!!!!!!!!!
!$$$      END DO  
!$$$ !$OMP END PARALLEL DO
!$$$
!$$$CASE(3)
  !$OMP PARALLEL DO    DEFAULT(SHARED) PRIVATE(I,J)    
     DO I=1,Ntri
        AU((I-1)*Nphi3+1:I*Nphi3)=-MATMUL(A(I,1:Nphi3,1:Nphi3),U((I-1)*Nphi3+1:I*Nphi3))

        DO J=1,3
           IF (Neigh(I,J).gt.0) then
              AU((I-1)*Nphi3+1:I*Nphi3)=AU((I-1)*Nphi3+1:I*Nphi3)-MATMUL(A(I,1:Nphi3,J*Nphi3+1:(J+1)*Nphi3)&
                   &,U((Neigh(I,J)-1)*Nphi3+1:Neigh(I,J)*Nphi3))
           end IF
        END DO
!!!!!!!!!!!!!
      END DO  
 !$OMP END PARALLEL DO

!$$$   END SELECT

  !$OMP PARALLEL DO    DEFAULT(SHARED) PRIVATE(I)    
     DO I=1,Ntri*Nphi
      U(i)=dt**2*AU(i)+2*U(i)-Uold(i)
   END DO
  !$OMP END PARALLEL DO

      F(:,1)=ValPhisource(:,1)*Valsource(N+1)
      U((Isource-1)*Nphi+1:Isource*Nphi)=U((Isource-1)*Nphi+1:Isource*Nphi)+dt**2*F(:,1)
!!! Absorbing boundary conditions
!$$$  if (iscla==1) then
!$$$  !$OMP PARALLEL DO    DEFAULT(SHARED) PRIVATE(I,Ftest,Info,J)    
!$$$     DO I=1,NbCla
!$$$        J=TriCla(I)
!$$$        Ftest=U((J-1)*Nphi+1:J*Nphi)+dt/2.D0*matmul(Bcla(I,1:Nphi,1:Nphi),Uold((J-1)*Nphi+1:J*Nphi))/sqrt(rho(J)*mu(J))
!$$$        !! On divise par (I+dt/2*C)
!$$$       CALL DGETRS( 'N',Nphi, 1,Bcla2(I,:,:), Nphi, IPIV(I,:),Ftest,Nphi,INFO )
!$$$        U((J-1)*Nphi+1:J*Nphi)=Ftest
!$$$     ENDDO
!$$$ !$OMP END PARALLEL DO
!$$$   end if
flush(22)

  Uold=Uinter
  Auold=AUinter
  AUinter=AU
  MUold=MUinter
  MUinter=MUnew
  Uinter=U

IF(Mod(N,200).eq.0) then 
  NBINST=NBINST+1
  call sub_writeunstruct(N)
END IF
END DO

 close(21)  
 close(22)
end program acoustique
