program acoustique_scalaire_vectoriel
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
  real*8,dimension (:,:,:),allocatable :: Bflusol
  integer, dimension (:,:),allocatable ::IPIV,IPIV2,IPIVflusol
  real*8, dimension(:,:), allocatable,save :: F,Fcomp
  integer :: ntfinal,recep,obstacle,Ineigh,Arete,Arete2,l,int_j,periodic
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
  ntfinal=600000

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
     alpha2=3
  CASE(2)
     Nphi=6
     alpha=5
  CASE(3)
      Nphi=10
     alpha=10
     alpha2=10
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
if(periodic.eq.1) then
 NPoints=2*Npoints
  allocate(Coor(Npoints,2))
  IF (CondBord==0) then
     iscla=0
     DO I=1,NPoints/2
        READ(20,*) toto,Coor(I,1),Coor(I,2)
        Coor(I+Npoints/2,1)=1D0-Coor(I,1)
        Coor(I+Npoints/2,2)=Coor(I,2)
     ENDDO
  else
     allocate(CondPoint(Npoints))
     DO I=1,NPoints/2
        READ(20,*) toto,Coor(I,1),Coor(I,2),CondPoint(I)
        Coor(I+Npoints/2,1)=1D0-Coor(I,1)
        Coor(I+Npoints/2,2)=Coor(I,2)
        CondPoint(I+Npoints/2)=CondPoint(I)
     ENDDO
     IF (maxval(CondPoint(1:Npoints))==3) then
        iscla=1
   else
        iscla=0
     end IF
  END IF
else

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
end if
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
if (periodic.eq.1) then
  DO I=1,NTri/2
     READ(20,*) toto,Tri(I,1),Tri(I,2),Tri(I,3),TYPE_MEDIA(I)
     Tri(Ntri+1-I,:)=Tri(I,3:1:-1)+Npoints/2
     TYPE_MEDIA(Ntri+1-I)=TYPE_MEDIA(I)
  ENDDO
else
  DO I=1,NTri
     READ(20,*) toto,Tri(I,1),Tri(I,2),Tri(I,3),Type_Media(I)
  ENDDO
end if
  close(20)
  write(6,*) 'Reading file ',  nommaillage(1:J)//".1.neigh"
  open(20,FILE=trim(nommaillage)//".1.neigh")
  read(20,*) NNeigh
if (periodic.eq.1) then
 NNeigh=2*NNeigh
end if
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
if (periodic.eq.1) then
  DO I=1,NTri/2
    READ(20,*) toto,Neigh(I,1),Neigh(I,2),Neigh(I,3)
!!$!! How to read Neigh : triangle Neigh(I,1) is a neighbourg of triangle I. Th\
ey are connected by edge Tri(I,2),Tri(I,3).
 DO J=1,3
    If ((Neigh(I,J).eq.-1)) then
       if ((coor(Tri(I,PhiEdge(J,1)),1).lt.1e-6)&
            &.and.(coor(Tri(I,PhiEdge(J,2)),1).lt.1e-6)) then
          Neigh(I,J)=Ntri+1-I
          Neigh(Ntri+1-I,4-J)=I
       elseif ((coor(Tri(I,PhiEdge(J,1)),1).gt..5-1e-6)&
           &.and.(coor(Tri(I,PhiEdge(J,2)),1).gt..5-1e-6)) then
          Neigh(I,J)=Ntri+1-I
          Neigh(Ntri+1-I,4-J)=I
       elseif(condbord.eq.1) then
             If (CondPoint(Tri(I,PhiEdge(J,1))).gt.0) then           
     Neigh(I,J)=-CondPoint(Tri(I,PhiEdge(J,1)))
          Neigh(Ntri+1-I,4-J)=-CondPoint(Tri(I,PhiEdge(J,1)))
             elseif (CondPoint(Tri(I,PhiEdge(J,2))).gt.0) then
                Neigh(I,J)=-CondPoint(Tri(I,PhiEdge(J,2)))
          Neigh(Ntri+1-I,4-J)=-CondPoint(Tri(I,PhiEdge(J,2)))
             END If
       end If
    else
          Neigh(Ntri+1-I,4-J)=Ntri+1-Neigh(I,J)
    end If
 end DO
end DO
else
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
end if
call sub_renum2
call sub_param_phys
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

  Call sub_allocate_scalaire_vectoriel
!!!Définition des formules d'interpolation 2D et 3D (noeuds et poids sur l'élément de référence
  CALL sub_defGL2D
  CALL sub_defGL1D

!!!On construit la matrice de masse et celle de raideur
  CALL sub_defmassmat

  CFL=10**8
Call sub_defstiffmat_elasto_acoustic_pp
write(6,*) 'CFL',CFL
Call sub_defstiffmat_elasto_acoustic_uxux
write(6,*) 'CFL',CFL
Call sub_defstiffmat_elasto_acoustic_uxuy
write(6,*) 'CFL',CFL
Call sub_defstiffmat_elasto_acoustic_uyux
write(6,*) 'CFL',CFL
Call sub_defstiffmat_elasto_acoustic_uyuy
write(6,*) 'CFL',CFL
Call sub_defstiffmat_elasto_acoustic_pux
write(6,*) 'CFL',CFL
Call sub_defstiffmat_elasto_acoustic_puy
write(6,*) 'CFL',CFL
Call sub_defstiffmat_elasto_acoustic_uxp
write(6,*) 'CFL',CFL
Call sub_defstiffmat_elasto_acoustic_uyp
write(6,*) 'CFL',CFL

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
    dt=cfl*0.125/10D0
  CASE(2)
    dt=cfl*156
  CASE(3)
!    dt=cfl*0.05
!    dt=cfl*0.048
    dt=cfl*0.048*0.5
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
	call sub_source_elasto_acoustic_localisation
	call sub_val_source
	call sub_rcv_localisation
	call sub_centre
if(Isource.le.Nflu+Nflusol) then
ValPhisource=Matmul(Minv,ValPhisource)/DFVec(Isource)
else
ValPhisource(:,1)=Matmul(Minv,ValPhisource(:,1))/DFVec(Isource)
ValPhisource(:,2)=Matmul(Minv,ValPhisource(:,2))/DFVec(Isource)
endif
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
 open(23,File="FILM/sismos_P.dat")
 open(24,File="FILM/sismos_Ux.dat")
 open(25,File="FILM/sismos_Uy.dat")
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
 !call sub_writemeshacouselasto
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
Pold=0.D0
P=0.D0
Uxold=0.D0
Ux=0.D0
Uyold=0.D0
Uy=0.D0
taille=size(Pold)
theta0 = 0.d0

!!les vecteurs AU, ...
allocate(APold(Nphi*(Nflu+Nflusol)))
allocate(APinter(Nphi*(Nflu+Nflusol)))
allocate(AUxold(Nphi*(Nflu+Nflusol)))
allocate(AUxinter(Nphi*(Nflu+Nflusol)))
allocate(AUyold(Nphi*(Nflu+Nflusol)))
allocate(AUyinter(Nphi*(Nflu+Nflusol)))
allocate(MPnew(taille))
allocate(MPold(taille))
allocate(MPinter(Nphi*(Nflu+Nflusol)))
allocate(MUxnew(taille))
allocate(MUxold(taille))
allocate(MUxinter(Nphi*(Nflu+Nflusol)))
allocate(MUynew(taille))
allocate(MUyold(taille))
allocate(MUyinter(Nphi*(Nflu+Nflusol)))
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
AP=0.D0
APold=0.D0
APinter=0.D0
AUx=0.D0
AUxold=0.D0
AUxinter=0.D0
AUy=0.D0
AUyold=0.D0
AUyinter=0.D0
MPnew=0.D0
MPold=0.D0
MPinter=0.D0
MUxnew=0.D0
MUxold=0.D0
MUxinter=0.D0
MUynew=0.D0
MUyold=0.D0
MUyinter=0.D0
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
write(6,*) 'toto'
        call DGETRF( Nphi, Nphi,Bcla2(I,:,:), Nphi, IPIV(I,:),INFO )
     ENDDO
  end if
allocate(IPIVflusol(Nflusol,Nphi))
allocate(Bflusol(Nflusol,Nphi,Nphi))
! On calcule et on factorise I+dt/2*Bcla par la factorisation LU
     Bflusol=0.D0      
     DO I=1,Nflusol
        DO J=1,Nphi
          Bflusol(I,J,J)=1.D0
        END DO
Arete=Corres_flusol(I)
K=Neigh(I+Nflu,Arete)-Nflu-Nflusol
Arete2=Corres_solflu(K)

Bflusol(I,:,:)=Bflusol(I,:,:)-dt**2/4D0&
     &*(matmul(A_pux(I,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),&
     &A_uxp(K,1:Nphi,(Arete2-1)*Nphi+1:Arete2*Nphi))&
     &+matmul(A_puy(I,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),&
     &A_uyp(K,1:Nphi,(Arete2-1)*Nphi+1:Arete2*Nphi)))
!!$Bflusol(I,:,:)=Bflusol(I,:,:)-dt/2D0*(matmul(A_pux(I,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),&
!!$     &A_uxp(K,1:Nphi,(Arete2-1)*Nphi+1:Arete2*Nphi))&
!!$     &+matmul(A_puy(I,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),A_uyp(K,1:Nphi,(Arete2-1)*Nphi+1:Arete2*Nphi)))
        call DGETRF( Nphi, Nphi,Bflusol(I,:,:), Nphi, IPIVflusol(I,:),INFO )
     ENDDO
!!! La boucle en temps (enfin)
DO N=0,ntfinal
     write(6,*) N

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!
     !!       TRIANGLE A L'INTERIEUR
     !!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               
     DO I=1,Nflu+Nflusol
        AP((I-1)*Nphi+1:I*Nphi)=-MATMUL(A_pp(I,1:Nphi,1:Nphi),P((I-1)*Nphi+1:I*Nphi))

        DO J=1,3
           IF ((Neigh(I,J).gt.0).and.(Neigh(I,J).le.Nflu+Nflusol)) then
 
               AP((I-1)*Nphi+1:I*Nphi)=AP((I-1)*Nphi+1:I*Nphi)-MATMUL(A_pp(I,1:Nphi,J*Nphi+1:(J+1)*Nphi)&
                   &,P((Neigh(I,J)-1)*Nphi+1:Neigh(I,J)*Nphi))
           end IF
        END DO
!!!! Les deux lignes suivantes servent juste pour le calcul de l'énergi
        MPnew((I-1)*Nphi+1:I*Nphi)=MATMUL(M(1:Nphi,1:Nphi),P((I-1)*Nphi+1:I*Nphi))*DFvec(I)
        MPold((I-1)*Nphi+1:I*Nphi)=MATMUL(M(1:Nphi,1:Nphi),Pold((I-1)*Nphi+1:I*Nphi))*DFvec(I)
!!!!!!!!!!!!!
      END DO  

      P(1:Nflu*Nphi)=dt**2*AP(1:Nflu*Nphi)+2*P(1:Nflu*Nphi)-Pold(1:Nflu*Nphi)

     
if (obstacle==1) then
!	U=U+dt**2*(F(:,1)*sin(omega_incident*N*dt)+0.D0*Fcomp(:,1)*sin(omega_incident*N*dt))
else
if(Isource.le.Nflu+Nflusol) then
        F(:,1)=ValPhisource(:,1)*Valsource(N+1)
	P((Isource-1)*Nphi+1:Isource*Nphi)=P((Isource-1)*Nphi+1:Isource*Nphi)+dt**2*F(:,1)
end if
endif
!     write(22,*) sum(((MPnew-MPinter)/dt)*((P-Pinter)/dt))-sum(AP*Pinter)

     DO I=1,Nsol+Nsolflu
        AUx((I-1)*Nphi+1:I*Nphi)=-MATMUL(A_uxux(I,1:Nphi,1:Nphi),Ux((I-1)*Nphi+1:I*Nphi))
        AUx((I-1)*Nphi+1:I*Nphi)=AUx((I-1)*Nphi+1:I*Nphi)-MATMUL(A_uxuy(I,1:Nphi,1:Nphi),Uy((I-1)*Nphi+1:I*Nphi))
        AUy((I-1)*Nphi+1:I*Nphi)=-MATMUL(A_uyuy(I,1:Nphi,1:Nphi),Uy((I-1)*Nphi+1:I*Nphi))
        AUy((I-1)*Nphi+1:I*Nphi)=AUy((I-1)*Nphi+1:I*Nphi)-MATMUL(A_uyux(I,1:Nphi,1:Nphi),Ux((I-1)*Nphi+1:I*Nphi))

        DO J=1,3
           IF ((Neigh(I+Nflu+Nflusol,J).gt.0).and.(Neigh(I+Nflu+Nflusol,J).gt.Nflu+Nflusol)) then
              K=Neigh(I+Nflu+Nflusol,J)-Nflu-Nflusol
              AUx((I-1)*Nphi+1:I*Nphi)=AUx((I-1)*Nphi+1:I*Nphi)-MATMUL(A_uxux(I,1:Nphi,J*Nphi+1:(J+1)*Nphi)&
                   &,Ux((K-1)*Nphi+1:K*Nphi))
              AUx((I-1)*Nphi+1:I*Nphi)=AUx((I-1)*Nphi+1:I*Nphi)-MATMUL(A_uxuy(I,1:Nphi,J*Nphi+1:(J+1)*Nphi)&
                   &,Uy((K-1)*Nphi+1:K*Nphi))
              AUy((I-1)*Nphi+1:I*Nphi)=AUy((I-1)*Nphi+1:I*Nphi)-MATMUL(A_uyuy(I,1:Nphi,J*Nphi+1:(J+1)*Nphi)&
                   &,Uy((K-1)*Nphi+1:K*Nphi))
              AUy((I-1)*Nphi+1:I*Nphi)=AUy((I-1)*Nphi+1:I*Nphi)-MATMUL(A_uyux(I,1:Nphi,J*Nphi+1:(J+1)*Nphi)&
                   &,Ux((K-1)*Nphi+1:K*Nphi))
           end IF
        END DO
!!!! Les deux lignes suivantes servent juste pour le calcul de l'énergi
!        MPnew((I-1)*Nphi+1:I*Nphi)=MATMUL(M(1:Nphi,1:Nphi),P((I-1)*Nphi+1:I*Nphi))*DFvec(I)
!        MPold((I-1)*Nphi+1:I*Nphi)=MATMUL(M(1:Nphi,1:Nphi),Pold((I-1)*Nphi+1:I*Nphi))*DFvec(I)
!!!!!!!!!!!!!
      END DO  

      Ux(Nsolflu*Nphi+1:(Nsolflu+Nsol)*Nphi)=dt**2*AUx(Nsolflu*Nphi+1:(Nsolflu+Nsol)*Nphi)+&
           &2*Ux(Nsolflu*Nphi+1:(Nsolflu+Nsol)*Nphi)-Uxold(Nsolflu*Nphi+1:(Nsolflu+Nsol)*Nphi)
     Uy(Nsolflu*Nphi+1:(Nsolflu+Nsol)*Nphi)=dt**2*AUy(Nsolflu*Nphi+1:(Nsolflu+Nsol)*Nphi)+&
          &2*Uy(Nsolflu*Nphi+1:(Nsolflu+Nsol)*Nphi)-Uyold(Nsolflu*Nphi+1:(Nsolflu+Nsol)*Nphi)

!I=Nflusol+10
!	Ux((I-1)*Nphi+1:I*Nphi)=Ux((I-1)*Nphi+1:I*Nphi)+dt**2*F(:,1)
if (obstacle==1) then
!	U=U+dt**2*(F(:,1)*sin(omega_incident*N*dt)+0.D0*Fcomp(:,1)*sin(omega_incident*N*dt))
else
if(Isource.gt.Nflu+Nflusol) then
        F(:,1)=ValPhisource(:,1)*Valsource(N+1)
	Ux((Isource-Nflu-Nflusol-1)*Nphi+1:(Isource-Nflu-Nflusol)*Nphi)&
      &=Ux((Isource-Nflu-Nflusol-1)*Nphi+1:(Isource-Nflu-Nflusol)*Nphi)+dt**2*F(:,1)
        F(:,1)=ValPhisource(:,2)*Valsource(N+1)
	Uy((Isource-Nflu-Nflusol-1)*Nphi+1:(Isource-Nflu-Nflusol)*Nphi)&
      &=Uy((Isource-Nflu-Nflusol-1)*Nphi+1:(Isource-Nflu-Nflusol)*Nphi)+dt**2*F(:,1)
end if
endif

DO J=1,Nflusol
I=J+Nflu
Arete=Corres_flusol(J)
K=Neigh(I,Arete)-Nflu-Nflusol
   AP((I-1)*Nphi+1:I*Nphi)=dt**2*AP((I-1)*Nphi+1:I*Nphi)+2*P((I-1)*Nphi+1:I*Nphi)-Pold((I-1)*Nphi+1:I*Nphi)
 AP((I-1)*Nphi+1:I*Nphi)=AP((I-1)*Nphi+1:I*Nphi)+dt/2D0*&
      &matmul(A_pux(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Uxold((K-1)*Nphi+1:K*Nphi))
 AP((I-1)*Nphi+1:I*Nphi)=AP((I-1)*Nphi+1:I*Nphi)+dt/2D0*&
      &matmul(A_puy(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Uyold((K-1)*Nphi+1:K*Nphi))
!!$ AP((I-1)*Nphi+1:I*Nphi)=AP((I-1)*Nphi+1:I*Nphi)-dt**2/2D0*&
!!$      &matmul(A_pux(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Uxold((K-1)*Nphi+1:K*Nphi))
!!$ AP((I-1)*Nphi+1:I*Nphi)=AP((I-1)*Nphi+1:I*Nphi)-dt**2/2D0*&
!!$      &matmul(A_puy(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Uyold((K-1)*Nphi+1:K*Nphi))
END DO
DO J=1,Nsolflu
I=J+Nflu+Nflusol
Arete=Corres_solflu(J)
K=Neigh(I,Arete)
   AUx((J-1)*Nphi+1:J*Nphi)=dt**2*AUx((J-1)*Nphi+1:J*Nphi)+2*Ux((J-1)*Nphi+1:J*Nphi)-Uxold((J-1)*Nphi+1:J*Nphi)
   AUy((J-1)*Nphi+1:J*Nphi)=dt**2*AUy((J-1)*Nphi+1:J*Nphi)+2*Uy((J-1)*Nphi+1:J*Nphi)-Uyold((J-1)*Nphi+1:J*Nphi)
 AUx((J-1)*Nphi+1:J*Nphi)=AUx((J-1)*Nphi+1:J*Nphi)+dt/2D0*&
      &matmul(A_uxp(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),&
      &Pold((K-1)*Nphi+1:K*Nphi))
 AUy((J-1)*Nphi+1:J*Nphi)=AUy((J-1)*Nphi+1:J*Nphi)+dt/2D0*&
      &matmul(A_uyp(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),&
      &Pold((K-1)*Nphi+1:K*Nphi))
!!$ AUx((J-1)*Nphi+1:J*Nphi)=AUx((J-1)*Nphi+1:J*Nphi)-&
!!$      &matmul(A_uxp(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Pold((K-1)*Nphi+1:K*Nphi))
!!$ AUy((J-1)*Nphi+1:J*Nphi)=AUy((J-1)*Nphi+1:J*Nphi)-&
!!$      &matmul(A_uyp(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Pold((K-1)*Nphi+1:K*Nphi))
!!$ AUx((J-1)*Nphi+1:J*Nphi)=AUx((J-1)*Nphi+1:J*Nphi)+2D0*&
!!$      &matmul(A_uxp(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),P((K-1)*Nphi+1:K*Nphi))
!!$ AUy((J-1)*Nphi+1:J*Nphi)=AUy((J-1)*Nphi+1:J*Nphi)+2D0*&
!!$     &matmul(A_uyp(J,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),P((K-1)*Nphi+1:K*Nphi))
END DO
     DO I=1,Nflusol
        J=I+Nflu
        Arete=Corres_flusol(I)
        K=Neigh(J,Arete)-Nflu-Nflusol
        Ftest=AP((J-1)*Nphi+1:J*Nphi)-dt/2D0*matmul(a_pux(I,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),AUx((K-1)*Nphi+1:K*Nphi))
        Ftest=Ftest-dt/2D0*matmul(a_puy(I,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),AUy((K-1)*Nphi+1:K*Nphi))

!!$        Ftest=AP((J-1)*Nphi+1:J*Nphi)-dt**2/2D0*matmul(a_pux(I,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),AUx((K-1)*Nphi+1:K*Nphi))
!!$        Ftest=Ftest-dt**2/2D0*matmul(a_puy(I,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),AUy((K-1)*Nphi+1:K*Nphi))
        !! On divise par (I+dt/2*C)
       CALL DGETRS( 'N',Nphi, 1,Bflusol(I,:,:), Nphi, IPIVflusol(I,:),Ftest,Nphi,INFO )
        P((J-1)*Nphi+1:J*Nphi)=Ftest
        Arete=Corres_solflu(K)
        Ux((K-1)*Nphi+1:K*Nphi)=AUx((K-1)*Nphi+1:K*Nphi)-dt/2D0*matmul(a_uxp(K,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Ftest)
        Uy((K-1)*Nphi+1:K*Nphi)=AUy((K-1)*Nphi+1:K*Nphi)-dt/2D0*matmul(a_uyp(K,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Ftest)
!!$        Ux((K-1)*Nphi+1:K*Nphi)=AUx((K-1)*Nphi+1:K*Nphi)-matmul(a_uxp(K,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Ftest)
!!$        Uy((K-1)*Nphi+1:K*Nphi)=AUy((K-1)*Nphi+1:K*Nphi)-matmul(a_uyp(K,1:Nphi,(Arete-1)*Nphi+1:Arete*Nphi),Ftest)
     ENDDO
  if (iscla==1) then
     DO I=1,NbCla
        J=TriCla(I)
        !! Attention, la j'ai mis la courbure
!        Ftest=P((J-1)*Nphi+1:J*Nphi)+(dt/2.D0-dt**2/4D0/2D0)*matmul(Bcla(I,1:Nphi,1:Nphi),Pold((J-1)*Nphi+1:J*Nphi))
        !! sans la courbure :
        Ftest=P((J-1)*Nphi+1:J*Nphi)+dt/2.D0*matmul(Bcla(I,1:Nphi,1:Nphi),Pold((J-1)*Nphi+1:J*Nphi))/sqrt(rho(J)*mu(J))
        !! On divise par (I+dt/2*C)
       CALL DGETRS( 'N',Nphi, 1,Bcla2(I,:,:), Nphi, IPIV(I,:),Ftest,Nphi,INFO )
        P((J-1)*Nphi+1:J*Nphi)=Ftest
     ENDDO
   end if
flush(22)

  Pold=Pinter
  APold=APinter
  APinter=AP
  MPold=MPinter
  MPinter=MPnew
  Pinter=P

  Uxold=Uxinter
  AUxold=AUxinter
  AUxinter=AUx
  MUxold=MUxinter
  MUxinter=MUxnew
  Uxinter=Ux

  Uyold=Uyinter
  AUyold=AUyinter
  AUyinter=AUy
  MUyold=MUyinter
  MUyinter=MUynew
  Uyinter=Uy

 !Calcul de l'énergie
write(6,*) mod(N*dt*omega_incident,2*pi),Nbinst2,maxval(Ux), minval(Ux)
if(obstacle==1) then
theta_snap=pi/2
IF(abs(N*dt-(2*pi*Nbinst2+theta_snap)/omega_incident).lt.dt) then  
   U_inc=(((2*pi*Nbinst2+theta_snap)/omega_incident-N*dt)*U&
        &+((N+1)*dt-(2*pi*Nbinst2+theta_snap)/omega_incident)*Uold)/dt 
   NBINST2=NBINST2+1 
   Nbinst=1
   call sub_writeunstruct(N)
end IF
else
IF(Mod(N,200).eq.0) then 
  NBINST=NBINST+1
	 call sub_writeunstruct_p
	 call sub_writeunstruct_ux
	 call sub_writeunstruct_uy
END IF
end if
! call sub_writeunstruct(N)
! call sub_write_sismos(N)

!else
!IF(Mod(N,100).eq.0) then 
!  NBINST=NBINST+1
!  call sub_writeunstruct(N)
!END IF

! call sub_writeunstruct(N)
 call sub_write_sismos_acoustic_elasto(N)
!!$IF(Mod(N,100).eq.0) then 
!!$  NBINST2=NBINST2+1
!!$  call sub_writeregular_vero
!!$END IF
 END DO

 close(22)
 close(23)  
 close(24)  
 close(25)  
end program acoustique_scalaire_vectoriel
