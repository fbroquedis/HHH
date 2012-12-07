SUBROUTINE sub_defmatD_4
  Use m_mat
  Use m_gen
  Use m_mesh
  Use m_condbord
  implicit none
  INTEGER :: I,INeigh,Node(2),J,K,L,L1,noeud,JJ,testcla,Q,II,Node1(3),JNeigh,L2,toto,L1E,L2E,JJ2,TOTOE
  REAL*8 :: Jfinv,N(2,2),V(3,2),cmax
  REAL*8,allocatable::DerivPhiIDerivPhiJ(:,:)
  REAL*8,allocatable::DerivPhiI(:,:)
  REAL*8,allocatable :: DerivPhi1D(:,:),DerivPhi0D(:,:,:)
  REAL*8 :: norm_point(2),h1,res,hmin,h2
  INTEGER,allocatable :: PhiEdge(:,:)
  
  allocate(PhiEdge(3,1+Order))
  allocate(DerivPhi1D(1+Order,NGL1D))
  allocate(DerivPhi0D(2,1+Order,1))

  SELECT CASE(Order)
  CASE(1)

     CALL DerivPhi1DOrder1(PtGL1D,NGL1D,DerivPhi1D)
     !On node 1
     CALL DerivPhi1DOrder1(0.D0,1,DerivPhi0D(1,:,:))
     !On node 2
     CALL DerivPhi1DOrder1(1.D0,1,DerivPhi0D(2,:,:))
     !Functions on Edge 23
     PhiEdge(1,1)=2
     PhiEdge(1,2)=3
     !Functions on Edge 31
     PhiEdge(2,1)=3
     PhiEdge(2,2)=1
     !Functions on Edge 12
     PhiEdge(3,1)=1
     PhiEdge(3,2)=2
 
 
  CASE(2)
     CALL DerivPhi1DOrder2(PtGL1D,NGL1D,DerivPhi1D)     
!On node 1 
     CALL DerivPhi1DOrder2(0.D0,1,DerivPhi0D(1,:,:))
!On node 2 
     CALL DerivPhi1DOrder2(1.D0,1,DerivPhi0D(2,:,:))
!Functions on Edge 23
PhiEdge(1,1)=2
PhiEdge(1,2)=5
PhiEdge(1,3)=3
!Functions on Edge 31
PhiEdge(2,1)=3
PhiEdge(2,2)=6
PhiEdge(2,3)=1
!Functions on Edge 12
PhiEdge(3,1)=1
PhiEdge(3,2)=4
PhiEdge(3,3)=2


  CASE(3)
     CALL DerivPhi1DOrder3(PtGL1D,NGL1D,DerivPhi1D)     
!On node 1 
     CALL DerivPhi1DOrder3(0.D0,1,DerivPhi0D(1,:,:))
!On node 2
     CALL DerivPhi1DOrder3(1.D0,1,DerivPhi0D(2,:,:))
!Functions on Edge 23
PhiEdge(1,1)=2
PhiEdge(1,2)=6
PhiEdge(1,3)=7
PhiEdge(1,4)=3
!Functions on Edge 31
PhiEdge(2,1)=3
PhiEdge(2,2)=8
PhiEdge(2,3)=9
PhiEdge(2,4)=1
!Functions on Edge 12
PhiEdge(3,1)=1
PhiEdge(3,2)=4
PhiEdge(3,3)=5
PhiEdge(3,4)=2

  END SELECT


allocate(DerivPhiI(2,1+Order))
DO L=1,2
   DO I=1,1+Order
            DerivPhiI(L,I)= sum(DerivPhi0D(L,I,:)*1) 
   END DO
END DO
allocate(DFVEC_arete(Narete_du_bord))
allocate(DerivPhiIDerivPhiJ(1+Order,1+Order))
DO I=1,1+Order
   DO J=1,1+Order
      DerivPhiIDerivPhiJ(I,J)=sum(DerivPhi1D(I,:)*DerivPhi1D(J,:)*wGL1D)
   END DO
END DO
NBcla=0
DO I=1,Ntri
   DO J=1,3
      IF (Neigh(I,J)==-3) THEN
         NbCLA=NbCLA+1
exit
      END IF
   ENDDO
ENDDO
write(6,*) Nbcla,Narete_du_bord
allocate(Dclabis(Narete_du_bord,1+order,3*(1+order)))
allocate(Ebis(nbcla,1+order,3*Nphi))
Dclabis=0.D0
Ebis=0.D0
Nbcla=0
DO I=1,NTri
   Node1 = Tri(I,:)
   !!Computation of vectors V12,V23 and V31
   V(1,1)=Coor(Node1(3),1)-Coor(Node1(2),1)
   V(1,2)=Coor(Node1(3),2)-Coor(Node1(2),2)
   V(2,1)=Coor(Node1(1),1)-Coor(Node1(3),1)
   V(2,2)=Coor(Node1(1),2)-Coor(Node1(3),2)
   V(3,1)=Coor(Node1(2),1)-Coor(Node1(1),1)
   V(3,2)=Coor(Node1(2),2)-Coor(Node1(1),2)
   TestCla=0
  DO JJ=1,3
!!! Searching Neighbor
      INeigh=Neigh(I,JJ)
!!! If there is a Neighbor
      IF(Ineigh.eq.-3) then
         If (TestCla==0) then
            TestCla=1
            NbCla=NbCla+1
         END If
           DO J=1,1+Order
               DO K=1,1+Order
                  res=DerivPhiIDerivPhiJ(J,K)
                  Dclabis(Nbcla,J,K)=Dclabis(Nbcla,J,K)+res/sqrt(sum(V(JJ,:)**2))/rho(I)
                  L=PhiEdge(JJ,K)
                  Ebis(Nbcla,J,L)=Ebis(Nbcla,J,L)+res/sqrt(sum(V(JJ,:)**2))/rho(I)
		end do
	  end do
        Node = Arete_du_bord(Nbcla,:)
	! construction des correspondances
        Corres_arete(nbcla,1)=1+Order
        Corres_arete_u(nbcla,1)=Nphi
        DO J=1,1+Order
           Corres_arete(nbcla,J+1)=(1+order)*(nbcla-1)+J
        END DO
        DO J=1,Nphi
           Corres_arete_u(nbcla,J+1)=Nphi*(I-1)+J
        END DO

	!!Computation of nodes V1,V2
     N(1,1)=Coor(Node(1),1)
     N(1,2)=Coor(Node(1),2)
     N(2,1)=Coor(Node(2),1)
     N(2,2)=Coor(Node(2),2)
     h1=sqrt((N(2,1)-N(1,1))**2+(N(2,2)-N(1,2))**2)
	!!Computation of normal vectors 
   IF (Coor(Node(1),1).gt.Coor(Node(2),1)) then
	 norm_point(1)=-1
	 norm_point(2)=1
	else 
	     norm_point(1)=1
	     norm_point(2)=-1
   end if

!! Computation of JFinv and DF
    DFVEC_arete(nbcla)=h1

      Jfinv=1/h1

     L1=1+Order
     L1E=Nphi
!! Searching Neighbor
     JNeigh=Neigh_bord(nbcla,1)
     IF(Jneigh.gt.0) then
        Node = Arete_du_bord(Jneigh,:)
        N(1,1)=Coor(Node(1),1)
        N(1,2)=Coor(Node(1),2)
        N(2,1)=Coor(Node(2),1)
        N(2,2)=Coor(Node(2),2)
        h2=sqrt((N(2,1)-N(1,1))**2+(N(2,2)-N(1,2))**2)
        hmin=min(h1,h2)
        cmax=max(1D0/rho(tricla(nbcla)),1D0/(tricla(Jneigh)))
      !! Recherche du numéro de l'arete voisine
        if (tri(tricla(Jneigh),1).eq.Node(1)) then
            if(tri(tricla(Jneigh),2).eq.Node(2)) then 
                 JJ2=3
            elseif (tri(tricla(Jneigh),3).eq.Node(2)) then
                     JJ2=2
            else
               stop
            end if
        elseif (tri(tricla(Jneigh),2).eq.Node(1)) then
            if(tri(tricla(Jneigh),1).eq.Node(2)) then 
                 JJ2=3
            elseif (tri(tricla(Jneigh),3).eq.Node(2)) then
                     JJ2=1
            else
               stop
            end if
        elseif (tri(tricla(Jneigh),3).eq.Node(1)) then
            if(tri(tricla(Jneigh),1).eq.Node(2)) then 
                 JJ2=2
            elseif (tri(tricla(Jneigh),2).eq.Node(2)) then
                     JJ2=1
            else
               stop
            end if
         else
            stop
         end if


        if (Neigh_bord(Jneigh,1).eq.nbcla) then
            toto=1
            L2=1+order
            L2E=Nphi
        else
            toto=1+order
            if (Neigh_bord(Jneigh,1).ne.-50) then
                 L2=2+2*order
                 L2E=2*Nphi
            else
                 L2=1+order
                 L2E=Nphi
            end if
        end if
	 if (Tri(TriCla(Jneigh),1).eq.nbcla) then
            totoE=1
        elseif (Tri(TriCla(Jneigh),2).eq.nbcla) then
            totoE=1+order
        else
            totoE=2+order
            
        end if
        Corres_arete(nbcla,1)=Corres_arete(nbcla,1)+1+Order
	Corres_arete_u(nbcla,1)=Corres_arete_u(nbcla,1)+Nphi
        DO J=1,1+Order
           Corres_arete(nbcla,J+1+L1)=(1+Order)*(JNeigh-1)+J
        END DO
        DO J=1,Nphi
           Corres_arete_u(nbcla,J+1+L1E)=Nphi*(Tricla(Jneigh)-1)+J
        End do
!! Computation of -{grad u}. [v]-{grad v}.[u]

	Do J=1,1+order
        res=(Jfinv*DerivPhiI(1,J)*norm_point(1))
        Dclabis(nbcla,J,1)=Dclabis(nbcla,J,1)+res/(2D0*rho(I))
	L=phiedge(JJ,1)
        Ebis(nbcla,J,L)=Ebis(nbcla,J,L)+res/(2D0*rho(I))        
	Dclabis(nbcla,1,J)=Dclabis(nbcla,1,J)+res/(2D0*rho(I))
	L=phiedge(JJ,J)
        Ebis(nbcla,1,L)=Ebis(nbcla,1,L)+res/(2D0*rho(I))

      

        res=-(Jfinv*DerivPhiI(1,J)*norm_point(1))
        Dclabis(nbcla,J,L1+toto)=Dclabis(nbcla,J,L1+toto)+res/(2D0*rho(I))
	L=phiedge(JJ2,toto)
        Ebis(nbcla,J,L1E+L)=Ebis(nbcla,J,L1E+L)+res/(2D0*rho(I))
        Dclabis(JNeigh,toto,L2+J)=Dclabis(JNeigh,toto,L2+J)+res/(2D0*rho(I))
	L=phiedge(JJ,J)
        Ebis(Jneigh,toto,L2E+L)=Ebis(Jneigh,toto,L2E+L)+res/(2D0*rho(I))
        end do

 !!alpha/hmin [u][v]
            res=alpha2*cmax**2/hmin            
            Dclabis(nbcla,1,1)=Dclabis(nbcla,1,1)+res
            L=phiedge(JJ,1)
            Ebis(nbcla,1,L)=Ebis(nbcla,1,L)+res 
            Dclabis(nbcla,1,L1+toto)=Dclabis(nbcla,1,L1+toto)-res
            L=phiedge(JJ2,toto)
            Ebis(nbcla,1,L1E+L)=Ebis(nbcla,1,L1E+L)-res 
    
     L1=L1+1+Order
     L1E=L1E+Nphi

end if
!!!!!!!!!!!!!!!!!!
!! 2eme noeud
!!!!!!!!!!!!!!!!!!

     JNeigh=Neigh_bord(nbcla,2)

!!! If there is a Neighbor
     IF(Jneigh.gt.0) then
        Node = Arete_du_bord(Jneigh,:)
        N(1,1)=Coor(Node(1),1)
        N(1,2)=Coor(Node(1),2)
        N(2,1)=Coor(Node(2),1)
        N(2,2)=Coor(Node(2),2)
        h2=sqrt((N(2,1)-N(1,1))**2+(N(2,2)-N(1,2))**2)
        hmin=min(h1,h2)
        cmax=max(1D0/rho(tricla(nbcla)),1D0/rho(tricla(Jneigh)))
      !! Recherche du numéro de l'arete voisine
        if (tri(tricla(Jneigh),1).eq.Node(1)) then
            if(tri(tricla(Jneigh),2).eq.Node(2)) then 
                 JJ2=3
            elseif (tri(tricla(Jneigh),3).eq.Node(2)) then
                     JJ2=2
            else
               stop
            end if
        elseif (tri(tricla(Jneigh),2).eq.Node(1)) then
            if(tri(tricla(Jneigh),1).eq.Node(2)) then 
                 JJ2=3
            elseif (tri(tricla(Jneigh),3).eq.Node(2)) then
                     JJ2=1
            else
               stop
            end if
        elseif (tri(tricla(Jneigh),3).eq.Node(1)) then
            if(tri(tricla(Jneigh),1).eq.Node(2)) then 
                 JJ2=2
            elseif (tri(tricla(Jneigh),2).eq.Node(2)) then
                     JJ2=1
            else
               stop
            end if
         else
            stop
         end if

        if (Neigh_bord(Jneigh,1).eq.nbcla) then
            toto=1
            L2=1+order
            L2E=Nphi
        else
            toto=1+order
            if (Neigh_bord(Jneigh,1).ne.-50) then
                 L2=2+2*order
                 L2E=2*Nphi
            else
                 L2=1+order
                 L2E=Nphi
            end if
        end if

        if (Tri(TriCla(Jneigh),1).eq.nbcla) then
            totoE=1
        elseif (Tri(TriCla(Jneigh),3).eq.nbcla) then
            totoE=1+order
        else 
           totoE=2+2*order 
       end if
        Corres_arete(nbcla,1)=Corres_arete(nbcla,1)+1+Order
        Corres_arete_u(nbcla,1)=Corres_arete_u(nbcla,1)+Nphi
        DO J=1,1+Order
           Corres_arete(nbcla,J+1+L1)=(1+Order)*(JNeigh-1)+J
        END DO
        Do J=1,Nphi
           Corres_arete_u(nbcla,J+1+L1E)=Nphi*(tricla(Jneigh)-1)+J
        END DO

	Do J=1,1+order
	res=(Jfinv*DerivPhiI(2,J)*norm_point(2))
        Dclabis(nbcla,J,1+order)=Dclabis(nbcla,J,1+order)+res/(2D0*rho(I))
        L=phiedge(JJ,1+order)
        Ebis(nbcla,J,L)=Ebis(nbcla,J,L)+res/(2D0*rho(I))
        Dclabis(nbcla,1+Order,J)=Dclabis(nbcla,1+Order,J)+res/(2D0*rho(I))
        L=phiedge(JJ,J)
        Ebis(nbcla,1+order,L)=Ebis(nbcla,1+order,L)+res/(2D0*rho(I)) 

        res=-(Jfinv*DerivPhiI(2,J)*norm_point(2))
        Dclabis(nbcla,J,L1+toto)=Dclabis(nbcla,J,L1+toto)+res/(2D0*rho(I))
        L=phiedge(JJ2,toto)
        Ebis(nbcla,J,L1E+L)=Ebis(nbcla,J,L1E+L)+res/(2D0*rho(I))
        Dclabis(JNeigh,toto,L2+J)=Dclabis(JNeigh,toto,L2+J)+res/(2D0*rho(I))
        L=phiedge(JJ,J)
        Ebis(Jneigh,toto,L2E+L)=Ebis(Jneigh,toto,L2E+L)+res/(2D0*rho(I)) 
	end do

 !!alpha/hmin [u][v]
      
        res=alpha2*cmax**2/hmin
        Dclabis(nbcla,1+order,1+order)=Dclabis(nbcla,1+order,1+order)+res
        L=phiedge(JJ,1+order)
        Ebis(nbcla,1+order,L)=Ebis(nbcla,1+order,L)+res
        Dclabis(nbcla,1+order,L1+toto)=Dclabis(nbcla,1+order,L1+toto)-res
        L=phiedge(JJ2,toto)
        Ebis(nbcla,1+order,L1E+L)=Ebis(nbcla,1+order,L1E+L)-res 
     
end if

	end if

  END DO
END DO

END SUBROUTINE sub_defmatD_4
