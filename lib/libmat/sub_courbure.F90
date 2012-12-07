subroutine sub_courbure
Use m_condbord  
 Use m_mat
  Use m_gen
  Use m_mesh  
 implicit none 
 real*8 :: courbure1,courbure2,aire_arete
 real*8 :: x1,y1,x,y,x2,y2,x3,y3
  real*8 :: x_centre,y_centre,d1,d2,d3,dtmp
  integer :: I,II,III,K,J,PhiEdge(3,2)
  
  
  
   PhiEdge(1,1)=2
     PhiEdge(1,2)=3
     PhiEdge(2,1)=3
     PhiEdge(2,2)=1
     PhiEdge(3,1)=1
     PhiEdge(3,2)=2
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!           CONSTRUCTION DES ARETES DU BORD
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!Nombre d'arêtes au bord égal au nombre de triangles au bord
  Narete_du_bord=0
  DO I=1,NTri
    IF (Condbord==1) then
       DO J=1,3
      
	   If  ( (CondPoint(Tri(I,PhiEdge(J,1))).eq.3) .and.  (CondPoint(Tri(I,PhiEdge(J,2))).eq.3) ) then
		Narete_du_bord=Narete_du_bord+1
	
	   end if
       END DO
    END IF
  ENDDO

!!Nombre de triangle qui ne sont pas sur le bord -3
Ntri_non_bord=Ntri-Narete_du_bord
write(6,*) Ntri_non_bord
write(6,*) Narete_du_bord
write(6,*) Ntri



allocate(Arete_du_bord(Narete_du_bord,2))
allocate(tri_bord(Narete_du_bord))
allocate(tri_non_bord(Ntri_non_bord))
allocate(test_bord(Ntri))

test_bord=0.D0

!! On stoke le numéro global des triangles sur le bord et on construit
!! le vecteur Arete_du_bord qui contient les 2 noeuds et on def test_bord

II=1
III=1
 DO I=1,NTri
    IF (Condbord==1) then
       DO J=1,3
	   If  ( (CondPoint(Tri(I,PhiEdge(J,1))).eq.3) .and.  (CondPoint(Tri(I,PhiEdge(J,2))).eq.3) ) then
		tri_bord(II)=I
		test_bord(I)=1		
		Arete_du_bord(II,1)=Tri(I,PhiEdge(J,1))
		Arete_du_bord(II,2)=Tri(I,PhiEdge(J,2))
		II=II+1
		
                exit
		
                                
           end if   
       END DO
    else 
    		    tri_non_bord(III)=I
                    III=III+1
    end if
    
  ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!              CONSTRUCTION DES VOISINS DU BORD
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(Neigh_bord(Narete_du_bord,2))

!!Pour chaque arete on donne le voisin de gauche et de droite
!!Si l'arete est au coin le numéro du voisin est -50
 If (Narete_du_bord.ne.0) then
 	DO I=1,Narete_du_bord-1
	   DO K=1,2
                If  (CondPoint(Arete_du_bord(I,K)).eq.3) then
			x1=coor(Arete_du_bord(I,K),1)			
			y1=coor(Arete_du_bord(I,K),2)
                        DO J=I+1,Narete_du_bord
                                x2=coor(Arete_du_bord(J,1),1)
				y2=coor(Arete_du_bord(J,1),2)
				x3=coor(Arete_du_bord(J,2),1)
                                y3=coor(Arete_du_bord(J,2),2)
				if ( (x1.eq.x2).and. (y1.eq.y2) ) then
                                        Neigh_bord(I,K)=J
					Neigh_bord(J,1)=I
                                        exit
					elseif  ( (x1.eq.x3).and. (y1.eq.y3) ) then
							Neigh_bord(I,K)=J
							Neigh_bord(J,2)=I
                                                        exit
				end if
			end do
		else 
			Neigh_bord(I,K)=-50
			write(6,*)'il y a un coin!!'
		end if
	   end do
	end do

	
end if

allocate(courbure(Narete_du_bord))
courbure=0.D0
If (Narete_du_bord.ne.0) then
!! On parcourt les aretes du bord
	DO I=1,Narete_du_bord
		!!coordonnees du premier point
		x1=coor(Arete_du_bord(I,1),1)
		y1=coor(Arete_du_bord(I,1),2)
		!! coordonnées du deuxieme point
		x2=coor(Arete_du_bord(I,2),1)
                y2=coor(Arete_du_bord(I,2),2)
		!! numero de l'arete voisine connectée par (I,1)
		J=Neigh_bord(I,1)
		!! numero de l'arete voisine connectée par (I,2)
		K=Neigh_bord(I,2)
		!!COURBURE AU PREMIER POINT
		!!coordonnées de (J,1)
		x3=coor(Arete_du_bord(J,1),1)
		y3=coor(Arete_du_bord(J,1),2)
		!! test si ce point correspond à (I,1)
		if ( (x1.eq.x3) .and. (y1.eq.y3) ) then
			!! le voisin est (J,2) donc on stocke ses coordonnées
			x3=coor(Arete_du_bord(J,2),1)
			y3=coor(Arete_du_bord(J,2),2)
			!! a priori il faut faire un test pour savoir s'ils sont alignés (mais bon c'est pas codé)
			if ( abs((y3-(y2-y1)*x3/(x2-x1)-y1+(y2-y1)*x1/(x2-x1))).le. 1e-6 ) then
			!	write(6,*) 'kk'
				!! les points sont alignés
				courbure1=0
			elseif ( (x3.eq.x1) .and. (x3.eq.x2) ) then
					!! les point sont alignés
					courbure1=0
			else
			!! on calcule les trois longueurs des cotes du triangle
                        d1=sqrt((x2-x1)**2+(y2-y1)**2)
                        d2=sqrt((x3-x1)**2+(y3-y1)**2)
                        d3=sqrt((x2-x3)**2+(y2-y3)**2)
                        !! on les range pour que d3>d2>d1
                        if (d1.ge.d2) then
                           dtmp=d1
                           d1=d2
                           d2=dtmp
                        endif
                        if (d1.ge.d3) then
                           dtmp=d3
                           d3=d2
                           d2=d1
                           d1=dtmp
                           elseif (d2.ge.d3) then
                              dtmp=d2
                              d2=d3
                              d3=dtmp
                        endif
                        aire_arete=sqrt((d3+d2+d1)*(d3+d2-d1)*(d1+d3-d2)*(d1-d3+d2))      
			!! on calcule la courbure aui est égale à 1/(rayon de ce cercle)
			courbure1=aire_arete/(d1*d2*d3)
			endif
		else !! le voisin est (J,1)
			!! a priori il faut faire un test pour savoir s'ils sont alignés (mais bon c'est pas codé)
			if ( abs((y3-(y2-y1)*x3/(x2-x1)-y1+(y2-y1)*x1/(x2-x1))).le. 1e-6 ) then
			!	write(6,*) 'kk'
				!! les points sont alignés
				courbure1=0
			elseif ( (x3.eq.x1) .and. (x3.eq.x2) ) then
					!! les point sont alignés
					courbure1=0
			else
			!! on calcule les trois longueurs des cotes du triangle
                        d1=sqrt((x2-x1)**2+(y2-y1)**2)
                        d2=sqrt((x3-x1)**2+(y3-y1)**2)
                        d3=sqrt((x2-x3)**2+(y2-y3)**2)
                        !! on les range pour que d3>d2>d1
                        if (d1.ge.d2) then
                           dtmp=d1
                           d1=d2
                           d2=dtmp
                        endif
                        if (d1.ge.d3) then
                           dtmp=d3
                           d3=d2
                           d2=d1
                           d1=dtmp
                           elseif (d2.ge.d3) then
                              dtmp=d2
                              d2=d3
                              d3=dtmp
                        endif
                        aire_arete=sqrt((d3+d2+d1)*(d3+d2-d1)*(d1+d3-d2)*(d1-d3+d2))
			!! on calcule la courbure aui est égale à 1/(rayon de ce cercle)
			courbure1=aire_arete/(d1*d2*d3)
			endif
		end if

		!!COURBURE AU DEUXIEME POINT
		!!coordonnées de (K,1)
		x3=coor(Arete_du_bord(K,1),1)
		y3=coor(Arete_du_bord(K,1),2)
		!! test si ce point correspond à (I,1)
		if ( (x2.eq.x3) .and. (y2.eq.y3) ) then
			!! le voisin est (K,2) donc on stocke ses coordonnées
			x3=coor(Arete_du_bord(K,2),1)
			y3=coor(Arete_du_bord(K,2),2)
			!! a priori il faut faire un test pour savoir s'ils sont alignés (mais bon c'est pas codé)
			if ( abs((y3-(y2-y1)*x3/(x2-x1)-y1+(y2-y1)*x1/(x2-x1))).le. 1e-6 ) then
				!write(6,*) 'kk'
				!! les points sont alignés
write(6,*) 'heho'
				courbure2=0
			elseif ( (x3.eq.x1) .and. (x3.eq.x2) ) then
					!! les point sont alignés
					courbure2=0
			else
			!! on calcule les trois longueurs des cotes du triangle
                        d1=sqrt((x2-x1)**2+(y2-y1)**2)
                        d2=sqrt((x3-x1)**2+(y3-y1)**2)
                        d3=sqrt((x2-x3)**2+(y2-y3)**2)
                        !! on les range pour que d3>d2>d1
                        if (d1.ge.d2) then
                           dtmp=d1
                           d1=d2
                           d2=dtmp
                        endif
                        if (d1.ge.d3) then
                           dtmp=d3
                           d3=d2
                           d2=d1
                           d1=dtmp
                           elseif (d2.ge.d3) then
                              dtmp=d2
                              d2=d3
                              d3=dtmp
                        endif
                        aire_arete=sqrt((d3+d2+d1)*(d3+d2-d1)*(d1+d3-d2)*(d1-d3+d2))
			!! on calcule la courbure aui est égale à 1/(rayon de ce cercle)
			courbure2=aire_arete/(d1*d2*d3)
			endif
		else !! le voisin est (K,1)
			!! a priori il faut faire un test pour savoir s'ils sont alignés (mais bon c'est pas codé)
			if ( abs((y3-(y2-y1)*x3/(x2-x1)-y1+(y2-y1)*x1/(x2-x1))).le. 1e-6 ) then
			!	write(6,*) 'kk'
				!! les points sont alignés

				courbure2=0
			elseif ( (x3.eq.x1) .and. (x3.eq.x2) ) then
					!! les point sont alignés
					courbure2=0
			else
			!! on calcule les trois longueurs des cotes du triangle
                        d1=sqrt((x2-x1)**2+(y2-y1)**2)
                        d2=sqrt((x3-x1)**2+(y3-y1)**2)
                        d3=sqrt((x2-x3)**2+(y2-y3)**2)
                        !! on les range pour que d3>d2>d1
                        if (d1.ge.d2) then
                           dtmp=d1
                           d1=d2
                           d2=dtmp
                        endif
                        if (d1.ge.d3) then
                           dtmp=d3
                           d3=d2
                           d2=d1
                           d1=dtmp
                           elseif (d2.ge.d3) then
                              dtmp=d2
                              d2=d3
                              d3=dtmp
                        endif
                        aire_arete=sqrt((d3+d2+d1)*(d3+d2-d1)*(d1+d3-d2)*(d1-d3+d2))      
			!! on calcule la courbure aui est égale à 1/(rayon de ce cercle)
			courbure2=aire_arete/(d1*d2*d3)
			endif
		end if
	courbure(I)=(courbure1+courbure2)/2.D0
!if (courbure(I).ne.0) then
!write(6,*) I, courbure(I),d1,d2,d3,aire_arete,courbure1,courbure2
!end if
	end do
end if
end
