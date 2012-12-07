SUBROUTINE creation_grille_diff_finies
 
  Use m_condbord    
  Use m_mat
  Use m_gen
  Use m_mesh
  Use m_output
  Use m_source
  Use m_inter
 ! Use m_rampant


  IMPLICIT NONE
  integer :: I,tmp
  real*8 :: x1,x2,y1,y2



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	 !!                                                                          !!   
	 !!                     LISTE DES SOMMETS DU BORD                            !!
	 !!                                                                          !!
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! liste_sommets_bord contient le numero global de tous les sommets des triangles du bord non répétés
  !! Pour cela on ne stocke qu'un noeud sur 2 mais vu qu'une arete peut etre stockée dans les 2 sens
  !! le noeud à stocker n'est pas forcément tout le temps le 1er de chaque arete ou le 2nd.
  allocate(liste_sommets_bord(Narete_du_bord))

  !! le 1er noeud correspond à la 1ere extrémité de la 1ere arete du bord
  liste_sommets_bord(1)=Arete_du_bord(1,1)
!  liste_sommets_bord(1)=Arete_du_bord(1,2)

  !! tmp sera le noeud à stocker mais il nous permet de vérifier s'il correspond à la 1ere ou la 2eme 
  !! extrémité de l'arête suivante
  tmp=Arete_du_bord(1,2)
!  tmp=Arete_du_bord(1,1)

  !! Maintenant que c'est initialisé, on parcourt tous les triangles du bord
  DO I=2,Narete_du_bord
	!! On regarde à quelle extrémité de l'arête I correspond tmp
	if (tmp==Arete_du_bord(I,1)) then
		liste_sommets_bord(I)=tmp
		tmp=Arete_du_bord(I,2)
	elseif (tmp==Arete_du_bord(I,2)) then
		liste_sommets_bord(I)=tmp
		tmp=Arete_du_bord(I,1)
	else
		stop
	endif
  END DO

write(6,*) 'Ok?'


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	 !!                                                                          !!   
	 !!                      LISTE DES NOEUDS POUR DF                            !!
	 !!                                                                          !!
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Le nombre de noeud sur le bord dépend de l'ordre des EF : si l'ordre choisit est 3, il nous
  !! faut considérer 2 noeuds entre chaque sommet
  !! liste_noeuds_DF contient les coordonnées de chaque noeuds et non pas le numéro des noeuds car
  !! les noeuds à rajouter n'ont pas de numéro global !!!
  !! En plus liste_noeuds_DF contient un marqueur qui précise si ce noeuds est issu d'un sommet ou non

  !! Vu que la liste des noeuds DF dépend de l'ordre on va considérer les 3 cas
  SELECT CASE (Order)

	!! On choisit de faire du P1
	CASE (1)
		!! Dans ce cas, il n'y a pas de noeuds intermédiaires liste_noeud_DF correspond 
		!! à liste_sommets_bord
		allocate(liste_noeuds_DF(Narete_du_bord,3))
		!! On parcourt chaque sommets du bord
		DO I=1,Narete_du_bord
			!! chaque noeud est un sommet donc les 2 1ere colonnes sont les coor des sommets
			liste_noeuds_DF(I,1)=Coor(liste_sommets_bord(I),1)
			liste_noeuds_DF(I,2)=Coor(liste_sommets_bord(I),2)
			!! vu que tous les noeuds sont des sommets, le marqueur est tout le temps égal à 1
			liste_noeuds_DF(I,3)=1
		END DO  
		 nombre_noeuds_DF=Narete_du_bord

	!! On choisit de faire du P2
	CASE (2)
		!! Dans ce cas il faut rajouter un noeud entre chaque sommet
		allocate(liste_noeuds_DF(2*Narete_du_bord,3))
		!! On parcourt chaque arete du bord
		DO I=1,Narete_du_bord
			!! le (2I-1)ème noeud correspond au sommet I
			liste_noeuds_DF(2*I-1,1)=Coor(liste_sommets_bord(I),1)
			liste_noeuds_DF(2*I-1,2)=Coor(liste_sommets_bord(I),2)
			liste_noeuds_DF(2*I-1,3)=1
			!! le (2I)ème correspond au milieu de l'arête initiale I
			!! On récupère les coordonnées des 2 extrémités de l'arête du bord I
			x1=Coor(Arete_du_bord(I,1),1)
			y1=Coor(Arete_du_bord(I,1),2) 
			x2=Coor(Arete_du_bord(I,2),1)
			y2=Coor(Arete_du_bord(I,2),2)
			!! On calcule le milieu
			liste_noeuds_DF(2*I,1)=x1+(x2-x1)/2
			liste_noeuds_DF(2*I,2)=y1+(y2-y1)/2
			liste_noeuds_DF(2*I,3)=0
		END DO
		 nombre_noeuds_DF=2*Narete_du_bord
	!! On choisit de faire du P3
	CASE (3)
		!! Dans ce cas il faut rajouter un noeud entre chaque sommet
		allocate(liste_noeuds_DF(3*Narete_du_bord,3))
		!! On parcourt chaque arete du bord
		DO I=1,Narete_du_bord
			!! le (3I-2)ème noeud correspond au sommet I
			liste_noeuds_DF(3*I-2,1)=Coor(liste_sommets_bord(I),1)
			liste_noeuds_DF(3*I-2,2)=Coor(liste_sommets_bord(I),2)
			liste_noeuds_DF(3*I-2,3)=1
			!! le (3I-1)ème et le (3I)ème correspondent au tiers et 2tiers de l'arête I
			!! On récupère les coordonnées des 2 extrémités de l'arête du bord I
			x1=Coor(Arete_du_bord(I,1),1)
			y1=Coor(Arete_du_bord(I,1),2) 
			x2=Coor(Arete_du_bord(I,2),1)
			y2=Coor(Arete_du_bord(I,2),2)
			!! On calcule les coordonnées et le marqueur = 0
			liste_noeuds_DF(3*I-1,1)=x1+(x2-x1)/3
			liste_noeuds_DF(3*I-1,2)=y1+(y2-y1)/3
			liste_noeuds_DF(3*I-1,3)=0

			liste_noeuds_DF(3*I,1)=x1+2*(x2-x1)/3
			liste_noeuds_DF(3*I,2)=y1+2*(y2-y1)/3
			liste_noeuds_DF(3*I,3)=0
		END DO
		 nombre_noeuds_DF=3*Narete_du_bord
	CASE DEFAULT
		write(6,*) 'Not implemented'
		stop	




  END SELECT


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	 !!                                                                          !!   
	 !!                          PAS D'ESPACE POUR DF                            !!
	 !!                                                                          !!
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  

  !! pas_espace_DF est le vecteur qui contient tous les pas d'espace :
  !!        pas_espace_DF(I)= distance entre liste_noeuds_DF(I) et liste_noeuds_DF(I+1)

  !! On récupère le nombre de noeuds_DF : on utilise size qui compte le nombre d'éléments de la matrice
  !! et on divise par le nombre de colonnes à savoir 3
 ! nombre_noeuds_DF=size(liste_noeuds_DF)/3
 
!write(6,*) nombre_noeuds_DF
  !! Initialisation du vecteur pas_espace_DF
  allocate(pas_espace_DF(nombre_noeuds_DF))

  DO I=1,nombre_noeuds_DF-1
	x1=liste_noeuds_DF(I,1)
	y1=liste_noeuds_DF(I,2)
	x2=liste_noeuds_DF(I+1,1)
	y2=liste_noeuds_DF(I+1,2)
	pas_espace_DF(I)=sqrt((x2-x1)**2+(y2-y1)**2)
!	write(6,*) I,x1,y1,x2,y2,pas_espace_DF(I)
  END DO
!stop
  x1=liste_noeuds_DF(nombre_noeuds_DF,1)
  y1=liste_noeuds_DF(nombre_noeuds_DF,2)
  x2=liste_noeuds_DF(1,1)
  y2=liste_noeuds_DF(1,2)
  pas_espace_DF(nombre_noeuds_DF)=sqrt((x2-x1)**2+(y2-y1)**2)
!write(6,*) pas_espace_DF
!stop

pas_espace_DF_max=0.D0
DO I=1,nombre_noeuds_DF
	if (pas_espace_DF(I)>pas_espace_DF_max) then  
		pas_espace_DF_max=pas_espace_DF(I)
	endif
ENDDO
END SUBROUTINE creation_grille_diff_finies
