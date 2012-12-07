SUBROUTINE sub_allocate
    Use OMP_LIB
    Use m_mat
    Use m_mesh
    Use m_gen
    Use m_inter

    implicit none
    integer :: J

    if (helmholtz.eq.0) then
        allocate(A(Nphi,4*Nphi,Ntri))
        allocate(Us(Nphi*Ntri,2))
        allocate(AU(Nphi*Ntri))
        allocate(Neigh(NTri,3))

        U    => Us(:,1)
        Uold => Us(:,2)
!        AU   => U_all(:,3)
        
!$OMP PARALLEL DO
        DO J=1,Ntri
            A(:,:,J)                  = 0.0d0
            Us((J-1)*Nphi+1:J*Nphi,:) = 0.0d0
            AU((J-1)*Nphi+1:J*Nphi)   = 0.0d0
            Neigh(J,:)                = 0
        END DO
        
        allocate(Corres_arete(Ntri,1+3*(1+order)))
        allocate(Corres_arete_u(Ntri,1+3*Nphi))
        allocate(DFVec(Ntri))


    else
    
        allocate(A(Nphi,4*Nphi,Ntri))
        allocate(Us(Nphi*Ntri,1))
        allocate(Neigh(NTri,3))

        U    => Us(:,1)
        
!$OMP PARALLEL DO
        DO J=1,Ntri
            A(:,:,J)               = 0.0d0
            U((J-1)*Nphi+1:J*Nphi) = 0.0d0
            Neigh(J,:)             = 0
        END DO
        
        allocate(Corres_arete(Ntri,1+3*(1+order)))
        allocate(Corres_arete_u(Ntri,1+3*Nphi))
        allocate(DFVec(Ntri))
        
    end if

END SUBROUTINE sub_allocate
