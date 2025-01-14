module particle_displacement_module
    use shape_module
    use read_data_module
    implicit none

    contains

    subroutine particleDisplacement(U, PARTICLES, KSIETA, ELEMS, gaprox, npe, npe_p, &
                                     UPART, FFORMA)
        implicit none
        ! Define types and parameters
        integer, intent(in) :: npe, npe_p
        integer, intent(in) :: gaprox
        real(8), intent(in) :: U(:)
        real, intent(in) :: PARTICLES(:,:)      ! Assuming shape (nparticles, 3)
        real(8), intent(in) :: KSIETA(:,:)         ! Assuming shape (nksieta, 3)
        real, intent(in) :: ELEMS(:,:)           ! Assuming shape (nelems, npe)
        
        real(8), intent(out) :: UPART(:)          ! To be allocated with appropriate size

        ! Local variables
        integer :: nparticles, nod, elem, node, current_node, part_i
        integer :: nksieta
        real(8) :: FFORMA(:,:)

        ! Temporary variables
        real(8), allocatable :: FI_nod(:,:)
        real(8), allocatable :: U_elem(:)
        real(8) :: ksi, eta
        real(8) :: result_fi(npe)
        integer :: line,column, dir, n, offset_u



        ! Remove last column of SIGMA_full and flatten to SIGMA_V
        ! allocate(SIGMA_M(size(SIGMA_full,1), size(SIGMA_full,2)-1))
        ! SIGMA_M = SIGMA_full(:,1:size(SIGMA_full,2)-1)
        ! allocate(SIGMA_V(size(SIGMA_M,1)*size(SIGMA_M,2)))
        ! SIGMA_V = reshape(transpose(SIGMA_M), shape(SIGMA_V))
        ! deallocate(SIGMA_M)

        ! Determine number of particles
        nparticles = size(PARTICLES, 1)
        nksieta = size(KSIETA, 1)

        ! Allocate UPART and SIGPART
        ! Assuming UPART_i and SIGPART_i are vectors of size 6*npe and 9*npe respectively


        ! Loop over each particle
        do part_i = 1, nparticles
            ! Initialize U_elem and SIG_elem

            allocate(U_elem(2*npe_p*npe))
            U_elem = 0.0d0

            ! Initialize FI matrices
            allocate(FI_nod(6, 2*npe_p*npe))
            FI_nod = 0.0d0

            ! Loop through each node of the particle (3 nodes)
            do nod = 1, 3
                offset_u = (nod - 1) * 2 * npe
                ! Get the corresponding element
                elem = int(KSIETA(int(PARTICLES(part_i, nod)) ,1))

                if (elem < 1 .or. elem > size(ELEMS,1)) then
                    print *, "Invalid element index for particle ", part_i, " node ", nod
                    stop
                endif

                ! Create displacement and stress vectors for the element
                do node = 1, npe
                    current_node = int(ELEMS(elem, node))

                    U_elem(offset_u+2*node-1:offset_u+2*node) = U(2*current_node-1:2*current_node)
                    

                end do

                ! Get ksi and eta coordinates
                ksi = KSIETA(int(PARTICLES(part_i, nod)), 2)
                eta = KSIETA(int(PARTICLES(part_i, nod)), 3)
                result_fi = fi(FFORMA, ksi, eta, gaprox)
                ! print *, 'result_fi:', result_fi
                ! stop
                ! Fill FI_nod matrix
                do dir = 1, 2
                    line = (nod - 1)*2 + dir
                    do n = 1, npe
                        column=(nod - 1)*2*npe + 2*(n-1) + dir
                        FI_nod(line, column) = result_fi(n)
                    end do
                end do


            end do


            UPART(6*(part_i-1)+1 : 6*part_i) = matmul(FI_nod, U_elem)


            ! SIGPART(part_i, :) = matmul(FI_sigma, SIG_elem)

            ! Deallocate temporary arrays
            deallocate(U_elem,  FI_nod)
        end do

    end subroutine particleDisplacement

    subroutine particleStress(PARTICLES,NODESP,U_PART,npart,npe_p,gaprox_p,FFORMA_P,SIGMA_P, ep_p)
        implicit none
        ! Define types and parameters
        integer, intent(in) :: npart, gaprox_p, npe_p, ep_p
        real(8), intent(in) :: U_PART(:)
        real, intent(in) :: PARTICLES(:,:)      ! Assuming shape (nparticles, 3)
        real(8), intent(in) :: NODESP(:,:)         ! Assuming shape (nnodes, 2)
        real(8), intent(inout) :: FFORMA_P(:,:)
        
        real(8), intent(out) :: SIGMA_P(:,:)        ! To be allocated with appropriate size

        integer :: iel, ih, ine,index
        real(8) :: ksi, eta, detjac
        real(8) ::  dxdksi, dxdeta, dydksi, dydeta, dudksi, dudeta, dvdksi, dvdeta, dksidx, dksidy, detadx, detady
        real(8) :: J(2,2), JINV(2,2)
        real(8) :: du(npe_p), dv(npe_p), result_fi_p(npe_p), result_dfidksi_p(npe_p), result_dfideta_p(npe_p)
        real(8) :: COORD_P(npe_p,2)
        real(8) :: cx(npe_p), cy(npe_p)
        real(8) :: d11, d12, d21, d22, d33, d44
        real(8) :: epsilon_x, epsilon_y, epsilon_xy, sigmax, sigmay, talxy
        real    :: h, bx, by
        ! Local 
        call shapeFunc(gaprox_p,npe_p,FFORMA_P,COORD_P)


        do iel=1,npart
            index=1
            do ine = 1, npe_p
                cx(ine) = NODESP(int(PARTICLES(iel, ine)), 1)
                cy(ine) = NODESP(int(PARTICLES(iel, ine)), 2)
                du(ine)=U_PART(2*int(PARTICLES(iel,ine)-1)+1)
                dv(ine)=U_PART(2*int(PARTICLES(iel,ine)-1)+2)
            end do
            
            do ih=1,npe_p
                ksi=COORD_P(ih,1)
                eta=COORD_P(ih,2)
    
                    
    
    
                result_fi_p=fi(FFORMA_P,ksi,eta,gaprox_p)
                result_dfidksi_p=dfidksi(FFORMA_P,ksi,eta,gaprox_p)
                result_dfideta_P=dfideta(FFORMA_P,ksi,eta,gaprox_p)
    
                dxdksi = dot_product(result_dfidksi_p, cx)
                dxdeta = dot_product(result_dfideta_p, cx)
                dydksi = dot_product(result_dfidksi_p, cy)
                dydeta = dot_product(result_dfideta_p, cy)
    
                dudksi = dot_product(result_dfidksi_p, du)
                dudeta = dot_product(result_dfideta_p, du)
                dvdksi = dot_product(result_dfidksi_p, dv)
                dvdeta = dot_product(result_dfideta_p, dv)
    
                J(1,1) = dxdksi
                J(1,2) = dxdeta
                J(2,1) = dydksi
                J(2,2) = dydeta
    
                ! call dgetrf(2, 2, J, 2, ipiv, info)
                ! call dgetri(2, J, 2, ipiv, work, lwork, info)
                ! JINV = J
    
    
                detjac = J(1,1)*J(2,2) - J(1,2)*J(2,1)
    
                JINV(1,1) =  J(2,2) / detjac
                JINV(1,2) = -J(1,2) / detjac
                JINV(2,1) = -J(2,1) / detjac
                JINV(2,2) =  J(1,1) / detjac
    
                dksidx=JINV(1,1)
                dksidy=JINV(1,2)
                detadx=JINV(2,1)
                detady=JINV(2,2)
    
                call materialProperties(PARTICLES, iel, npe_p, ep_p, d11, d12, d21, d22, d33, d44, h, bx, by)

                epsilon_x = dudksi*dksidx + dudeta*detadx
                epsilon_y = dvdksi*dksidy + dvdeta*detady
                epsilon_xy = 0.5*(dudksi*dksidy+dudeta*detady+dvdksi*dksidx+dvdeta*detadx)
    
                sigmax=d11*epsilon_x+d12*epsilon_y
                sigmay=d21*epsilon_x+d22*epsilon_y
                talxy=d33*epsilon_xy
                
                
    
                SIGMA_P(int(PARTICLES(iel,ih)),:) = SIGMA_P(int(PARTICLES(iel,ih)),:) + [sigmax, sigmay, talxy,1.0d0]
                ! stop
            end do
        end do 

        end subroutine particleStress
        

end module particle_displacement_module