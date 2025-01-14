module writing_acadview_module
    implicit none

    contains

    subroutine appending(SIGMA, ELEMS, NODES, U, SIGMA_P, PARTICLES, NODESP, U_PART, nelems, nnodes, npe, gaprox)
        implicit none

        ! Input parameters
        integer, intent(in) :: nelems, nnodes, npe, gaprox
        real(8), intent(in) :: SIGMA(:,:)
        real(8), intent(in) :: SIGMA_P(:,:)
        real(8), intent(in) :: U(:)
        real(8), intent(in) :: U_PART(:)
        real,intent(in) ::  ELEMS(:,:)
        real,intent(in) :: PARTICLES(:,:)
        
        real(8), intent(in) :: NODES(:,:)
        real(8), intent(in) :: NODESP(:,:)

        
        ! Local variables
        integer :: i, j,total_nodes,total_elements
        
        ! File handling
        integer :: unit_number, io_status
        character(len=100) :: path
        integer,allocatable :: list(:)
        

        total_nodes=nnodes+size(NODESP,1)
        total_elements=nelems+size(PARTICLES,1)

        ! Open a file for writing
        path = 'output.ogl'
        open(newunit=unit_number, file=path, status='replace', action='write', iostat=io_status)
        if (io_status /= 0) then
            print *, "Error opening file."
            stop
        end if

        ! Write header information
        write(unit_number, '(A)') 'arquivo de entrada'
        write(unit_number, '(A)') 'n.nos n.elems n.listas'
        write(unit_number, '(A)') '#'
        write(unit_number, '(I10, I10, I10)') total_nodes , total_elements , 5
        write(unit_number, '(A)') 'coordx coordy coordz deslx desly deslz'
        write(unit_number, '(A)') '#'

        ! Write NODES
        do i = 1, nnodes
            write(unit_number, '(F18.8 , F18.8 , F18.8 , F18.8 , F18.8 , F18.8)') &
                NODES(i, 1), NODES(i, 2), 0.0, 0.0, 0.0, 0.0
        end do

                ! Write NODESP (Particle Nodes)
        do i = 1, size(NODESP, 1)
            write(unit_number, '(F18.8, F18.8, F18.8, F18.8, F18.8, F18.8)') &
                NODESP(i, 1), NODESP(i, 2), 0.0d0, 0.0d0, 0.0d0, 0.0d0
        end do


        ! Write element information
        write(unit_number, '(A)') 'tpelem (1-barra/2-triang/3-quad) grauaprox nó1 nó2...nó_n group'
        write(unit_number, '(A)') '#'

        ! Write ELEMS
        do i = 1, nelems
            allocate(list(npe + 3))
            list(1) = 2
            list(2) = gaprox
            do j = 1, npe
                list(j + 2) = int(ELEMS(i, j))
            end do
            list(npe + 3) = 0
            write(unit_number, '(13I10)') list
            deallocate(list)
        end do


                ! Write PARTICLES
        do i = 1, size(PARTICLES, 1)
            allocate(list(3 + 3))
            list(1) = 2                                 ! Element type
            list(2) = 1                                 ! Approximation degree for particles
            do j = 1, 3
                list(j + 2) = int(PARTICLES(i, j)) + nnodes
            end do
            list(6) = 1                                 ! Group
            write(unit_number, '(6I10)') list
            deallocate(list)
        end do


        ! Write lists
        write(unit_number, '(A)') 'listas'
        write(unit_number, '(A)') '#'
        write(unit_number, '(A)') 'desl.x'

        ! Write displacement in x-direction
        do i = 1, nnodes
            write(unit_number, '(4F18.8)') &
                U(2 * i - 1), U(2 * i), 0.0, U(2 * i - 1)
        end do


        ! print *, 'U_PART:', U_PART
        ! STOP
        do i = 1, size(NODESP, 1)
            write(unit_number, '(4F18.8)') &
                U_PART(2 * i - 1), U_PART(2 * i), 0.0d0, U_PART(2 * i - 1)
        end do


        ! Write displacement in y-direction
        write(unit_number, '(A)') '#'
        write(unit_number, '(A)') 'desl.y'
        do i = 1, nnodes
            write(unit_number, '(4F18.8)') &
                U(2 * i - 1), U(2 * i), 0.0, U(2 * i)
        end do
        

        do i = 1, size(NODESP, 1)
            write(unit_number, '(4F18.8)') &
                U_PART(2 * i - 1), U_PART(2 * i), 0.0d0, U_PART(2 * i)
        end do

        ! Write sigma.x
        write(unit_number, '(a)') '#'
        write(unit_number, '(a)') 'sigma.x'
        do i = 1, nnodes
            write(unit_number, '(4F18.8)') &
                U(2 * i - 1), U(2 * i), 0.0, SIGMA(i, 1) / SIGMA(i, 4)
        end do

        do i = 1, size(NODESP, 1)
            write(unit_number, '(4F18.8)') &
                U_PART(2 * i - 1), U_PART(2 * i), 0.0d0, SIGMA_P(i, 1)/SIGMA_P(i,4)
        end do


        ! Write sigma.y
        write(unit_number, '(a)') '#'
        write(unit_number, '(a)') 'sigma.y'
        do i = 1, nnodes
            write(unit_number, '(4F18.8)') &
                U(2 * i - 1), U(2 * i), 0.0, SIGMA(i, 2) / SIGMA(i, 4)
        end do

        do i = 1, size(NODESP, 1)
            write(unit_number, '(4F18.8)') &
                U_PART(2 * i - 1), U_PART(2 * i), 0.0d0, SIGMA_P(i, 2)/SIGMA_P(i,4)
        end do
        

        ! Write tal.xy
        write(unit_number, '(a)') '#'
        write(unit_number, '(a)') 'tal.xy'
        do i = 1, nnodes
            write(unit_number, '(4F18.8)') &
                U(2 * i - 1), U(2 * i), 0.0, SIGMA(i, 3) / SIGMA(i, 4)
        end do


        do i = 1, size(NODESP, 1)
            write(unit_number, '(4F18.8)') &
                U_PART(2 * i - 1), U_PART(2 * i), 0.0d0, SIGMA_P(i, 3)/SIGMA_P(i,4)
        end do
        

        ! Close the file
        close(unit_number)
        print *, "File '", trim(path), "' created successfully."

end subroutine appending
    
end module writing_acadview_module