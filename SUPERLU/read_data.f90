module read_data_module
    implicit none
contains
    subroutine read_input_file(filename,nnodes, nmats, nthick, nelem3, nelem6, nelem10, &
                                nelems, nforces, npressures, ndisplas, &
                                NODES, MATS, THICKS, ELEMS, LOAD, VINC, &
                                gaprox, npe)
        ! Declare arguments
        character(len=100), intent(in) :: filename
        integer, intent(out) :: nnodes, nmats, nthick, nelem3, nelem6, nelem10
        integer, intent(out) :: nelems, nforces, npressures, ndisplas, gaprox, npe
        real, allocatable :: MATS(:,:), THICKS(:), ELEMS(:,:), LOAD_aux(:,:), LOAD(:,:),VINC_aux(:,:) ,VINC(:,:)
        
        real(8), allocatable :: NODES(:,:)

        integer :: i, ierr,iel,ino, material_index,thick_index,component_i,components_count ,vinc_i,vinc_count
        
        real :: young, ni
        character(len=100) :: line  ! To read full lines
        character(len=10) :: type_DOF  ! To read full lines
        ! Open the input file
        open(unit=10, file=filename, status="old", action="read", iostat=ierr)
        if (ierr /= 0) then
            print *, "Error opening file"
            stop
        end if

        ! Read and extract the number of nodes
        read(10, *) line, nnodes
        read(10, *) line, nmats
        read(10, *) line, nthick

        read(10, *) line, nelem3
        read(10, *) line, nelem6
        read(10, *) line, nelem10

        nelems = nelem3 + nelem6 + nelem10

        read(10, *) line, nforces
        read(10, *) line, npressures
        read(10, *) line, ndisplas

        read(10, *) 
        read(10, *)

        allocate (NODES(nnodes, 2))

        do i = 1, nnodes
            read(10, *) line, NODES(i, 1), NODES(i, 2)
        end do

        read(10, *) 
        read(10, *)

        allocate (MATS(nmats, 2))

        do i = 1, nmats
            read(10, *) line, MATS(i, 1), MATS(i, 2)
        end do

        read(10, *) 
        read(10, *)

        allocate (THICKS(nthick))
        do i = 1, nthick
            read(10, *) line, THICKS(i)
        end do

        if (nelem3 > 0) then
            npe = 3
            gaprox = 1 
        end if

        if (nelem6 > 0) then
            npe = 6
            gaprox = 2
        end if

        if (nelem10 > 0) then
            npe = 10
            gaprox = 3
        end if

        read(10, *) 
        read(10, *)

        allocate (ELEMS(nelems, npe + 5))  ! Fixed allocate statement
        do iel = 1, nelems
            read(10, *) line, (ELEMS(iel, ino), ino = 1, npe + 2)  ! Explicitly declare ino
            ! YOUNG MODULUS
                   ! Store the material index for reuse
            material_index = int(ELEMS(iel, npe + 1))
            
            ! YOUNG MODULUS
            young = MATS(material_index, 1)
            ELEMS(iel, npe + 1) = young
            
            ! POISSON RATIO
            ni = MATS(material_index, 2)
            ELEMS(iel, npe + 2) = ni

            ! THICKNESS
            ! Weirdly, if I try to use THICKS(int(ELEMS(iel, npe + 3))), It doenst work
            thick_index = int(ELEMS(iel, npe + 3))
            ELEMS(iel, npe + 3) = THICKS(material_index)
            ! BODY FORCES
            ELEMS(iel, npe + 4) = 0.0
            ELEMS(iel, npe + 5) = 0.0
        end do
        
        read(10, *)
        read(10, *)
        
        allocate (LOAD_aux(nforces, 3))
        components_count=0
        do i = 1, nforces
            read(10, *) line, LOAD_aux(i, 1), LOAD_aux(i, 2), LOAD_aux(i, 3)
            ! Component X
            if (LOAD_aux(i, 2) /= 0.0) then
                components_count = components_count + 1
            end if
            ! Component Y
            if (LOAD_aux(i, 3) /= 0.0) then
                components_count = components_count + 1
            end if
        end do
        allocate (LOAD(components_count, 3))
        component_i=0
        do i = 1, nforces
            if (LOAD_aux(i, 2) /= 0.0) then
                component_i = component_i + 1
                LOAD(component_i, 1) = LOAD_aux(i, 1)
                LOAD(component_i, 2) = 1
                LOAD(component_i, 3) = LOAD_aux(i, 2)
            end if
            if (LOAD_aux(i, 3) /= 0.0) then
                component_i = component_i + 1
                LOAD(component_i, 1) = LOAD_aux(i, 1)
                LOAD(component_i, 2) = 2
                LOAD(component_i, 3) = LOAD_aux(i, 3)
            end if
        end do

        read(10, *)
        read(10, *)
        allocate (VINC_aux(ndisplas, 3))
        vinc_count=0
        do i = 1, ndisplas
            read(10, *) line, VINC_aux(i, 1), type_DOF, VINC_aux(i, 3)
            print *, "test: ", type_DOF
            ! Component X
            if (trim(type_DOF)=="BOTH") then
                VINC_aux(i, 2) = 3 ! 3 MEANS RESTRICTED BOTH X AND Y
                vinc_count = vinc_count + 2
            end if
            if (trim(type_DOF)=="X") then
                VINC_aux(i, 2) = 1 ! 1 MEANS RESTRICTED X
                vinc_count = vinc_count + 1
            end if
            if (trim(type_DOF)=="Y") then
                VINC_aux(i, 2) = 2 ! 2 MEANS RESTRICTED Y
                vinc_count = vinc_count + 1
            end if
        end do
        allocate (VINC(vinc_count, 4))
        vinc_i=0
        do i = 1, ndisplas
            if (VINC_aux(i, 2) == 1) then
                vinc_i = vinc_i + 1
                VINC(vinc_i, 1) = VINC_aux(i, 1)
                VINC(vinc_i, 2) = 1
                VINC(vinc_i, 3) = VINC_aux(i, 3)
                VINC(vinc_i, 4) = 1
            end if
            if (VINC_aux(i, 2) == 2) then
                vinc_i = vinc_i + 1
                VINC(vinc_i, 1) = VINC_aux(i, 1)
                VINC(vinc_i, 2) = 2
                VINC(vinc_i, 3) = VINC_aux(i, 3)
                VINC(vinc_i, 4) = 1
            end if
            if (VINC_aux(i, 2) == 3) then
                vinc_i = vinc_i + 1

                VINC(vinc_i, 1) = VINC_aux(i, 1)
                VINC(vinc_i, 2) = 1
                VINC(vinc_i, 3) = VINC_aux(i, 3)
                VINC(vinc_i, 4) = 1

                vinc_i = vinc_i + 1

                VINC(vinc_i, 1) = VINC_aux(i, 1)
                VINC(vinc_i, 2) = 2
                VINC(vinc_i, 3) = VINC_aux(i, 3)
                VINC(vinc_i, 4) = 1
            end if
        end do

        ! Close the file after reading all data
        close(10)
    end subroutine read_input_file



    !===========================================!
    ! Subroutine to compute material properties !
    !===========================================!
    subroutine materialProperties(ELEMS, elem, npe, ep, d11, d12, d21, d22, d33, d44, h, bx, by)
    implicit none
    ! Arguments
    real, intent(in) :: ELEMS(:,:)     ! 2D array containing element properties
    integer, intent(in) :: elem        ! Element index
    integer, intent(in) :: npe       ! Node index
    integer, intent(in) :: ep          ! Element property flag (0 or other)

    real(8), intent(out) :: d11, d12, d21, d22, d33, d44
    real, intent(out) :: h, bx, by     ! Output variables: thickness, body forces bx and by

    real :: young, ni  ! Young's modulus, Poisson's ratio

    ! Extract material properties from the ELEMS matrix for the specified element and node
    young = ELEMS(elem, npe + 1)
    ni = ELEMS(elem, npe + 2)
    h = ELEMS(elem, npe + 3)
    bx = ELEMS(elem, npe + 4)         ! Body force bx
    by = ELEMS(elem, npe + 5)         ! Body force by
    ! print *, "Young's modulus = ",young,"element",elem
    ! print *, "Young's modulus = ",young,"element",elem
    ! print *, "Poisson's ratio = ", ni
    ! print *, "Thickness = ", h
    ! Compute the material properties based on ep value
    if (ep == 0) then
        d11 = young / (1.0 - ni**2)
        d12 = d11 * ni
        d21 = d12
        d22 = d11
        d33 = d11 * (1.0 - ni)/2
        d44 = d33
    else
        d11 = (1.0 - ni) * young / ((1.0 + ni) * (1.0 - 2.0 * ni)) 
        d12 = d11 * ni
        d21 = d12
        d22 = d11
        d33 = d11 * (1.0 - ni)
        d44 = d33
    end if
end subroutine materialProperties

end module read_data_module
