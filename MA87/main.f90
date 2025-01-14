program main
    use read_data_module  ! Import the module containing the subroutine
    use hammer_module
    use shape_module
    use assembly_module
    use particle_generator_module
    use sigma_module
    use particle_displacement_module
    use writing_acadview_module
    use sparse

    implicit none
    !================================================!
    !================= READING DATA =================!
    !================================================!
    character(len=100) :: input_file
    integer :: nnodes,nmats,nthick, nelem3, nelem6, nelem10, nelems,nforces,npressures,ndisplas, gaprox, npe
    real, allocatable :: MATS(:,:) , THICKS(:), ELEMS(:,:), LOAD(:,:), VINC(:,:)
    real(8),allocatable :: NODES(:,:)

    !==============================================!
    !============== MATERIAL PROPERTIES ===========!
    !==============================================!

    integer :: ep

    !================================================!
    !================= HAMMER POINTS ================!
    !================================================!

    real(8) :: hammer_points(3, 12)
    integer :: nph

    !================================================!
    !=============== SHAPE FUNCTIONS ================!
    !================================================!

    real(8), dimension(:,:), allocatable :: FFORMA, COORD


    !================================================!
    !================== PARTICLES ===================!
    !================================================!
    ! real:: radius,fv,young_p,ni_p,h_p,bx_p,by_p
    real,allocatable:: PARTICLES(:,:)

    real(8),allocatable:: NODESP(:,:),KSIETA(:,:)

    real(8),dimension(:,:),allocatable:: FFORMA_P,COORD_P
    integer :: npart, npe_p,gaprox_p,ep_p
    character(len=100) :: particle_input_file

    !================================================!
    !=================== ASSEMBLY  ==================!
    !================================================!
    real(8),dimension(:),allocatable::FGLOBAL
    type(sparse_matrix):: KGLOBAL_sp
    ! real(8) :: sum_KGLOBAL

    real(8), allocatable :: U(:)

    !================================================!
    !================ STRESS CALCULATION ============!
    !================================================!
    real(8),dimension(:,:),allocatable:: SIGMA


    !================================================!
    !================ PARTICLE DISPLACEMENT =========!
    !================================================!
    real(8),dimension(:,:),allocatable::  SIGMA_P
    real(8),dimension(:),allocatable:: UPART 
    !========================================================================================================!
    !__________________________________________|==  NON-OPTIMIZED ==|________________________________________!
    !__________________________________________|= COMPILE COMMAND  =|________________________________________!
    !__________________________________________|====================|________________________________________!
    !_    gfortran -fcheck=all -g -fopenmp -ffree-line-length-none main.f90 read_data.f90 hammer.f90 shape.f90 particle_generator.f90 assembly.f90 sigma.f90 particle_displacement.f90 writing_acadview.f90 Dependent_BLAS_package.for HSL_Dependencies_01.for MA87_common90.f90 MA87_ddeps90.f90 MA87_fakemetis.f MA87d.f90 MA87_INTERFACE.f90 sparse_set.f90 -o program.out -lopenblas -llapack -lblas          ___!
    !========================================================================================================!
    !========================================================================================================!


    !========================================================================================================!
    !__________________________________________|====  OPTIMIZED ====|________________________________________!
    !__________________________________________|= COMPILE COMMAND  =|________________________________________!
    !__________________________________________|====================|________________________________________!
    !_    gfortran -O3 -ffast-math -fopenmp -ffree-line-length-none main.f90 read_data.f90 hammer.f90 shape.f90 particle_generator.f90 assembly.f90 sigma.f90 particle_displacement.f90 writing_acadview.f90 Dependent_BLAS_package.for HSL_Dependencies_01.for MA87_common90.f90 MA87_ddeps90.f90 MA87_fakemetis.f MA87d.f90 MA87_INTERFACE.f90 sparse_set.f90 -o program.out -lopenblas -llapack -lblas         _____!
    !========================================================================================================!
    !========================================================================================================!

    ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||!

    !========================================================================================================!
    !___________________________________|==== FIRST TIME COMPILING  ====|____________________________________!
    !========================================================================================================!

    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none Dependent_BLAS_package.for
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none HSL_Dependencies_01.for
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none MA87_common90.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none MA87_ddeps90.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none MA87_fakemetis.f
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none MA87d.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none MA87_INTERFACE.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none sparse_set.f90

    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none main.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none read_data.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none hammer.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none shape.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none particle_generator.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none assembly.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none sigma.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none particle_displacement.f90
    !       gfortran -c -fcheck=all -g -fopenmp -ffree-line-length-none writing_acadview.f90

    !========================================================================================================!
    !========================================================================================================!
    


    ! Call the subroutine to read the input file and store values in the declared variables
    ep  = 0
    input_file = "fi1.txt"
    call read_input_file(input_file,nnodes, nmats, nthick, nelem3, nelem6, nelem10, &
                                nelems, nforces, npressures, ndisplas, &
                                NODES, MATS, THICKS, ELEMS, LOAD, VINC, &
                                gaprox, npe)


    call hammer(nph, hammer_points)



    ! Allocate FFORMA based on your needs (example: 3xsize_)
    allocate(FFORMA(npe,npe))
    allocate(COORD(npe,2))
    call shapeFunc(gaprox, npe, FFORMA, COORD)
    


    ! young_p=1000.0
    ! ni_p=0.25
    ! h_p=1.0
    ! bx_p=0.0
    ! by_p=0.0
    ! radius=0.1
    ! fv=0.1
    gaprox_p=1
    npe_p=3
    particle_input_file="pf1.txt"
    ! call generateParticles(NODES,ELEMS, npe, gaprox,radius,fv,young_p,ni_p,h_p,bx_p,by_p,NODESP,PARTICLES,KSIETA)

    call readParticles(particle_input_file, PARTICLES, NODESP, KSIETA)
    npart=size(PARTICLES,1)

    
    allocate(FGLOBAL(2*nnodes))
    
    !=========================================================!
    !=========== PREPARING TO USE SPARSE MATRIX ==============!
    !=========================================================!
    !|                                                       |!
    ! CALL prepare_to_use(KGLOBAL_sp,2*nnodes,((2*npe*npe+npe)*nelems+(((2*npe_p*npe)**2+2*npe_p*npe)/2)*npart))  !
    CALL prepare_to_use(KGLOBAL_sp,2*nnodes,((2*npe*npe+npe)*nelems+((2*npe_p*npe)**2)*npart)) !
    print*,((2*npe*npe+npe)*nelems)
    print*,(((2*npe_p*npe)**2)*npart)
    ! stop


    !|                                                       |!
    !=========================================================!

    call stiffnessMatrix(ELEMS, NODES, LOAD, FFORMA, gaprox, npe, nelems, nnodes, nph, hammer_points, ep, FGLOBAL, KGLOBAL_sp)

    print *, "======================================="
    Print *, "| GLOBAL STIFFNESS OF ELEMENTS DONE ! |"
    print *, "======================================="

    allocate(FFORMA_P(npe_p,npe_p))
    allocate(COORD_P(npe_p,2))
    call shapeFunc(gaprox_p, npe_p, FFORMA_P, COORD_P)


    call particleStiffnessMatrix (KGLOBAL_sp,FGLOBAL,ELEMS,PARTICLES, NODESP, KSIETA, VINC, FFORMA,FFORMA_P,gaprox, &
    gaprox_p,npe,npe_p,npart, nph, hammer_points, ep, nnodes)

    print *, "======================================"
    Print *, "|  CONTRIBUTION OF PARTICLES DONE !  |"
    print *, "======================================"


    allocate(U(2*nnodes))

    call solve_system_of_equation(KGLOBAL_sp, FGLOBAL, U)

    print *, "===================================="
    Print *, "|   SYSTEM OF EQUATIONS SOLVED!    |"
    print *, "===================================="

    ! ! Allocate SIGMA based on your needs (example: 3x3)
    allocate(SIGMA(nnodes,4))
    ! ! Initialize SIGMA to zero
    SIGMA = 0.0
    ! ! Call the subroutine to calculate the stresses
    call stressPiola(ELEMS,NODES,U,nelems,npe,gaprox,ep,FFORMA,SIGMA)

    print *, "===================================="
    print *, "|       STRESSES CALCULATED !        |"
    print *, "===================================="


    allocate(UPART(npart*2*npe_p))
    ! allocate(SIGPART(npart, 3*npe_p))
    UPART = 0.0
    ! SIGPART = 0.0

    call particleDisplacement(U, PARTICLES, KSIETA, ELEMS, gaprox, npe, npe_p, UPART, FFORMA)

    allocate(SIGMA_P(npart*3, 4))
    SIGMA_P = 0.0
    ep_p=0
    call particleStress(PARTICLES,NODESP,UPART,npart,npe_p,gaprox_p,FFORMA_P,SIGMA_P,ep_p)

    print *, "============================================"
    Print *, "|   PARTICLE DISPLACEMENTS CALCULATED !    |"
    print *, "============================================"


    ! allocate(SIGMA_P(npart*3, 3))
    ! SIGMA_P = 0.0
    ! do i = 1, npart
    !     do j = 1, 3
    !         SIGMA_P((i - 1) * 3 + j, :) = SIGPART(i, (j - 1) * 3 + 1:(j - 1) * 3 + 3)
    !     end do
    ! end do
    ! deallocate(SIGPART)


    call appending(SIGMA, ELEMS, NODES, U, SIGMA_P, PARTICLES, NODESP, UPART, nelems, nnodes, npe, gaprox)
    print *, "===================================="
    Print *, "|      OUTPUT FILE CREATED !       |"
    print *, "===================================="


    print *, "===================================="
    Print *, "|========    FINISHED !    ========|"
    print *, "===================================="


end program main
