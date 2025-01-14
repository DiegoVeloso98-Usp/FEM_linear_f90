module assembly_module
    use shape_module
    use read_data_module

    !_________________!
    !=================!
        USE SPARSE
    !=================!
    !_________________!


    implicit none


contains

subroutine stiffnessMatrix(ELEMS, NODES, LOAD, FFORMA, gaprox, npe, nelems, nnodes, nph, hammer, ep, &
FGLOBAL, KGLOBAL_sp)
    implicit none


    integer, intent(in) :: npe, nelems, nnodes, nph, gaprox,ep
    real(8), intent(inout) ::  hammer(3, nph), FFORMA(:,:)
    real, intent(in) :: ELEMS(nelems, npe+5), LOAD(:,:)   !, VINC(:,:)

    real(8),  intent(in) :: NODES(:,:)

    real(8), intent(out) ::  FGLOBAL(2*nnodes)

    !______________________________________!
    !                                      !
       type(sparse_matrix) :: KGLOBAL_sp   !
       integer::indexes_sparse(2*npe)      !
    !                                      !
    !______________________________________!

    ! Local variables
    integer :: iel, ih, ine, idir, in_
    real(8) :: ksi, eta, peso, detjac,  dxdksi, dxdeta, dydksi, dydeta
    real(8) :: J(2,2), JINV(2,2), DX(2,4), DY(2,4), DFI(4, 2*npe) , MXX(2,2), MXY(2,2), MYX(2,2), MYY(2,2)
    real(8) :: cx(npe), cy(npe), KLOCAL(2*npe, 2*npe), FLOCAL(2*npe)
    real(8) :: d11, d12, d21, d22, d33, d44
    real:: h, bx, by
    real(8), dimension(4, 4) :: temp1, temp2, temp3, temp4, sum_terms

    real(8)::result_fi(npe),result_dfidksi(npe),result_dfideta(npe)
    real(8), allocatable :: MATRIX_result_fi(:,:), MATRIX_result_dfidksi(:,:), MATRIX_result_dfideta(:,:)


    ! ! Initialize global matrices
    FGLOBAL = 0.0d0

    ! Shape functions don't change in the loop, so initialize them once
    allocate(MATRIX_result_fi(nph, npe))
    allocate(MATRIX_result_dfidksi(nph, npe))
    allocate(MATRIX_result_dfideta(nph, npe))
    MATRIX_result_fi = 0.0d0
    MATRIX_result_dfidksi = 0.0d0
    MATRIX_result_dfideta = 0.0d0
    do ih = 1, nph
        ksi = hammer(1, ih)
        eta = hammer(2, ih)
        MATRIX_result_fi(ih, :) = fi(FFORMA, ksi, eta, gaprox)
        MATRIX_result_dfidksi(ih, :) = dfidksi(FFORMA, ksi, eta, gaprox)
        MATRIX_result_dfideta(ih, :) = dfideta(FFORMA, ksi, eta, gaprox)
    end do

    ! Loop over elements
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(KGLOBAL_sp,FGLOBAL,NODES,ELEMS,FFORMA,hammer,ep,gaprox,nph,npe,nelems,nnodes,MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta)   
    !$OMP DO
    do iel = 1, nelems
        print *, "Element: ", iel

        call materialProperties(ELEMS, iel, npe, ep, d11, d12, d21, d22, d33, d44, h, bx, by)
        MXX(1, 1) = d11
        MXX(1, 2) = 0
        MXX(2, 1) = 0
        MXX(2, 2) = d33

        MXY(1, 1) = 0
        MXY(1, 2) = d12
        MXY(2, 1) = d33
        MXY(2, 2) = 0

        MYX(1, 1) = 0
        MYX(1, 2) = d44
        MYX(2, 1) = d21
        MYX(2, 2) = 0

        MYY(1, 1) = d44
        MYY(1, 2) = 0
        MYY(2, 1) = 0
        MYY(2, 2) = d22
        
        
        ! Extract element node coordinates
        do ine = 1, npe
            cx(ine) = NODES(int(ELEMS(iel, ine)), 1)
            cy(ine) = NODES(int(ELEMS(iel, ine)), 2)
        end do

        ! Loop over integration points (nph)
        KLOCAL = 0.0d0
        do ih = 1, nph
            peso = hammer(3, ih)

            result_fi      = MATRIX_result_fi(ih, :)
            result_dfidksi = MATRIX_result_dfidksi(ih, :)
            result_dfideta = MATRIX_result_dfideta(ih, :)

            ! Calculate Jacobian matrix and its determinant
            dxdksi = dot_product(result_dfidksi, cx)
            dxdeta = dot_product(result_dfideta, cx)
            dydksi = dot_product(result_dfidksi, cy)
            dydeta = dot_product(result_dfideta, cy)
            

            J(1,1) = dxdksi
            J(1,2) = dxdeta
            J(2,1) = dydksi
            J(2,2) = dydeta

            ! Compute determinant and inverse of Jacobian matrix
            detjac = J(1,1)*J(2,2) - J(1,2)*J(2,1)

            
            if (detjac == 0.0d0) then
                print *, "Warning: Singular Jacobian in element ", iel
                stop
            endif

            ! Compute inverse of Jacobian directly (optimized)
            JINV(1,1) =  J(2,2) / detjac
            JINV(1,2) = -J(1,2) / detjac
            JINV(2,1) = -J(2,1) / detjac
            JINV(2,2) =  J(1,1) / detjac


            ! Partial derivatives of ksi, eta with respect to x, y
            DX(1, 1) = JINV(1, 1)
            DX(1, 2) = 0.0d0
            DX(1, 3) = JINV(2,1)
            DX(1, 4) = 0.0d0
            DX(2, 1) = 0.0d0
            DX(2, 2) = JINV(1, 1)
            DX(2, 3) = 0.0d0
            DX(2, 4) = JINV(2, 1)

            DY(1, 1) = JINV(1,2)
            DY(1, 2) = 0.0d0
            DY(1, 3) = JINV(2, 2)
            DY(1, 4) = 0.0d0
            DY(2, 1) = 0.0d0
            DY(2, 2) = JINV(1,2)
            DY(2, 3) = 0.0d0
            DY(2, 4) = JINV(2, 2)

            DFI = 0.0d0
            do ine = 1, npe
                DFI(1, 2*ine-1) = result_dfidksi(ine)
                DFI(2, 2*ine)   = result_dfidksi(ine)
                DFI(3, 2*ine-1) = result_dfideta(ine)
                DFI(4, 2*ine)   = result_dfideta(ine)
            end do


            FLOCAL = 0.0d0

            ! Compute nodal forces and local stiffness matrix
            do ine = 1, npe
                FLOCAL(2*ine-1) = h * bx * detjac * peso * result_fi(ine)
                FLOCAL(2*ine)   = h * by * detjac * peso * result_fi(ine)
            end do


        ! Step-by-step multiplication of terms
            temp1 = matmul(transpose(DX), matmul(MXX, DX))           ! DX.T @ MXX @ DX
            temp2 = matmul(transpose(DX), matmul(MXY, DY))           ! DX.T @ MXY @ DY
            temp3 = matmul(transpose(DY), matmul(MYX, DX))           ! DY.T @ MYX @ DX
            temp4 = matmul(transpose(DY), matmul(MYY, DY))           ! DY.T @ MYY @ DY

            ! Summing up the terms
            sum_terms = temp1 + temp2 + temp3 + temp4


            KLOCAL = KLOCAL+ h * matmul(transpose(DFI), matmul(sum_terms, DFI)) * detjac * peso


            ! Assemble into global matrices
            indexes_sparse=0
            !$OMP CRITICAL
            do in_ = 1, npe
                do idir = 1, 2
                    FGLOBAL(2*(int(ELEMS(iel, in_))-1)+idir) = FGLOBAL(2*(int(ELEMS(iel, in_))-1)+idir) + FLOCAL(2*(in_-1)+idir)
                    indexes_sparse(2*(in_-1)+idir)=2*(int(ELEMS(iel, in_))-1)+idir

                    if (indexes_sparse(2*(in_-1)+idir) > 704880) then
                        write(*,*)
                        write(*,*) "ERRO"
                        write(*,*) iel,indexes_sparse(2*(in_-1)+idir),in_,idir
                        stop
                    endif
                end do
            end do
            !$OMP END CRITICAL

        end do


        !$OMP CRITICAL
        call add_matrix(KGLOBAL_sp,KLOCAL,indexes_sparse,2*npe)
        !$OMP END CRITICAL
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    ! call assemble_sparse_matrix(KGLOBAL_sp,timeit=.false.)
    
    ! Apply boundary conditions and forces
    call globalForce(LOAD, FGLOBAL)
    ! call boundaryCondition(VINC, KGLOBAL_sp, FGLOBAL,nnodes)


end subroutine stiffnessMatrix


subroutine particleStiffnessMatrix&
(KGLOBAL_sp,FGLOBAL,ELEMS,PARTICLES, NODESP, KSIETA, VINC, FFORMA,FFORMA_P,gaprox, &
gaprox_p,npe,npe_p,npart, nph, hammer, ep,nnodes)
    implicit none

    integer, intent(in) :: npe, nph, gaprox,ep,nnodes
    integer, intent(in) :: gaprox_p, npe_p, npart
    real(8), intent(inout) ::  hammer(3, nph), FFORMA(:,:),FFORMA_P(:,:)
    real, intent(in) :: ELEMS(:,:), VINC(:,:)
    real, intent(in) :: PARTICLES(:,:)

    real(8), intent(in) :: NODESP(:,:), KSIETA(:,:)

    real(8), intent(inout) :: FGLOBAL(:)

    !______________________________________ !
    !                                       !
       type(sparse_matrix) :: KGLOBAL_sp    !
       integer::indexes_sparse(2*npe_p*npe) ! [==  SIZE OF MATRIX KG  ==]
    !                                       !
    !______________________________________ !

    ! Local variables
    integer :: iel, ih, ine,index,i,j_,k
    real(8) :: ksi, eta, peso, detjac,  dxdksi, dxdeta, dydksi, dydeta
    real(8) :: J(2,2), JINV(2,2), DX(2,4), DY(2,4), DFI(4, 2*npe_p) , MXX(2,2), MXY(2,2), MYX(2,2), MYY(2,2)
    real(8) :: cx(npe_p), cy(npe_p), KLOCAL(2*npe_p, 2*npe_p), FLOCAL(2*npe_p)
    real(8) :: d11, d12, d21, d22, d33, d44
    real:: h, bx, by

    real(8), dimension(4, 4) :: temp1, temp2, temp3, temp4, sum_terms

    real(8),dimension(npe_p) ::result_fi_p,result_dfidksi_p,result_dfideta_p
    real(8), allocatable :: MATRIX_result_fi(:,:), MATRIX_result_dfidksi(:,:), MATRIX_result_dfideta(:,:)

    real(8)::vfi1(npe),vfi2(npe),vfi3(npe)
    real(8), dimension(2*npe*npe_p,2*npe*npe_p)::KG
    real(8), dimension(6,3*2*npe)::MFI
    real(8)::ksi1,ksi2,ksi3,eta1,eta2,eta3




        ! Shape functions don't change in the loop, so initialize them once
    allocate(MATRIX_result_fi(nph, npe_p))
    allocate(MATRIX_result_dfidksi(nph, npe_p))
    allocate(MATRIX_result_dfideta(nph, npe_p))
    MATRIX_result_fi = 0.0d0
    MATRIX_result_dfidksi = 0.0d0
    MATRIX_result_dfideta = 0.0d0
    do ih = 1, nph
        ksi = hammer(1, ih)
        eta = hammer(2, ih)
        MATRIX_result_fi(ih, :) = fi(FFORMA_P, ksi, eta, gaprox_p)
        MATRIX_result_dfidksi(ih, :) = dfidksi(FFORMA_P, ksi, eta, gaprox_p)
        MATRIX_result_dfideta(ih, :) = dfideta(FFORMA_P, ksi, eta, gaprox_p)
    end do


    do iel = 1, npart
        Print *, "Particle: ", iel
        call materialProperties(PARTICLES, iel, npe_p, ep, d11, d12, d21, d22, d33, d44, h, bx, by)
        MXX(1, 1) = d11
        MXX(1, 2) = 0
        MXX(2, 1) = 0
        MXX(2, 2) = d33

        MXY(1, 1) = 0
        MXY(1, 2) = d12
        MXY(2, 1) = d33
        MXY(2, 2) = 0

        MYX(1, 1) = 0
        MYX(1, 2) = d44
        MYX(2, 1) = d21
        MYX(2, 2) = 0

        MYY(1, 1) = d44
        MYY(1, 2) = 0
        MYY(2, 1) = 0
        MYY(2, 2) = d22


        ! Extract element node coordinates
        do ine = 1, npe_p
            cx(ine) = NODESP(int(PARTICLES(iel, ine)), 1)
            cy(ine) = NODESP(int(PARTICLES(iel, ine)), 2)
        end do


        KG = 0.0d0
        do ih = 1, nph

            peso = hammer(3, ih)

            result_fi_p      = MATRIX_result_fi(ih, :)
            result_dfidksi_p = MATRIX_result_dfidksi(ih, :)
            result_dfideta_p = MATRIX_result_dfideta(ih, :)


            ! Calculate Jacobian matrix and its determinant
            dxdksi = dot_product(result_dfidksi_p, cx)
            dxdeta = dot_product(result_dfideta_p, cx)
            dydksi = dot_product(result_dfidksi_p, cy)
            dydeta = dot_product(result_dfideta_p, cy)
            

            J(1,1) = dxdksi
            J(1,2) = dxdeta
            J(2,1) = dydksi
            J(2,2) = dydeta

            ! Compute determinant and inverse of Jacobian matrix
            detjac = J(1,1)*J(2,2) - J(1,2)*J(2,1)

            
            ! Print *, "Testing 3 !"
            if (detjac == 0.0d0) then
                print *, "Warning: Singular Jacobian in element ", iel
                stop
            endif

            ! Inverse of Jacobian matrix
            ! Compute inverse of Jacobian directly (optimized)
            JINV(1,1) =  J(2,2) / detjac
            JINV(1,2) = -J(1,2) / detjac
            JINV(2,1) = -J(2,1) / detjac
            JINV(2,2) =  J(1,1) / detjac

            ! Partial derivatives of ksi, eta with respect to x, y
            DX(1, 1) = JINV(1, 1)
            DX(1, 2) = 0.0d0
            DX(1, 3) = JINV(2,1)
            DX(1, 4) = 0.0d0
            DX(2, 1) = 0.0d0
            DX(2, 2) = JINV(1, 1)
            DX(2, 3) = 0.0d0
            DX(2, 4) = JINV(2, 1)

            DY(1, 1) = JINV(1,2)
            DY(1, 2) = 0.0d0
            DY(1, 3) = JINV(2, 2)
            DY(1, 4) = 0.0d0
            DY(2, 1) = 0.0d0
            DY(2, 2) = JINV(1,2)
            DY(2, 3) = 0.0d0
            DY(2, 4) = JINV(2, 2)


            DFI = 0.0d0
            do ine = 1, npe_p
                DFI(1, 2*ine-1) = result_dfidksi_p(ine)
                DFI(2, 2*ine)   = result_dfidksi_p(ine)
                DFI(3, 2*ine-1) = result_dfideta_p(ine)
                DFI(4, 2*ine)   = result_dfideta_p(ine)
            end do



            ! Local stiffness matrix
            KLOCAL = 0.0d0
            FLOCAL = 0.0d0

            ! Compute nodal forces and local stiffness matrix
            do ine = 1, npe_p
                FLOCAL(2*ine-1) = h * bx * detjac * peso * result_fi_p(ine)
                FLOCAL(2*ine)   = h * by * detjac * peso * result_fi_p(ine)
            end do



            ! Step-by-step multiplication of terms
            temp1 = matmul(transpose(DX), matmul(MXX, DX))           ! DX.T @ MXX @ DX
            temp2 = matmul(transpose(DX), matmul(MXY, DY))           ! DX.T @ MXY @ DY
            temp3 = matmul(transpose(DY), matmul(MYX, DX))           ! DY.T @ MYX @ DX
            temp4 = matmul(transpose(DY), matmul(MYY, DY))           ! DY.T @ MYY @ DY
            
            ! Summing up the terms

            sum_terms = temp1 + temp2 + temp3 + temp4

            KLOCAL = h * matmul(transpose(DFI), matmul(sum_terms, DFI)) * detjac * peso
            
            ksi1=KSIETA(int(PARTICLES(iel,1)),2)
            ksi2=KSIETA(int(PARTICLES(iel,2)),2)
            ksi3=KSIETA(int(PARTICLES(iel,3)),2)

            eta1=KSIETA(int(PARTICLES(iel,1)),3)
            eta2=KSIETA(int(PARTICLES(iel,2)),3)
            eta3=KSIETA(int(PARTICLES(iel,3)),3)

            vfi1=fi(FFORMA,ksi1,eta1,gaprox)
            vfi2=fi(FFORMA,ksi2,eta2,gaprox)
            vfi3=fi(FFORMA,ksi3,eta3,gaprox)

            index=1
            do i=1,npe_p
                do j_=1,npe
                    do k=1,2

                        indexes_sparse(index)=(2*int(ELEMS(int(KSIETA(int(PARTICLES(iel,i)),1)),j_)-1)+k)
                        if (indexes_sparse(index) > 704880) then
                            write(*,*) iel,indexes_sparse(index),i,j_,k
                            stop
                        endif
                        index=index+1
                    end do
                end do
            end do

            do i=1,npe
                MFI(1,2*(i-1)+1)=vfi1(i)
                MFI(2,2*(i-1)+2)=vfi1(i)
                MFI(3,2*npe+2*(i-1)+1)=vfi2(i)
                MFI(4,2*npe+2*(i-1)+2)=vfi2(i)
                MFI(5,4*npe+2*(i-1)+1)=vfi3(i)
                MFI(6,4*npe+2*(i-1)+2)=vfi3(i)
            end do

            KG=KG+matmul(transpose(MFI),matmul(KLOCAL,MFI))

        
        end do
        call add_full_matrix(KGLOBAL_sp,KG,indexes_sparse,2*npe_p*npe)
    end do

    call assemble_sparse_matrix(KGLOBAL_sp,timeit=.false.)


    ! Apply boundary conditions and forces
    call boundaryCondition(VINC, KGLOBAL_sp, FGLOBAL,nnodes)

end subroutine particleStiffnessMatrix

subroutine globalForce(LOAD, FGLOBAL)
    implicit none
    real, intent(in) :: LOAD(:,:)
    real(8), intent(inout) :: FGLOBAL(:)
    integer :: i, ino, idir
    real:: value
    print*, "Applying global forces"
    do i = 1, size(LOAD, 1)
        ino=int(LOAD(i, 1))
        idir=int(LOAD(i, 2))
        value=LOAD(i, 3)
        FGLOBAL(2*(ino-1)+idir)=FGLOBAL(2*(ino-1)+idir)+value
    end do
end subroutine globalForce

subroutine boundaryCondition(VINC, KGLOBAL, FGLOBAL,nnodes)
    implicit none
    real, intent(in) :: VINC(:,:)
    real(8), intent(inout) ::  FGLOBAL(:)
    type(sparse_matrix), intent(inout):: KGLOBAL
    integer,intent(in):: nnodes
    integer :: i, ino, idir, type_,ngl,il
    real(8):: value_ !, nlg
    

    ! print*, "VINC(1,:): ", VINC(1,:)
    ! print*, "VINC(2,:): ", VINC(2,:)
    ! print*, "VINC(3,:): ", VINC(3,:)
    ! print*, "VINC(4,:): ", VINC(4,:)

    Print*, "Applying boundary conditions"
    do i = 1, size(VINC, 1)
        ino=int(VINC(i, 1))
        idir=int(VINC(i, 2))
        value_=VINC(i, 3)
        type_=int(VINC(i, 4))
        
        !_________________!
        ngl=2*(ino-1)+idir
        !_________________!
        
        if (type_==1 .and. value_==0) then
            FGLOBAL(ngl) = 0.d0
            call set_value_to_row(KGLOBAL,ngl,0.0d0)
            call set_value_to_col(KGLOBAL,ngl,0.0d0)
            print*, ngl
            
            call set_value_in_term(KGLOBAL,ngl,ngl,1.0d0)
            
            ! KGLOBAL(:,2*(ino-1)+idir)=0.0d0
            ! KGLOBAL(2*(ino-1)+idir,:)=0.0d0
            ! KGLOBAL(2*(ino-1)+idir,2*(ino-1)+idir)=1.0d0
            ! FGLOBAL(2*(ino-1)+idir)=value
        else if (type_==1 .and. value_/=0) then
            do il=1,2*nnodes
                FGLOBAL(il)=FGLOBAL(il)-value_*get_term_value(KGLOBAL,min(il,ngl),max(il,ngl))
            end do
            FGLOBAL(ngl)=value_
            call set_value_to_row(KGLOBAL,ngl,0.0d0)
            call set_value_to_col(KGLOBAL,ngl,0.0d0)
            call set_value_in_term(KGLOBAL,ngl,ngl,1.0d0)
        ! else
        !     call set_value_in_term(KGLOBAL,ngl,ngl,value_)
        end if
        
    end do
    Print*, "Boundary conditions applied ! "
end subroutine boundaryCondition

end module assembly_module