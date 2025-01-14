module ma87_interface
  
    implicit none
  
    contains

    subroutine ma87solver(ne,n,rows,columns,values,vector)
        
        use hsl_ma87_double
        use hsl_mc68_double
        use hsl_mc69_double
        !input and output variables
        integer(4),intent(in) :: ne,n ! ne - non-null terms in the sparse simetric matrix; n - matrix order (nxn)
        integer(4),dimension(ne),intent(in) :: rows,columns !vectors containing row and column indexes
        real(8),dimension(ne),intent(in) :: values !vector containing matrix terms values
        real(8),dimension(n),intent(inout) :: vector !vector of the equation system and vector where the solution will be returned
        !ma87 computation variables  
  
        type(mc68_control) :: control68
        type(mc68_info)    :: info68
        type(ma87_keep)    :: keep
        type(ma87_control) :: control
        type(ma87_info)    :: info

        integer(4) :: nrhs,lmap,flag
        real(8),dimension(:),allocatable :: val 
        integer(4),dimension(:),allocatable  :: ptr,row,order,map

        ! convert to hsl standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(hsl_matrix_real_sym_psdef,n,n,ne,rows,columns, &
                                ptr,row,flag,val_in=values,val_out=val,lmap=lmap,map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        
        ! call mc68 to find a fill reducing ordering (1=amd)
        allocate(order(n))
        call mc68_order(1,n,ptr,row,order,control68,info68)
        call stop_on_bad_flag("mc68_order",info68%flag)

        ! analyse
        call ma87_analyse(n,ptr,row,order,keep,control,info)
        call stop_on_bad_flag("analyse",info%flag)

        ! factor
        call ma87_factor(n,ptr,row,val,order,keep,control,info)
        call stop_on_bad_flag("factor",info%flag)

        ! solve
        call ma87_solve(vector,order,keep,control,info)
        call stop_on_bad_flag("solve",info%flag)

        ! finalize
        call ma87_finalise(keep, control)
        
        deallocate(val,ptr,row,order,map)
        
        contains
    
        subroutine stop_on_bad_flag(context,flag)
        
            character(len=*), intent(in) :: context
            integer(4),intent(in) :: flag
        
            if(flag.eq.0) return
            write(*,*) "Failure during ", context, " with flag = ", flag
            read(*,*)
            stop
            
        end subroutine stop_on_bad_flag
        
end subroutine ma87solver

end module ma87_interface
