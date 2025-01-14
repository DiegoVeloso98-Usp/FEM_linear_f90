module particle_generator_module

    implicit none
    contains

    !  I want to read a .txt file to populate the PARTICLES, NODESP and KSIETA
    subroutine readParticles(filename, PARTICLES, NODESP, KSIETA)
        implicit none

        ! Input parameters
        character(len=*), intent(in) :: filename
        real, allocatable, intent(out) :: PARTICLES(:,:)
        real(8), allocatable, intent(out) :: NODESP(:,:)
        real(8), allocatable, intent(out) :: KSIETA(:,:)

        ! Local variables
        integer :: nnodes_p, nparticles
        real :: young_p , ni_p , bx_p , by_p , h_p
        integer :: i, unit
        real(8) :: x, y, ksi, eta
        integer :: element

        ! Open the file for reading
        bx_p = 0.0
        by_p = 0.0


        unit = 10
        open(unit, file=filename, status='old', action='read')

        ! Read header information
        read(unit, *) nnodes_p, nparticles, young_p,ni_p, h_p 

        ! Allocate arrays based on the number of nodes and particles
        allocate(NODESP(nnodes_p, 2))
        allocate(KSIETA(nnodes_p, 3))
        allocate(PARTICLES(nparticles, 8))

        ! Read node data
        do i = 1, nnodes_p
            read(unit, *) x, y, element, ksi, eta
            NODESP(i, 1) = x
            NODESP(i, 2) = y
            KSIETA(i, 1) = real(element,kind=8)
            KSIETA(i, 2) = ksi
            KSIETA(i, 3) = eta
        end do

        ! Read particle data
        do i = 1, nparticles
            read(unit, *) PARTICLES(i, 1:3)
            PARTICLES(i, 4) = young_p
            PARTICLES(i, 5) = ni_p
            PARTICLES(i, 6) = h_p
            PARTICLES(i, 7) = bx_p
            PARTICLES(i, 8) = by_p
        end do

        close(unit)
    end subroutine readParticles






    subroutine generateParticles(NODES,ELEMS, npe, gaprox,radius,fv,young,ni,h,bx,by,NODESP,PARTICLES,KSIETA)

        real, dimension(:,:), intent(in) :: ELEMS

        real(8), dimension(:,:), intent(in) :: NODES

        real, intent(in) :: radius,fv,young,ni,bx,by,h
        integer,intent(in) :: gaprox,npe

        real, dimension(:,:),allocatable:: PARTICLES

        real(8), dimension(:,:),allocatable:: NODESP,KSIETA

        integer :: i,j
        real(8) :: xi,yi,xf,yf,dx,dy,xp,yp,area,areap,fraction,r,xc,yc,alfa
        real(8) :: x1,y1,x2,y2,x3,y3
        integer :: index,particle_count,nnodesp,nelems,appended
        real(8) :: ksi, eta

        real :: pi
        pi = 3.1415
        


        ! calculate the maximum and minimum x and y coordinates of the nodes
        xi = minval(NODES(:,1))
        yi = minval(NODES(:,2))
        
        xf = maxval(NODES(:,1))
        yf = maxval(NODES(:,2))

        area= (xf-xi)*(yf-yi)

        xi=xi+radius
        yi=yi+radius

        xf=xf-radius
        yf=yf-radius
        
        dx=xf-xi
        dy=yf-yi

        areap=3*(radius**2)*sqrt(3.0)/4
        fraction=0
        index=0
        
        allocate(PARTICLES(floor(fv*area/areap)+1,8))
        allocate(NODESP(3*(floor(fv*area/areap)+1),2))
        particle_count=0
        do 
            if (fraction>=fv*area) exit
            
            call random_number(r)
            xc=xi+dx*r
            call random_number(r)
            yc=yi+dy*r
            call random_number(r)
            alfa=(2*pi/3)*r
            

            NODESP(index+1,:)=(/xc+radius*cos(alfa),yc+radius*sin(alfa)/)
            NODESP(index+2,:)=(/xc+radius*cos(alfa+2*pi/3),yc+radius*sin(alfa+2*pi/3)/)
            NODESP(index+3,:)=(/xc+radius*cos(alfa+4*pi/3),yc+radius*sin(alfa+4*pi/3)/)
            index=index+3 
            particle_count=particle_count+1
            PARTICLES(particle_count,:)=(/real(index-2),real(index-1),real(index),young,ni,h,bx,by/)


            fraction=fraction+areap

        end do
        
        allocate(KSIETA(size(NODESP,1),3))
        KSIETA=0.0d0
        index=0
        nnodesp=size(NODESP,1)
        nelems=size(ELEMS,1)
        do j=1,nnodesp
            
            xp=NODESP(j,1)
            yp=NODESP(j,2)
            appended=0

            inner_loop: do i=1,nelems
                ! stop
                x1=NODES(int(ELEMS(i,1)),1)
                y1=NODES(int(ELEMS(i,1)),2)

                x2=NODES(int(ELEMS(i,gaprox+1)),1)
                y2=NODES(int(ELEMS(i,gaprox+1)),2)

                x3=NODES(int(ELEMS(i,npe)),1)
                y3=NODES(int(ELEMS(i,npe)),2)

                ksi=-((x3*y1-xp*y1-x1*y3+xp*y3+x1*yp-x3*yp)/(x2*y1-x3*y1-x1*y2+x3*y2+x1*y3-x2*y3))
                eta=-((x2*y1-xp*y1-x1*y2+xp*y2+x1*yp-x2*yp)/(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3))

                if (ksi >= 0.0d0 .and. ksi <= 1.0d0 .and. eta >= 0.0d0 .and. eta <= 1.0d0 &
                    .and. (1.0d0 - ksi - eta) >= 0.0d0 .and. (1.0d0 - ksi - eta) <= 1.0d0) then
                    index=index+1
                    KSIETA(index,:) = (/real(i,kind=8),ksi,eta/)

                    appended=1

                    exit inner_loop 
                end if
                ! if (i==nelems .and. appended==0) then
                !     print*,"ERROR: PARTICLE NOT FOUND"
                !     print*,"xp: ",xp
                !     print*,"yp: ",yp
                ! end if

            end do inner_loop
        end do


    end subroutine generateParticles

end module particle_generator_module
                


