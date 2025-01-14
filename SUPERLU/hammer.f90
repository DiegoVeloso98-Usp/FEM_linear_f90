module hammer_module
    implicit none
    contains

    subroutine hammer(nph, hammer_points)
        ! Declare outputs
        integer, intent(out) :: nph
        real(8), intent(out) :: hammer_points(3, 12)

        ! Assign the number of points
        nph = 12

        ! Assign values to the hammer array in row-major order (like Python)
        hammer_points(1,:) = [0.501426509658179, 0.249286745170910, 0.249286745170910, &
            0.873821971016996, 0.063089014491502, 0.063089014491502, &
            0.053145049844816, 0.310352451033785, 0.636502499121399, &
            0.310352451033785, 0.636502499121399, 0.053145049844816]

        hammer_points(2,:) = [0.249286745170910, 0.249286745170910, 0.501426509658179, &
            0.063089014491502, 0.063089014491502, 0.873821971016996, &
            0.310352451033785, 0.636502499121399, 0.053145049844816, &
            0.053145049844816, 0.310352451033785, 0.636502499121399]

        hammer_points(3,:) = [0.116786275726379/2.0, 0.116786275726379/2.0, 0.116786275726379/2.0, &
            0.050844906370207/2.0, 0.050844906370207/2.0, 0.050844906370207/2.0, &
            0.082851075618374/2.0, 0.082851075618374/2.0, 0.082851075618374/2.0, &
            0.082851075618374/2.0, 0.082851075618374/2.0, 0.082851075618374/2.0]
    end subroutine hammer


end module hammer_module